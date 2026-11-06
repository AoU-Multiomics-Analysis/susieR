#!/usr/bin/env python3
"""Check WDL types and simulate cloud localization before Bash rendering.

These tests run the real task command with a recording Rscript replacement.
They do not run Docker, fit models, or submit a Terra workflow.
"""
import json
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile
import unittest

import WDL

ROOT = pathlib.Path(sys.argv.pop(1)).resolve() if len(sys.argv) > 1 else pathlib.Path(__file__).resolve().parents[1]
WRITERS = {'write_lines', 'write_tsv', 'write_map', 'write_json'}
GENE_URI_INPUTS = {
    'gene_rds_uris': 'GeneRds',
    'gene_parquet_uris': 'GeneParquets',
    'gene_tsv_uris': 'GeneTsvs',
    'gene_log_uris': 'GeneLogs',
    'gene_provenance_uris': 'GeneProvenance',
}
GENE_URI_FLAGS = {'--' + name.replace('_', '-') + '-file': name for name in GENE_URI_INPUTS}
LIST_FLAGS = ('--summary-files', '--credible-set-files', '--index-files', '--diagnostic-files', '--chromosomes-file') + tuple(GENE_URI_FLAGS)
FILE_NAMES = {'phenotype_bed', 'qtl_results', 'gene_map', 'worklist', 'covariates', 'sample_list',
              'ancestry', 'variant_list', 'additional_genotypes', 'expected_outcomes'}
OPTIONAL_FILES = {'gene_map', 'ancestry', 'variant_list', 'additional_genotypes'}


def walk(node):
    yield node
    for child in node.children:
        yield from walk(child)


def workflow_writers(workflow):
    return [node for node in walk(workflow)
            if isinstance(node, WDL.Expr.Apply) and node.function_name in WRITERS]


def has_file_type(typ):
    if isinstance(typ, WDL.Type.File):
        return True
    if isinstance(typ, WDL.Type.Array):
        return has_file_type(typ.item_type)
    if isinstance(typ, WDL.Type.StructInstance):
        return any(has_file_type(t) for t in typ.members.values())
    return False


class CommandStdLib(WDL.StdLib.Base):
    # Evaluate task-command writers after the test has localized File values.
    def _virtualize_filename(self, filename):
        return filename


class PgenWdlTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path = ROOT / 'workflows/SusieR2.wdl'
        if not cls.path.is_file():
            raise AssertionError('Missing SuSiE 2.0 WDL')
        cls.doc = WDL.load(str(cls.path))
        cls.tasks = {task.name: task for task in cls.doc.tasks}

    def test_workflow_types_and_defaults(self):
        self.assertEqual(self.doc.effective_wdl_version, '1.0')
        inputs = {decl.name: decl for decl in self.doc.workflow.inputs}
        chrom_type = inputs['ChromosomePgens'].type
        self.assertIsInstance(chrom_type, WDL.Type.Array)
        self.assertEqual(str(chrom_type.item_type), 'ChromosomePgen')
        self.assertEqual({k: str(v) for k, v in chrom_type.item_type.members.items()},
                         {'chromosome': 'String', 'pgen': 'File', 'pvar': 'File', 'psam': 'File'})
        self.assertIsNone(inputs['LegacyMissingCode'].expr)
        self.assertEqual(inputs['BatchSize'].expr.eval(WDL.Env.Bindings(), WDL.StdLib.Base('1.0')).value, 50)
        self.assertIsInstance(inputs['MemoryDebug'].type, WDL.Type.Boolean)
        self.assertFalse(inputs['MemoryDebug'].expr.eval(WDL.Env.Bindings(), WDL.StdLib.Base('1.0')).value)
        batch_call = next(node for node in walk(self.doc.workflow)
                          if isinstance(node, WDL.Tree.Call) and node.name == 'RunPgenSusieBatch')
        for enabled in (False, True):
            env = WDL.Env.Bindings().bind('MemoryDebug', WDL.Value.Boolean(enabled))
            self.assertEqual(batch_call.inputs['memory_debug'].eval(env, WDL.StdLib.Base('1.0')).value,
                             enabled)
        for name in ('GeneMap', 'Ancestry', 'VariantList', 'AdditionalGenotypes'):
            self.assertIsInstance(inputs[name].type, WDL.Type.File)
            self.assertTrue(inputs[name].type.optional)
        for task in self.doc.tasks:
            # No pre-command declarations can turn File inputs into strings.
            self.assertFalse(task.postinputs)
        planner_files = {d.name for d in self.tasks['PlanPgenSusie'].inputs if has_file_type(d.type)}
        self.assertEqual(planner_files, {'phenotype_bed', 'qtl_results', 'gene_map'})
        for name in ('summary_files', 'credible_set_files', 'index_files', 'diagnostic_files'):
            decl = next(d for d in self.tasks['AggregatePgenSusie'].inputs if d.name == name)
            self.assertIsInstance(decl.type.item_type, WDL.Type.File)
        # Output paths remain typed, including per-gene artifacts and diagnostics.
        outputs = {d.name: str(d.type) for d in self.doc.workflow.outputs}
        for name in ('GeneRds', 'GeneParquets', 'GeneTsvs', 'GeneLogs', 'GeneProvenance',
                     'PlinkLogs', 'ChromosomeGeneManifests'):
            self.assertEqual(outputs[name], 'Array[File]')

    def test_gene_manifest_metadata_keeps_completed_output_uris(self):
        task = self.tasks['AggregatePgenSusie']
        declarations = {decl.name: decl for decl in task.inputs}
        call = next(node for node in walk(self.doc.workflow)
                    if isinstance(node, WDL.Tree.Call) and node.name == 'AggregatePgenSusie')
        outputs = {decl.name: decl for decl in task.outputs}
        self.assertEqual(str(outputs['ChromosomeGeneManifests'].type), 'Array[File]')
        for name, output_name in GENE_URI_INPUTS.items():
            with self.subTest(metadata=name):
                declaration = declarations[name]
                self.assertEqual(str(declaration.type), 'Array[String]')
                uris = ['gs://bucket/outputs/first ' + name,
                        "gs://bucket/outputs/second 'quoted' $(touch PWNED) `touch PWNED` " + name]
                source_type = WDL.Type.Array(WDL.Type.Array(WDL.Type.File()))
                source = WDL.Value.from_json(source_type, [[uris[0]], [uris[1]]])
                env = WDL.Env.Bindings().bind('RunPgenSusieBatch.' + output_name, source)
                # File outputs become text only for output-manifest metadata.
                # They are not files that this task opens or localizes.
                evaluated = call.inputs[name].eval(env, WDL.StdLib.Base('1.0'))
                self.assertIsInstance(evaluated.type.item_type, WDL.Type.File)
                metadata = evaluated.coerce(declaration.type)
                self.assertEqual(metadata.json, uris)

                def reject_localization(file):
                    self.fail('Output URI metadata must not be localized: ' + file.value)

                rewritten = WDL.Value.rewrite_paths(metadata, reject_localization)
                self.assertEqual(rewritten.json, uris)

    def test_example_inputs_match_wdl_types(self):
        inputs = json.loads((ROOT / 'examples/inputs/SusieR2.inputs.json').read_text())
        WDL.values_from_json(inputs, self.doc.workflow.available_inputs,
                             self.doc.workflow.required_inputs, namespace='SusieR2')

    def test_batch_indices_select_matching_typed_chromosome_files(self):
        # Input order differs from genomic order; two batches reuse chromosome 1.
        typ = next(d.type for d in self.doc.workflow.inputs if d.name == 'ChromosomePgens')
        chromosomes = [{'chromosome': chrom, **{ext: 'gs://bucket/chr' + chrom + '.' + ext
                                               for ext in ('pgen', 'pvar', 'psam')}}
                       for chrom in ('2', '1')]
        files = WDL.Value.from_json(typ, chromosomes)
        env = WDL.Env.Bindings().bind('ChromosomePgens', files)
        env = env.bind('PlanPgenSusie.ChromosomeIndices',
                       WDL.Value.from_json(WDL.Type.Array(WDL.Type.Int()), [1, 1, 0]))
        call = next(node for node in walk(self.doc.workflow)
                    if isinstance(node, WDL.Tree.Call) and node.name == 'RunPgenSusieBatch')
        index = next(node for node in walk(self.doc.workflow)
                     if isinstance(node, WDL.Tree.Decl) and node.name == 'chromosome_index')
        self.assertIsInstance(index.type, WDL.Type.Int)
        for i, chromosome in enumerate(('1', '1', '2')):
            iteration_env = env.bind('i', WDL.Value.Int(i))
            chromosome_index = index.expr.eval(iteration_env, WDL.StdLib.Base('1.0'))
            iteration_env = iteration_env.bind(index.name, chromosome_index)
            selected = call.inputs['pgen_files'].eval(iteration_env, WDL.StdLib.Base('1.0'))
            self.assertEqual(selected.value['chromosome'].value, chromosome)
            for ext in ('pgen', 'pvar', 'psam'):
                self.assertIsInstance(selected.value[ext], WDL.Value.File)
                self.assertEqual(selected.value[ext].value, 'gs://bucket/chr' + chromosome + '.' + ext)

    def test_no_workflow_file_writers(self):
        self.assertEqual(workflow_writers(self.doc.workflow), [])

    def test_command_file_lists_preserve_generated_file_types(self):
        # Cromwell maps a generated File to its container path at the final
        # command placeholder. A String wrapper would hide that File mapping.
        writers = []
        for task in self.doc.tasks:
            for part in task.command.parts:
                if not isinstance(part, WDL.Expr.Placeholder):
                    continue
                for node in walk(part.expr):
                    if isinstance(node, WDL.Expr.Apply) and node.function_name == 'write_lines':
                        writers.append((task.name, str(node.arguments[0])))
                        self.assertIs(part.expr, node)
                        self.assertIsInstance(part.expr.type, WDL.Type.File)
                        name = str(node.arguments[0])
                        metadata = name in GENE_URI_INPUTS
                        expected_item_type = WDL.Type.String if task.name == 'PlanPgenSusie' or metadata else WDL.Type.File
                        self.assertIsInstance(node.arguments[0].type.item_type, expected_item_type)
        self.assertCountEqual(writers, [
            ('PlanPgenSusie', 'chromosomes'),
            ('AggregatePgenSusie', 'summary_files'),
            ('AggregatePgenSusie', 'credible_set_files'),
            ('AggregatePgenSusie', 'index_files'),
            ('AggregatePgenSusie', 'diagnostic_files'),
            *[('AggregatePgenSusie', name) for name in GENE_URI_INPUTS],
        ])

    def test_static_guard_rejects_writers_in_all_workflow_scopes(self):
        # Test the guard itself: nested workflow expressions must not escape it.
        for writer, expression in (
            ('write_lines', 'write_lines(["x"])'),
            ('write_tsv', 'write_tsv([["x"]])'),
            ('write_map', 'write_map({"x":"y"})'),
            ('write_json', 'write_json(["x"])'),
        ):
            for scope in ('declaration', 'scatter', 'conditional', 'call', 'output'):
                source = 'version 1.0\ntask T { input { File f } command { true } }\nworkflow W {\n'
                if scope == 'declaration':
                    source += ' File bad = ' + expression
                elif scope == 'scatter':
                    source += ' scatter (i in [1]) { File bad = ' + expression + ' }'
                elif scope == 'conditional':
                    source += ' if (true) { File bad = ' + expression + ' }'
                elif scope == 'call':
                    source += ' call T { input: f = ' + expression + ' }'
                else:
                    source += ' output { File bad = ' + expression + ' }'
                source += '\n}'
                with self.subTest(writer=writer, scope=scope):
                    parsed = WDL.parse_document(source)
                    self.assertEqual([n.function_name for n in workflow_writers(parsed.workflow)], [writer])

    def inputs(self, task, optional):
        raw = {}
        for decl in task.inputs:
            name, typ = decl.name, decl.type
            if name == 'pgen_files':
                value = {'chromosome': '2', 'pgen': 'gs://bucket/a.pgen',
                         'pvar': 'gs://bucket/a.pvar', 'psam': 'gs://bucket/a.psam'}
            elif name in OPTIONAL_FILES and not optional:
                value = None
            elif isinstance(typ, WDL.Type.File):
                value = 'gs://bucket/' + name
            elif isinstance(typ, WDL.Type.Array) and isinstance(typ.item_type, WDL.Type.File):
                value = ['gs://bucket/' + name + '/one', 'gs://bucket/' + name + '/two']
            elif name in GENE_URI_INPUTS:
                value = ['gs://bucket/outputs/' + name + '/one',
                         "gs://bucket/outputs/" + name + "/two 'quoted' $(touch PWNED) `touch PWNED`"]
            elif name == 'chromosomes':
                # Invalid chromosome IDs test shell quoting; the R planner validates IDs.
                value = ['2', "1 'one' 'two' $(touch PWNED) `touch PWNED`"]
            elif isinstance(typ, WDL.Type.Boolean):
                value = optional
            elif isinstance(typ, WDL.Type.Int):
                value = {'batch_size': 100, 'cis_distance': 1000000, 'base_seed': 1,
                         'plink_memory_mb': 1024}.get(name, 2)
            elif isinstance(typ, WDL.Type.Float):
                value = 0.0
            elif name == 'missing_code':
                value = '-1' if optional else 'NA'
            else:
                value = "qval 'one' 'two' $(touch PWNED) `touch PWNED`"
            raw[name] = WDL.Value.from_json(typ, value)
        return raw

    def execute(self, task, optional=True, broken=None):
        with tempfile.TemporaryDirectory(prefix='susier-wdl-') as temporary:
            cwd = pathlib.Path(temporary)
            raw = self.inputs(task, optional)
            mapping = {}

            def localize(file):
                if file.value not in mapping:
                    dest = cwd / ("input 'one' 'two' $(touch PWNED) `touch PWNED` " + str(len(mapping)))
                    dest.write_text('fixture\n')
                    mapping[file.value] = str(dest)
                return mapping[file.value]

            localized = {name: WDL.Value.rewrite_paths(value, localize) for name, value in raw.items()}
            for name in GENE_URI_INPUTS.keys() & raw.keys():
                self.assertEqual(localized[name].json, raw[name].json)
                self.assertTrue(all(uri not in mapping for uri in raw[name].json))
            if broken:
                # Mutate a File only after constructing the proper localized environment.
                name = broken[0]
                if broken[1] == 'uri':
                    localized[name] = raw[name]
                else:
                    pathlib.Path(localized[name].value).unlink()
            env = WDL.Env.Bindings()
            for name, value in localized.items():
                env = env.bind(name, value)
            command = task.command.eval(env, CommandStdLib('1.0', write_dir=str(cwd))).value
            (cwd / 'command.sh').write_text(command)
            bindir = cwd / 'bin'
            bindir.mkdir()
            stub = bindir / 'Rscript'
            stub.write_text('#!' + sys.executable + '\n' + '''import json, pathlib, sys
args = sys.argv[1:]
lists = {}
for i, arg in enumerate(args):
    if arg in ''' + repr(LIST_FLAGS) + ''':
        lists[arg] = pathlib.Path(args[i+1]).read_text().splitlines()
pathlib.Path('capture.json').write_text(json.dumps({'args': args, 'lists': lists}))
''')
            stub.chmod(0o755)
            result = subprocess.run(['bash', 'command.sh'], cwd=cwd,
                                    env={**os.environ, 'PATH': str(bindir) + os.pathsep + os.environ['PATH']},
                                    text=True, capture_output=True)
            self.assertFalse((cwd / 'PWNED').exists(), command)
            captured = cwd / 'capture.json'
            if broken:
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(captured.exists(), 'Computation must not start with an unresolved input')
                self.assertIn('ERROR', result.stderr)
                if broken[1] == 'uri':
                    self.assertIn('localization', result.stderr)
                return
            self.assertEqual(result.returncode, 0, result.stderr + '\n' + command)
            payload = json.loads(captured.read_text())
            args = payload['args']
            self.assertNotIn('gs://', ' '.join(args))
            if task.name == 'PlanPgenSusie':
                self.assertEqual(args[0], '/opt/r/scripts/plan_pgen_susie.R')
                self.assertEqual(payload['lists']['--chromosomes-file'], localized['chromosomes'].json)
                self.assertEqual(args[args.index('--batch-size') + 1], '100')
            elif task.name == 'RunPgenSusieBatch':
                self.assertEqual(args[0], '/opt/r/scripts/run_pgen_susie_batch.R')
                self.assertEqual('--memory-debug' in args, optional)
                self.assertIn('--missing-code=' + localized['missing_code'].value, args)
                self.assertNotIn('--missing-code', args)
                for name, value in localized['pgen_files'].value.items():
                    if name != 'chromosome':
                        self.assertEqual(args[args.index('--' + name) + 1], value.value)
            else:
                self.assertEqual(args[0], '/opt/r/scripts/aggregate_pgen_susie.R')
                for flag, paths in payload['lists'].items():
                    if flag in GENE_URI_FLAGS:
                        self.assertEqual(paths, raw[GENE_URI_FLAGS[flag]].json)
                    else:
                        self.assertEqual(paths, localized[flag[2:].replace('-', '_')].json)
                        for path in paths:
                            self.assertTrue(pathlib.Path(path).is_file())
                for flag, name in GENE_URI_FLAGS.items():
                    self.assertIn(flag, args)
                    self.assertEqual(payload['lists'][flag], raw[name].json)
            for name, value in localized.items():
                if name in FILE_NAMES:
                    flag = '--' + name.replace('_', '-')
                    if name in OPTIONAL_FILES and not optional:
                        self.assertNotIn(flag, args)
                    else:
                        self.assertEqual(args[args.index(flag) + 1], value.value)
            if task.name != 'AggregatePgenSusie':
                flag = '--top-phenotype-pvalue-column'
                self.assertEqual(args[args.index(flag) + 1], localized['top_phenotype_pvalue_column'].value)
                self.assertEqual('--select-top-phenotype-per-cluster' in args, optional)
            return payload

    def test_optparse_accepts_rendered_negative_missing_code(self):
        rscript = shutil.which('Rscript')
        if not rscript:
            self.skipTest('Rscript is not installed; command-rendering assertions still run')
        payload = self.execute(self.tasks['RunPgenSusieBatch'], optional=True)
        argument = next(arg for arg in payload['args'] if arg.startswith('--missing-code='))
        script = 'library_path <- Sys.getenv("BENCH_LIBRARY", unset = "")\nif (nzchar(library_path)) .libPaths(c(library_path, .libPaths()))\nif (!requireNamespace("optparse", quietly = TRUE)) quit(status = 77L)\nparser <- optparse::OptionParser(option_list = list(\n  optparse::make_option("--missing-code", type = "character")))\nparsed <- optparse::parse_args(parser, args = commandArgs(trailingOnly = TRUE),\n                              print_help_and_exit = FALSE)\nstopifnot(identical(parsed[["missing-code"]], "-1"))\n'
        result = subprocess.run([rscript, '-', argument], input=script, text=True, capture_output=True)
        if result.returncode == 77:
            self.skipTest('R optparse is not installed; command-rendering assertions still run')
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_all_task_commands_after_cloud_localization(self):
        for task in self.doc.tasks:
            for optional in (False, True):
                with self.subTest(task=task.name, optional=optional):
                    self.execute(task, optional)

    def test_all_tasks_reject_unresolved_cloud_files_and_missing_inputs(self):
        for task in self.doc.tasks:
            name = 'expected_outcomes' if task.name == 'AggregatePgenSusie' else 'phenotype_bed'
            for failure in ('uri', 'missing'):
                with self.subTest(task=task.name, failure=failure):
                    self.execute(task, broken=(name, failure))

    def test_aggregate_file_lists_reject_unresolved_cloud_paths(self):
        for name in ('summary_files', 'credible_set_files', 'index_files', 'diagnostic_files'):
            with self.subTest(file_list=name):
                self.execute(self.tasks['AggregatePgenSusie'], broken=(name, 'uri'))


if __name__ == '__main__':
    unittest.main()
