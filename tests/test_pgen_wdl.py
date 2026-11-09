#!/usr/bin/env python3
"""Check WDL types, output collection, and simulated cloud file transfer.

These tests run the real task command with a recording Rscript replacement.
The transfer test copies selected output files, removes the source directory,
and localizes those files for the next task. It does not run Cromwell upload
code, Docker, model fits, or a Terra workflow.
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
PGEN_METADATA_FLAGS = {'--' + name.replace('_', '-') + '-file': name for name in
                       ('pgen_chromosomes', 'pgen_uris', 'pvar_uris', 'psam_uris')}
METADATA_INPUTS = set(GENE_URI_INPUTS) | set(PGEN_METADATA_FLAGS.values())
LIST_FLAGS = ('--summary-files', '--credible-set-files', '--index-files', '--diagnostic-files') + tuple(GENE_URI_FLAGS) + tuple(PGEN_METADATA_FLAGS)
FILE_NAMES = {'pgen_manifest', 'phenotype_bed', 'qtl_results', 'gene_map', 'worklist', 'covariates', 'sample_list',
              'ancestry', 'variant_list', 'additional_genotypes', 'genotype_qc_report', 'expected_outcomes', 'pgen', 'pvar', 'psam'}
OPTIONAL_FILES = {'gene_map', 'sample_list', 'ancestry', 'variant_list', 'additional_genotypes', 'genotype_qc_report'}


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

    def _devirtualize_filename(self, filename):
        return filename


class TaskOutputStdLib(CommandStdLib):
    """Evaluate real task outputs with Bash glob order in a fixture directory."""

    def __init__(self, directory):
        super().__init__('1.0')
        self.directory = directory
        self.collected_files = set()
        self.glob = WDL.StdLib.StaticFunction(
            'glob', [WDL.Type.String()], WDL.Type.Array(WDL.Type.File()), self.collect_glob)

    def _devirtualize_filename(self, filename):
        return str(self.directory / filename)

    def collect_glob(self, pattern):
        result = subprocess.run(
            ['bash', '-c', 'shopt -s nullglob; paths=( $1 ); '
             'if (( ${#paths[@]} )); then printf "%s\\0" "${paths[@]}"; fi',
             '--', pattern.value], cwd=self.directory, check=True, capture_output=True)
        paths = [self.directory / name for name in result.stdout.decode().split('\0') if name]
        self.collected_files.update(paths)
        return WDL.Value.Array(WDL.Type.File(), [WDL.Value.File(str(path)) for path in paths])


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
        self.assertIsInstance(inputs['PgenManifest'].type, WDL.Type.File)
        self.assertIsNone(inputs['PgenManifest'].expr)
        self.assertNotIn('ChromosomePgens', inputs)
        required = {name for name, decl in inputs.items()
                    if decl.expr is None and not decl.type.optional}
        self.assertEqual(required, {'PgenManifest', 'PhenotypeBed', 'QtlResults',
                                    'Covariates', 'LegacyMissingCode'})
        chrom_type = next(node.type for node in walk(self.doc.workflow)
                          if isinstance(node, WDL.Tree.Decl) and node.name == 'ChromosomePgens')
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
        for name in ('SampleList', 'GeneMap', 'Ancestry', 'VariantList', 'AdditionalGenotypes'):
            self.assertIsInstance(inputs[name].type, WDL.Type.File)
            self.assertTrue(inputs[name].type.optional)
        for task in self.doc.tasks:
            # No pre-command declarations can turn File inputs into strings.
            self.assertFalse(task.postinputs)
        planner_files = {d.name for d in self.tasks['PlanPgenSusie'].inputs if has_file_type(d.type)}
        self.assertEqual(planner_files, {'pgen_manifest', 'phenotype_bed', 'qtl_results', 'gene_map'})
        planner_outputs = {decl.name: decl for decl in self.tasks['PlanPgenSusie'].outputs}
        self.assertIsInstance(planner_outputs['NormalizedManifest'].type, WDL.Type.File)
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

    def test_planner_collects_dynamic_files_with_glob(self):
        # A list of filenames does not cause the GCP backend to upload them.
        outputs = {decl.name: decl for decl in self.tasks['PlanPgenSusie'].outputs
                   if isinstance(decl.type, WDL.Type.Array)
                   and isinstance(decl.type.item_type, WDL.Type.File)}
        self.assertEqual(set(outputs), {'Worklists', 'PhenotypeBeds'})
        for name, declaration in outputs.items():
            with self.subTest(output=name):
                self.assertIsInstance(declaration.expr, WDL.Expr.Apply)
                self.assertEqual(declaration.expr.function_name, 'glob',
                                 'Dynamic task files need explicit glob collection before upload')

    def test_example_inputs_match_wdl_types(self):
        inputs = json.loads((ROOT / 'examples/inputs/SusieR2.inputs.json').read_text())
        WDL.values_from_json(inputs, self.doc.workflow.available_inputs,
                             self.doc.workflow.required_inputs, namespace='SusieR2')
        self.assertNotIn('SusieR2.ChromosomePgens', inputs)
        explicit = {**inputs, 'SusieR2.SampleList': 'gs://bucket/optional-samples.tsv'}
        WDL.values_from_json(explicit, self.doc.workflow.available_inputs,
                             self.doc.workflow.required_inputs, namespace='SusieR2')
        absent = {name: value for name, value in inputs.items() if name != 'SusieR2.SampleList'}
        WDL.values_from_json(absent, self.doc.workflow.available_inputs,
                             self.doc.workflow.required_inputs, namespace='SusieR2')

    def test_shared_qc_defaults_and_threshold_forwarding(self):
        inputs = {decl.name: decl for decl in self.doc.workflow.inputs}
        self.assertIn('MAC', inputs)
        self.assertIsInstance(inputs['MAC'].type, WDL.Type.Int)
        self.assertEqual(inputs['MAC'].expr.eval(WDL.Env.Bindings(), WDL.StdLib.Base('1.0')).value, 0)
        self.assertIn('MinCallRate', inputs)
        self.assertEqual(inputs['MinCallRate'].expr.eval(WDL.Env.Bindings(), WDL.StdLib.Base('1.0')).value, 0.95)
        qc_call = next(node for node in walk(self.doc.workflow)
                       if isinstance(node, WDL.Tree.Call) and node.name == 'FilterPgen')
        call = next(node for node in walk(self.doc.workflow)
                    if isinstance(node, WDL.Tree.Call) and node.name == 'RunPgenSusieBatch')
        for threshold in (0, 5):
            with self.subTest(mac=threshold):
                env = WDL.Env.Bindings().bind('MAC', WDL.Value.Int(threshold))
                forwarded = qc_call.inputs['mac'].eval(env, WDL.StdLib.Base('1.0'))
                self.assertEqual(forwarded.value, threshold)
                self.assertEqual(call.inputs['mac'].eval(env, WDL.StdLib.Base('1.0')).value, 0,
                                 'Batches must not repeat shared MAC filtering')
                payload = self.execute(self.tasks['FilterPgen'], overrides={'mac': forwarded})
                self.assertIn('--mac=' + str(threshold), payload['args'])

    def manifest_files(self):
        rows = next(node for node in walk(self.doc.workflow)
                    if isinstance(node, WDL.Tree.Decl) and node.name == 'pgen_rows')
        bundle = next(node for node in walk(self.doc.workflow)
                      if isinstance(node, WDL.Tree.Decl) and node.name == 'chromosome_pgen')
        array_type = next(node.type for node in walk(self.doc.workflow)
                          if isinstance(node, WDL.Tree.Decl) and node.name == 'ChromosomePgens')
        self.assertIsInstance(rows.expr, WDL.Expr.Apply)
        self.assertEqual(rows.expr.function_name, 'read_tsv')
        values = [[chrom, *['gs://bucket/' + chrom + " 'quoted' \"double\" $(touch PWNED) `touch PWNED`." + ext
                            for ext in ('pgen', 'pvar', 'psam')]] for chrom in ('2', '1')]
        with tempfile.TemporaryDirectory(prefix='susier-manifest-rows-') as temporary:
            normalized = pathlib.Path(temporary) / 'chromosome_pgens.tsv'
            normalized.write_text(''.join('\t'.join(row) + '\n' for row in values))
            env = WDL.Env.Bindings().bind('PlanPgenSusie.UniquePgenManifest',
                                          WDL.Value.File(str(normalized)))
            parsed = rows.expr.eval(env, CommandStdLib('1.0'))
        self.assertEqual(parsed.json, values)
        files = []
        for row in parsed.value:
            env = WDL.Env.Bindings().bind('pgen_row', row)
            value = bundle.expr.eval(env, WDL.StdLib.Base('1.0')).coerce(bundle.type)
            self.assertEqual(value.value['chromosome'].value, row.value[0].value)
            for field, original in zip(('pgen', 'pvar', 'psam'), row.value[1:]):
                self.assertIsInstance(value.value[field], WDL.Value.File)
                self.assertEqual(value.value[field].value, original.value)
            files.append(value)
        return WDL.Value.Array(array_type.item_type, files)

    def test_manifest_paths_become_files_before_batch_localization(self):
        files = self.manifest_files()
        inputs = {decl.name: decl for decl in self.tasks['RunPgenSusieBatch'].inputs}
        for name in ('pgen', 'pvar', 'psam'):
            self.assertIsInstance(inputs[name].type, WDL.Type.File)
        for file_bundle in files.value:
            self.execute(self.tasks['RunPgenSusieBatch'], overrides={
                name: file_bundle.value[name] for name in ('pgen', 'pvar', 'psam')})

    def test_shared_qc_reuses_files_and_reports_across_batches(self):
        sources = self.manifest_files()
        source_scatter = next(node for node in walk(self.doc.workflow)
                              if isinstance(node, WDL.Tree.Scatter) and node.variable == 'pgen_row')
        condition = next(node for node in walk(source_scatter) if isinstance(node, WDL.Tree.Conditional))
        qc_calls = [node for node in walk(source_scatter)
                    if isinstance(node, WDL.Tree.Call) and node.name == 'FilterPgen']
        self.assertEqual(len(qc_calls), 1)
        self.assertFalse(any(isinstance(node, WDL.Tree.Call) and node.name == 'FilterPgen'
                             for node in walk(self.batch_scatter())),
                         'QC must not run inside the batch scatter')
        file_declarations = {node.name: node for node in walk(source_scatter)
                             if isinstance(node, WDL.Tree.Decl) and node.name.startswith('prepared_')}
        batch_call = next(node for node in walk(self.doc.workflow)
                          if isinstance(node, WDL.Tree.Call) and node.name == 'RunPgenSusieBatch')
        index = next(node for node in walk(self.batch_scatter())
                     if isinstance(node, WDL.Tree.Decl) and node.name == 'chromosome_index')
        for mac, rate, enabled in ((0, 0.0, False), (0, .95, True), (5, 0.0, True)):
            with self.subTest(mac=mac, min_call_rate=rate):
                switches = WDL.Env.Bindings().bind('MAC', WDL.Value.Int(mac))
                switches = switches.bind('MinCallRate', WDL.Value.Float(rate))
                self.assertEqual(condition.expr.eval(switches, WDL.StdLib.Base('1.0')).value, enabled)
                prepared, reports = [], []
                for i, source in enumerate(sources.value):
                    env = switches.bind('chromosome_pgen', source)
                    for role, output in (('pgen', 'Pgen'), ('pvar', 'Pvar'), ('psam', 'Psam')):
                        uri = 'gs://qc/filtered_%d.%s' % (i, role)
                        env = env.bind('FilterPgen.' + output,
                                       WDL.Value.File(uri) if enabled else WDL.Value.Null())
                        declaration = file_declarations['prepared_' + role]
                        selected = declaration.expr.eval(env, WDL.StdLib.Base('1.0')).coerce(declaration.type)
                        self.assertIsInstance(selected, WDL.Value.File)
                        self.assertEqual(selected.value, uri if enabled else source.value[role].value)
                        env = env.bind(declaration.name, selected)
                    declaration = file_declarations['prepared_files']
                    prepared.append(declaration.expr.eval(env, WDL.StdLib.Base('1.0')).coerce(declaration.type))
                    reports.append('gs://qc/report_%d.json' % i if enabled else None)
                env = WDL.Env.Bindings().bind('ChromosomePgens', WDL.Value.Array(sources.type.item_type, prepared))
                env = env.bind('FilterPgen.Report', WDL.Value.from_json(
                    WDL.Type.Array(WDL.Type.File(optional=True)), reports))
                for name, values, typ in (
                    ('Worklists', ['gs://plan/a.tsv', 'gs://plan/b.tsv', 'gs://plan/c.tsv'], WDL.Type.File()),
                    ('PhenotypeBeds', ['gs://plan/a.bed', 'gs://plan/b.bed', 'gs://plan/c.bed'], WDL.Type.File()),
                    ('BatchPgenIndices', [0, 0, 1], WDL.Type.Int()),
                    ('Chromosomes', ['1', '2', '3'], WDL.Type.String()),
                ):
                    env = env.bind('PlanPgenSusie.' + name, WDL.Value.from_json(WDL.Type.Array(typ), values))
                scatter = self.batch_scatter()
                for batch, expected in zip(scatter.expr.eval(env, WDL.StdLib.Base('1.0')).value, (0, 0, 1)):
                    iteration = env.bind(scatter.variable, batch)
                    iteration = iteration.bind(index.name, index.expr.eval(iteration, WDL.StdLib.Base('1.0')))
                    for role in ('pgen', 'pvar', 'psam'):
                        actual = batch_call.inputs[role].eval(iteration, WDL.StdLib.Base('1.0'))
                        self.assertIsInstance(actual, WDL.Value.File)
                        self.assertEqual(actual.value, prepared[expected].value[role].value)
                    actual_report = batch_call.inputs['genotype_qc_report'].eval(iteration, WDL.StdLib.Base('1.0'))
                    self.assertEqual(actual_report.json, reports[expected])
                    self.assertEqual(batch_call.inputs['mac'].eval(iteration, WDL.StdLib.Base('1.0')).value, 0)

    def test_batch_indices_select_matching_typed_chromosome_files(self):
        # Input order differs from genomic order; two batches reuse chromosome 1.
        files = self.manifest_files()
        env = WDL.Env.Bindings().bind('ChromosomePgens', files)
        env = env.bind('PlanPgenSusie.BatchPgenIndices',
                       WDL.Value.from_json(WDL.Type.Array(WDL.Type.Int()), [1, 1, 0]))
        env = env.bind('PlanPgenSusie.Chromosomes',
                       WDL.Value.from_json(WDL.Type.Array(WDL.Type.String()), ['1', '1', '2']))
        for output, extension in (('Worklists', 'tsv'), ('PhenotypeBeds', 'bed')):
            env = env.bind('PlanPgenSusie.' + output, WDL.Value.from_json(
                WDL.Type.Array(WDL.Type.File()),
                ['gs://bucket/plan/batch_%06d.%s' % (i, extension) for i in (1, 2, 3)]))
        call = next(node for node in walk(self.doc.workflow)
                    if isinstance(node, WDL.Tree.Call) and node.name == 'RunPgenSusieBatch')
        scatter = self.batch_scatter()
        batches = scatter.expr.eval(env, WDL.StdLib.Base('1.0'))
        self.assertEqual(len(batches.value), 3)
        index = next(node for node in walk(self.doc.workflow)
                     if isinstance(node, WDL.Tree.Decl) and node.name == 'chromosome_index')
        self.assertIsInstance(index.type, WDL.Type.Int)
        for i, (batch, chromosome) in enumerate(zip(batches.value, ('1', '1', '2')), 1):
            iteration_env = env.bind(scatter.variable, batch)
            chromosome_index = index.expr.eval(iteration_env, WDL.StdLib.Base('1.0'))
            iteration_env = iteration_env.bind(index.name, chromosome_index)
            for ext in ('pgen', 'pvar', 'psam'):
                selected = call.inputs[ext].eval(iteration_env, WDL.StdLib.Base('1.0'))
                self.assertIsInstance(selected, WDL.Value.File)
                self.assertEqual(selected.value,
                                 files.value[(1, 1, 0)[i - 1]].value[ext].value)
                self.assertIn('gs://bucket/' + chromosome + ' ', selected.value)
            for name, extension in (('worklist', 'tsv'), ('phenotype_bed', 'bed')):
                selected_file = call.inputs[name].eval(iteration_env, WDL.StdLib.Base('1.0'))
                self.assertIsInstance(selected_file, WDL.Value.File)
                self.assertEqual(selected_file.value, 'gs://bucket/plan/batch_%06d.%s' % (i, extension))

    def batch_scatter(self):
        return next(node for node in walk(self.doc.workflow)
                    if isinstance(node, WDL.Tree.Scatter)
                    and any(isinstance(child, WDL.Tree.Call) and child.name == 'RunPgenSusieBatch'
                            for child in walk(node)))

    def test_batch_scatter_rejects_unequal_planner_arrays(self):
        scatter = self.batch_scatter()
        values = {'Worklists': ['gs://bucket/one.tsv', 'gs://bucket/two.tsv'],
                  'PhenotypeBeds': ['gs://bucket/one.bed', 'gs://bucket/two.bed'],
                  'BatchPgenIndices': [1, 0], 'Chromosomes': ['1', '2']}
        for shortened in values:
            env = WDL.Env.Bindings()
            for name, value in values.items():
                item_type = WDL.Type.Int() if name == 'BatchPgenIndices' else (
                    WDL.Type.String() if name == 'Chromosomes' else WDL.Type.File())
                env = env.bind('PlanPgenSusie.' + name, WDL.Value.from_json(
                    WDL.Type.Array(item_type), value[:1] if name == shortened else value))
            with self.subTest(shortened=shortened), self.assertRaises(WDL.Error.EvalError):
                scatter.expr.eval(env, WDL.StdLib.Base('1.0'))

    def test_planner_files_survive_simulated_export_and_localization(self):
        # This checks our output expressions and File values, not Cromwell's
        # cloud uploader. The static glob guard checks the cloud-backend rule.
        with tempfile.TemporaryDirectory(prefix='susier-plan-transfer-') as temporary:
            root = pathlib.Path(temporary)
            producer = root / 'producer'
            batches = producer / 'plan/batches'
            batches.mkdir(parents=True)
            chromosomes = ['1'] * 8 + ['2'] * 4
            indices = [1] * 8 + [0] * 4  # Manifest order is chromosome 2, then 1.
            expected = {}
            for number in (12, 3, 11, 2, 10, 1, 9, 4, 8, 5, 7, 6):
                batch_id = 'batch_%06d' % number
                chromosome = chromosomes[number - 1]
                outcome_id = 'p%02d' % number
                expected[batch_id] = {
                    'tsv': 'batch_id\tchromosome\tphenotype_id\n%s\t%s\t%s\n' %
                           (batch_id, chromosome, outcome_id),
                    'bed': '#chr\tstart\tend\tID\tS1\nchr%s\t10\t11\t%s\t0.%02d\n' %
                           (chromosome, outcome_id, number)}
                for extension, content in expected[batch_id].items():
                    (batches / (batch_id + '.' + extension)).write_text(content)
            for name, extension in (('batch_worklists.txt', 'tsv'), ('batch_phenotype_beds.txt', 'bed')):
                (producer / 'plan' / name).write_text(''.join(
                    'plan/batches/batch_%06d.%s\n' % (number, extension)
                    for number in range(1, 13)))
            (producer / 'plan/batch_pgen_indices.json').write_text(json.dumps(indices))
            (producer / 'plan/batch_chromosomes.txt').write_text(''.join(chrom + '\n' for chrom in chromosomes))
            outputs = {decl.name: decl for decl in self.tasks['PlanPgenSusie'].outputs}
            stdlib = TaskOutputStdLib(producer)
            evaluated = {name: outputs[name].expr.eval(WDL.Env.Bindings(), stdlib).coerce(outputs[name].type)
                         for name in ('Worklists', 'PhenotypeBeds', 'BatchPgenIndices', 'Chromosomes')}
            for name, extension in (('Worklists', 'tsv'), ('PhenotypeBeds', 'bed')):
                self.assertEqual([pathlib.Path(value.value).name for value in evaluated[name].value],
                                 ['batch_%06d.%s' % (number, extension) for number in range(1, 13)])
                list_name = 'batch_worklists.txt' if name == 'Worklists' else 'batch_phenotype_beds.txt'
                self.assertEqual((producer / 'plan' / list_name).read_text().splitlines(),
                                 [str(pathlib.Path(value.value).relative_to(producer))
                                  for value in evaluated[name].value])
            self.assertEqual(evaluated['BatchPgenIndices'].json, indices)
            storage = root / 'storage'
            storage.mkdir()
            exported = {}
            cloud_values = {}
            for name in ('Worklists', 'PhenotypeBeds'):
                def export(file):
                    source = producer / file.value
                    self.assertIn(source, stdlib.collected_files,
                                  'The output is a listed path but was not selected for upload by glob')
                    uri = 'gs://test-bucket/call-PlanPgenSusie/glob-' + name + '/' + source.name
                    stored = storage / (name + '-' + source.name)
                    shutil.copyfile(source, stored)
                    exported[uri] = stored
                    return uri
                cloud_values[name] = WDL.Value.rewrite_paths(evaluated[name], export)
                self.assertIsInstance(cloud_values[name].type.item_type, WDL.Type.File)
                self.assertTrue(all('/plan/batches/' not in value.value for value in cloud_values[name].value))
            shutil.rmtree(producer)
            self.assertFalse(producer.exists())
            consumer = root / 'consumer'
            consumer.mkdir()
            files = self.manifest_files()
            env = WDL.Env.Bindings().bind('ChromosomePgens', files)
            for name, value in {**cloud_values, 'BatchPgenIndices': evaluated['BatchPgenIndices'],
                                'Chromosomes': evaluated['Chromosomes']}.items():
                env = env.bind('PlanPgenSusie.' + name, value)
            scatter = self.batch_scatter()
            planned = scatter.expr.eval(env, WDL.StdLib.Base('1.0'))
            call = next(node for node in walk(scatter)
                        if isinstance(node, WDL.Tree.Call) and node.name == 'RunPgenSusieBatch')
            index = next(node for node in walk(scatter)
                         if isinstance(node, WDL.Tree.Decl) and node.name == 'chromosome_index')
            self.assertEqual(len(planned.value), 12)
            for number, batch in enumerate(planned.value, 1):
                iteration = env.bind(scatter.variable, batch)
                chromosome_index = index.expr.eval(iteration, WDL.StdLib.Base('1.0'))
                self.assertEqual(chromosome_index.value, indices[number - 1])
                iteration = iteration.bind(index.name, chromosome_index)
                for name in ('pgen', 'pvar', 'psam'):
                    genotype = call.inputs[name].eval(iteration, WDL.StdLib.Base('1.0'))
                    self.assertIsInstance(genotype, WDL.Value.File)
                    self.assertEqual(genotype.value, files.value[indices[number - 1]].value[name].value)
                for name, extension in (('worklist', 'tsv'), ('phenotype_bed', 'bed')):
                    value = call.inputs[name].eval(iteration, WDL.StdLib.Base('1.0'))
                    self.assertIsInstance(value, WDL.Value.File)
                    def localize(file):
                        source = exported[file.value]
                        destination = consumer / pathlib.Path(file.value).name
                        shutil.copyfile(source, destination)
                        return str(destination)
                    localized = WDL.Value.rewrite_paths(value, localize)
                    self.assertEqual(pathlib.Path(localized.value).read_text(),
                                     expected['batch_%06d' % number][extension])

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
                        metadata = name in METADATA_INPUTS
                        expected_item_type = WDL.Type.String if metadata else WDL.Type.File
                        self.assertIsInstance(node.arguments[0].type.item_type, expected_item_type)
        self.assertCountEqual(writers, [
            ('AggregatePgenSusie', 'summary_files'),
            ('AggregatePgenSusie', 'credible_set_files'),
            ('AggregatePgenSusie', 'index_files'),
            ('AggregatePgenSusie', 'diagnostic_files'),
            *[('AggregatePgenSusie', name) for name in GENE_URI_INPUTS],
            *[('AggregatePgenSusie', name) for name in PGEN_METADATA_FLAGS.values()],
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
            if name in OPTIONAL_FILES and not optional:
                value = None
            elif isinstance(typ, WDL.Type.File):
                value = 'gs://bucket/' + name
            elif isinstance(typ, WDL.Type.Array) and isinstance(typ.item_type, WDL.Type.File):
                value = ['gs://bucket/' + name + '/one', 'gs://bucket/' + name + '/two']
            elif name == 'pgen_chromosomes':
                value = ['1', '2']
            elif name in METADATA_INPUTS:
                value = ['gs://bucket/outputs/' + name + '/one',
                         "gs://bucket/outputs/" + name + "/two 'quoted' $(touch PWNED) `touch PWNED`"]
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

    def execute(self, task, optional=True, broken=None, overrides=None, plan_failure=None):
        with tempfile.TemporaryDirectory(prefix='susier-wdl-') as temporary:
            cwd = pathlib.Path(temporary)
            raw = self.inputs(task, optional)
            raw.update(overrides or {})
            mapping = {}

            def localize(file):
                if file.value not in mapping:
                    dest = cwd / ("input 'one' 'two' $(touch PWNED) `touch PWNED` " + str(len(mapping)))
                    dest.write_text('fixture\n')
                    mapping[file.value] = str(dest)
                return mapping[file.value]

            localized = {name: WDL.Value.rewrite_paths(value, localize) for name, value in raw.items()}
            for name in METADATA_INPUTS & raw.keys():
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
            stub.write_text('#!' + sys.executable + '\n' + '''import json, os, pathlib, sys
args = sys.argv[1:]
lists = {}
for i, arg in enumerate(args):
    if arg in ''' + repr(LIST_FLAGS) + ''':
        lists[arg] = pathlib.Path(args[i+1]).read_text().splitlines()
pathlib.Path('capture.json').write_text(json.dumps({'args': args, 'lists': lists}))
if args[0].endswith('/filter_pgen_susie.R'):
    qc = pathlib.Path(args[args.index('--output-dir') + 1])
    qc.mkdir()
    for name in ('filtered.pgen', 'filtered.pvar', 'filtered.psam', 'variant_qc.tsv.gz', 'qc.json'):
        (qc / name).write_text('fixture\\n')
if args[0].endswith('/plan_pgen_susie.R'):
    plan = pathlib.Path(args[args.index('--output-dir') + 1])
    batches = plan / 'batches'
    batches.mkdir(parents=True)
    for number in (3, 1, 2):
        batch_id = 'batch_%06d' % number
        (batches / (batch_id + '.tsv')).write_text('batch_id\\tphenotype_id\\n' + batch_id + '\\tp' + str(number) + '\\n')
        (batches / (batch_id + '.bed')).write_text('#chr\\tstart\\tend\\tID\\tS1\\nchr1\\t10\\t11\\tp' + str(number) + '\\t1\\n')
    worklists = [str(batches / ('batch_%06d.tsv' % number)) for number in (1, 2, 3)]
    beds = [str(batches / ('batch_%06d.bed' % number)) for number in (1, 2, 3)]
    failure = os.environ.get('SUSIER_TEST_PLAN_FAILURE', '')
    if failure == 'missing_bed':
        pathlib.Path(beds[1]).unlink()
    elif failure == 'extra_bed':
        (batches / 'batch_999999.bed').write_text('extra\\n')
    elif failure == 'mismatched_stems':
        replacement = batches / 'batch_000099.bed'
        pathlib.Path(beds[1]).rename(replacement)
        beds = sorted([beds[0], beds[2], str(replacement)])
    elif failure == 'reversed_lists':
        worklists.reverse()
    elif failure == 'empty_lists':
        worklists, beds = [], []
    elif failure == 'empty_outputs':
        for path in batches.iterdir():
            path.unlink()
    elif failure:
        raise ValueError('Unknown planner fixture failure: ' + failure)
    (plan / 'batch_worklists.txt').write_text(''.join(path + '\\n' for path in worklists))
    (plan / 'batch_phenotype_beds.txt').write_text(''.join(path + '\\n' for path in beds))
    (plan / 'batch_chromosomes.txt').write_text('1\\n1\\n1\\n')
    (plan / 'batch_chromosome_indices.json').write_text('[0, 0, 0]')
    (plan / 'expected_outcomes.tsv').write_text('batch_id\\tphenotype_id\\n')
    (plan / 'batch_index.tsv').write_text('batch_id\\tchromosome\\tworklist\\tphenotype_bed\\n')
    (plan / 'chromosome_pgens.tsv').write_text('1\\tgs://bucket/a.pgen\\tgs://bucket/a.pvar\\tgs://bucket/a.psam\\n')
    (plan / 'unique_pgens.tsv').write_text('1\\tgs://bucket/a.pgen\\tgs://bucket/a.pvar\\tgs://bucket/a.psam\\n')
    (plan / 'batch_pgen_indices.json').write_text('[0, 0, 0]')
''')
            stub.chmod(0o755)
            result = subprocess.run(['bash', 'command.sh'], cwd=cwd,
                                    env={**os.environ, 'PATH': str(bindir) + os.pathsep + os.environ['PATH'],
                                         'SUSIER_TEST_PLAN_FAILURE': plan_failure or ''},
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
            if plan_failure:
                self.assertNotEqual(result.returncode, 0, 'Invalid planner outputs must stop the task')
                self.assertTrue(captured.exists(), 'The planner must run before output checks')
                self.assertIn('ERROR', result.stderr)
                self.assertIn('CHECK', result.stdout)
                self.assertNotIn('COMPLETE PGEN batch planning', result.stdout)
                return
            self.assertEqual(result.returncode, 0, result.stderr + '\n' + command)
            payload = json.loads(captured.read_text())
            args = payload['args']
            self.assertNotIn('gs://', ' '.join(args))
            if task.name == 'PlanPgenSusie':
                self.assertIn('CHECK', result.stdout)
                self.assertIn('COMPLETE PGEN batch planning', result.stdout)
                self.assertEqual(args[0], '/opt/r/scripts/plan_pgen_susie.R')
                self.assertNotIn('--chromosomes-file', args)
                self.assertEqual(args[args.index('--pgen-manifest') + 1], localized['pgen_manifest'].value)
                self.assertEqual(args[args.index('--batch-size') + 1], '100')
            elif task.name == 'RunPgenSusieBatch':
                self.assertEqual(args[0], '/opt/r/scripts/run_pgen_susie_batch.R')
                self.assertEqual('--memory-debug' in args, optional)
                self.assertIn('--missing-code=' + localized['missing_code'].value, args)
                self.assertNotIn('--missing-code', args)
            elif task.name == 'FilterPgen':
                self.assertEqual(args[0], '/opt/r/scripts/filter_pgen_susie.R')
                self.assertIn('--mac=' + str(localized['mac'].value), args)
                rate = next(arg.split('=', 1)[1] for arg in args if arg.startswith('--min-call-rate='))
                self.assertEqual(float(rate), localized['min_call_rate'].value)
            else:
                self.assertEqual(args[0], '/opt/r/scripts/aggregate_pgen_susie.R')
                for flag, paths in payload['lists'].items():
                    if flag in GENE_URI_FLAGS or flag in PGEN_METADATA_FLAGS:
                        self.assertEqual(paths, raw[{**GENE_URI_FLAGS, **PGEN_METADATA_FLAGS}[flag]].json)
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
            if task.name in ('PlanPgenSusie', 'RunPgenSusieBatch'):
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
            name = {'AggregatePgenSusie': 'expected_outcomes', 'FilterPgen': 'pgen'}.get(task.name, 'phenotype_bed')
            for failure in ('uri', 'missing'):
                with self.subTest(task=task.name, failure=failure):
                    self.execute(task, broken=(name, failure))

    def test_aggregate_file_lists_reject_unresolved_cloud_paths(self):
        for name in ('summary_files', 'credible_set_files', 'index_files', 'diagnostic_files'):
            with self.subTest(file_list=name):
                self.execute(self.tasks['AggregatePgenSusie'], broken=(name, 'uri'))

    def test_optional_sample_list_rejects_unresolved_or_missing_file(self):
        for failure in ('uri', 'missing'):
            with self.subTest(failure=failure):
                self.execute(self.tasks['RunPgenSusieBatch'], broken=('sample_list', failure))

    def test_genotype_files_reject_unresolved_or_missing_file(self):
        for name in ('pgen', 'pvar', 'psam'):
            for failure in ('uri', 'missing'):
                with self.subTest(file=name, failure=failure):
                    self.execute(self.tasks['RunPgenSusieBatch'], broken=(name, failure))

    def test_planner_rejects_unresolved_or_missing_manifest(self):
        for failure in ('uri', 'missing'):
            with self.subTest(failure=failure):
                self.execute(self.tasks['PlanPgenSusie'], broken=('pgen_manifest', failure))

    def test_planner_rejects_incomplete_or_misordered_output_files(self):
        for failure in ('missing_bed', 'extra_bed', 'mismatched_stems',
                        'reversed_lists', 'empty_lists', 'empty_outputs'):
            with self.subTest(failure=failure):
                self.execute(self.tasks['PlanPgenSusie'], plan_failure=failure)


if __name__ == '__main__':
    unittest.main()
