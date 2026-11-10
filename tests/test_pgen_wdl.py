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
LIST_FLAGS = ('--worklists-file', '--summary-files', '--full-susie-files', '--credible-set-files', '--index-files', '--diagnostic-files') + tuple(GENE_URI_FLAGS) + tuple(PGEN_METADATA_FLAGS)
FILE_NAMES = {'pgen_manifest', 'phenotype_bed', 'qtl_results', 'gene_map', 'worklist', 'covariates', 'sample_list',
              'ancestry', 'variant_list', 'additional_genotypes', 'genotype_qc_report', 'prepared_batch_manifest', 'expected_outcomes', 'outcome_index', 'pgen', 'pvar', 'psam'}
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
             'for path in "${paths[@]}"; do '
             'if [[ -f "$path" ]]; then printf "%s\\0" "$path"; fi; done',
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

    def test_preparation_runs_once_per_unique_input_and_batches_receive_shards(self):
        self.assertIn('PreparePgenBatches', self.tasks,
                      'Extract batch PGEN files before the SuSiE batch tasks')
        self.assertNotIn('FilterPgen', self.tasks,
                         'The preparation task must perform shared QC and extraction together')
        call = next(node for node in walk(self.doc.workflow)
                    if isinstance(node, WDL.Tree.Call) and node.name == 'PreparePgenBatches')
        self.assertEqual(str(call.inputs['worklists'].type), 'Array[File]')
        batch = next(node for node in walk(self.doc.workflow)
                     if isinstance(node, WDL.Tree.Call) and node.name == 'RunPgenSusieBatch')
        for role, output in (('pgen', 'BatchPgens'), ('pvar', 'BatchPvars'),
                             ('psam', 'BatchPsams'), ('prepared_batch_manifest', 'Manifests')):
            self.assertIn('PreparePgenBatches.' + output, str(batch.inputs[role]))
            self.assertIsInstance(batch.inputs[role].type, WDL.Type.File)

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
        for name in ('summary_files', 'full_susie_files', 'credible_set_files', 'index_files', 'diagnostic_files'):
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

    def test_full_susie_files_remain_typed_from_batch_collection_to_aggregate(self):
        batch_outputs = {decl.name: decl for decl in self.tasks['RunPgenSusieBatch'].outputs}
        self.assertIn('FullSusieParquetFiles', batch_outputs,
                      'Collect full per-gene results as files for aggregation')
        declaration = batch_outputs['FullSusieParquetFiles']
        self.assertEqual(str(declaration.type), 'Array[File]')
        self.assertIsInstance(declaration.expr, WDL.Expr.Apply)
        self.assertEqual(declaration.expr.function_name, 'glob')
        with tempfile.TemporaryDirectory(prefix='susier-full-output-') as temporary:
            directory = pathlib.Path(temporary)
            selected = []
            for chromosome, gene in (('chr2', 'gene_b'), ('chr1', 'gene_a')):
                gene_dir = directory / 'batch/genes' / chromosome
                gene_dir.mkdir(parents=True)
                for suffix in ('.full_susie.parquet', '.susie.parquet', '.credible_sets.parquet'):
                    path = gene_dir / (gene + suffix)
                    path.write_text('fixture\n')
                    if suffix == '.full_susie.parquet':
                        selected.append(str(path))
            stdlib = TaskOutputStdLib(directory)
            collected = declaration.expr.eval(WDL.Env.Bindings(), stdlib)
            self.assertEqual(collected.json, sorted(selected))
            self.assertEqual({str(path) for path in stdlib.collected_files}, set(selected))

        task = self.tasks['AggregatePgenSusie']
        inputs = {decl.name: decl for decl in task.inputs}
        self.assertIn('full_susie_files', inputs)
        self.assertEqual(str(inputs['full_susie_files'].type), 'Array[File]')
        call = next(node for node in walk(self.doc.workflow)
                    if isinstance(node, WDL.Tree.Call) and node.name == 'AggregatePgenSusie')
        self.assertIn('full_susie_files', call.inputs)
        uris = ['gs://bucket/batch_1/gene_a.full_susie.parquet',
                'gs://bucket/batch_2/gene_b.full_susie.parquet']
        source = WDL.Value.from_json(WDL.Type.Array(WDL.Type.Array(WDL.Type.File())),
                                     [[uris[0]], [], [uris[1]]])
        env = WDL.Env.Bindings().bind('RunPgenSusieBatch.FullSusieParquetFiles', source)
        actual = call.inputs['full_susie_files'].eval(env, WDL.StdLib.Base('1.0'))
        self.assertIsInstance(actual.type, WDL.Type.Array)
        self.assertIsInstance(actual.type.item_type, WDL.Type.File)
        self.assertEqual(actual.json, uris)
        rewritten = WDL.Value.rewrite_paths(actual, lambda file: '/localized/' + pathlib.Path(file.value).name)
        self.assertEqual(rewritten.json, ['/localized/gene_a.full_susie.parquet',
                                          '/localized/gene_b.full_susie.parquet'])

        aggregate_outputs = {decl.name: decl for decl in task.outputs}
        self.assertIn('FullSuSiEParquet', aggregate_outputs)
        output = aggregate_outputs['FullSuSiEParquet']
        self.assertEqual(str(output.type), 'File')
        self.assertEqual(output.expr.eval(WDL.Env.Bindings(), WDL.StdLib.Base('1.0')).value,
                         'aggregate/full_susie.parquet')
        workflow_outputs = {decl.name: decl for decl in self.doc.workflow.outputs}
        self.assertIn('AggregateFullSuSiEParquet', workflow_outputs)
        output = workflow_outputs['AggregateFullSuSiEParquet']
        self.assertEqual(str(output.type), 'File')
        env = WDL.Env.Bindings().bind('AggregatePgenSusie.FullSuSiEParquet',
                                      WDL.Value.File('gs://bucket/aggregate/full_susie.parquet'))
        self.assertEqual(output.expr.eval(env, WDL.StdLib.Base('1.0')).value,
                         'gs://bucket/aggregate/full_susie.parquet')

    def saved_full_aggregation(self):
        path = ROOT / 'workflows/AggregateFullSusie.wdl'
        self.assertTrue(path.is_file(), 'Saved full outputs need an aggregation-only workflow')
        return WDL.load(str(path))

    def test_saved_full_aggregation_has_one_typed_aggregate_call(self):
        doc = self.saved_full_aggregation()
        self.assertEqual(doc.effective_wdl_version, '1.0')
        self.assertEqual(doc.workflow.name, 'AggregateFullSusie')
        self.assertEqual([task.name for task in doc.tasks], ['AggregateFullSusieFiles'])
        self.assertEqual(workflow_writers(doc.workflow), [])
        calls = [node for node in walk(doc.workflow) if isinstance(node, WDL.Tree.Call)]
        self.assertEqual([call.name for call in calls], ['AggregateFullSusieFiles'])
        workflow_inputs = {decl.name: decl for decl in doc.workflow.inputs}
        task = doc.tasks[0]
        task_inputs = {decl.name: decl for decl in task.inputs}
        self.assertFalse(task.postinputs, 'Keep File types until command rendering')
        writers = []
        command_text = ''.join(part for part in task.command.parts if isinstance(part, str))
        for part in task.command.parts:
            if isinstance(part, WDL.Expr.Placeholder):
                for node in walk(part.expr):
                    if isinstance(node, WDL.Expr.Apply) and node.function_name in WRITERS:
                        writers.append(node)
                        self.assertIs(part.expr, node)
                        self.assertIsInstance(node.arguments[0].type.item_type, WDL.Type.File)
        self.assertEqual([(node.function_name, str(node.arguments[0])) for node in writers],
                         [('write_lines', 'full_susie_files')])
        for forbidden in ('plink', 'run_pgen_susie_batch.R', 'susie_rss'):
            self.assertNotIn(forbidden, command_text)
        values = {
            'outcome_index': WDL.Value.File('gs://bucket/saved/outcome_index.tsv'),
            'full_susie_files': WDL.Value.from_json(WDL.Type.Array(WDL.Type.File()),
                                                   ['gs://bucket/saved/gene_a.full_susie.parquet'])}
        for name, value in values.items():
            with self.subTest(input=name):
                self.assertEqual(str(workflow_inputs[name].type), str(value.type).rstrip('+'))
                self.assertEqual(str(task_inputs[name].type), str(value.type).rstrip('+'))
                actual = calls[0].inputs[name].eval(WDL.Env.Bindings().bind(name, value),
                                                   WDL.StdLib.Base('1.0'))
                self.assertTrue(has_file_type(actual.type))
                self.assertEqual(actual.json, value.json)
        for name, expected in (('memory_gb', 8), ('disk_gb', 50),
                               ('docker_image', 'ghcr.io/aou-multiomics-analysis/susier/pgen:main')):
            self.assertEqual(workflow_inputs[name].expr.eval(WDL.Env.Bindings(),
                                                            WDL.StdLib.Base('1.0')).value, expected)
        task_outputs = {decl.name: decl for decl in task.outputs}
        workflow_outputs = {decl.name: decl for decl in doc.workflow.outputs}
        self.assertEqual(set(task_outputs), {'FullSuSiEParquet', 'TaskStdout', 'TaskStderr'})
        self.assertEqual(set(workflow_outputs), set(task_outputs))
        for name, output in workflow_outputs.items():
            self.assertIsInstance(output.type, WDL.Type.File)
            value = WDL.Value.File('gs://bucket/aggregate/' + name)
            actual = output.expr.eval(WDL.Env.Bindings().bind('AggregateFullSusieFiles.' + name, value),
                                      WDL.StdLib.Base('1.0'))
            self.assertEqual(actual.json, value.json)
        self.assertEqual(task_outputs['FullSuSiEParquet'].expr.eval(WDL.Env.Bindings(),
                                                                 WDL.StdLib.Base('1.0')).value,
                         'aggregate/full_susie.parquet')

    def test_saved_full_aggregation_localizes_files_and_rejects_cloud_uris(self):
        task = self.saved_full_aggregation().tasks[0]
        payload = self.execute(task)
        self.assertEqual([arg for arg in payload['args'][1:] if arg.startswith('--')],
                         ['--full-only', '--outcome-index', '--full-susie-files', '--output-dir'])
        self.assertEqual(set(payload['lists']), {'--full-susie-files'})
        for name in ('outcome_index', 'full_susie_files'):
            for failure in ('uri', 'missing'):
                with self.subTest(input=name, failure=failure):
                    self.execute(task, broken=(name, failure))

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

    def test_preparation_collects_each_shard_and_optional_qc_output(self):
        outputs = {decl.name: decl for decl in self.tasks['PreparePgenBatches'].outputs}
        with tempfile.TemporaryDirectory(prefix='susier-prepared-outputs-') as temporary:
            directory = pathlib.Path(temporary)
            batches = directory / 'prepared/batches'
            batches.mkdir(parents=True)
            for number in (12, 2, 1):
                for extension in ('pgen', 'pvar', 'psam', 'preparation.json'):
                    (batches / ('batch_%06d.' % number + extension)).write_text('%d:%s\n' % (number, extension))
            stdlib = TaskOutputStdLib(directory)
            exported = {}
            cloud_outputs = {}
            storage = directory / 'storage'
            storage.mkdir()
            for output, extension in (('BatchPgens', 'pgen'), ('BatchPvars', 'pvar'),
                                      ('BatchPsams', 'psam'), ('Manifests', 'preparation.json')):
                declaration = outputs[output]
                self.assertEqual(str(declaration.type), 'Array[File]')
                self.assertEqual(declaration.expr.function_name, 'glob')
                values = declaration.expr.eval(WDL.Env.Bindings(), stdlib)
                self.assertEqual([pathlib.Path(value.value).name for value in values.value],
                                 ['batch_%06d.%s' % (number, extension) for number in (1, 2, 12)])
                self.assertTrue(all(pathlib.Path(value.value) in stdlib.collected_files
                                    for value in values.value))
                def export(file):
                    source = pathlib.Path(file.value)
                    stored = storage / source.name
                    shutil.copyfile(source, stored)
                    uri = 'gs://test-bucket/call-PreparePgenBatches/glob/' + source.name
                    exported[uri] = stored
                    return uri
                cloud_outputs[extension] = WDL.Value.rewrite_paths(values, export)
            shutil.rmtree(batches)
            consumer = directory / 'consumer'
            consumer.mkdir()
            for extension, values in cloud_outputs.items():
                self.assertIsInstance(values.type.item_type, WDL.Type.File)
                def localize(file):
                    destination = consumer / pathlib.Path(file.value).name
                    shutil.copyfile(exported[file.value], destination)
                    return str(destination)
                localized = WDL.Value.rewrite_paths(values, localize)
                self.assertEqual([pathlib.Path(file.value).read_text() for file in localized.value],
                                 ['%d:%s\n' % (number, extension) for number in (1, 2, 12)])
            for output in ('Counts', 'Reports'):
                self.assertEqual(outputs[output].expr.eval(WDL.Env.Bindings(), stdlib).json, [])
            qc = directory / 'prepared/qc'
            qc.mkdir()
            (qc / 'variant_qc.tsv.gz').write_text('counts\n')
            (qc / 'qc.json').write_text('{}\n')
            for output in ('Counts', 'Reports'):
                self.assertEqual(len(outputs[output].expr.eval(WDL.Env.Bindings(), stdlib).value), 1)
            for mac, rate, enabled in ((0, 0.0, False), (0, .95, True), (5, 0.0, True)):
                env = WDL.Env.Bindings().bind('mac', WDL.Value.Int(mac))
                env = env.bind('min_call_rate', WDL.Value.Float(rate))
                for output, role in (('Pgen', 'pgen'), ('Pvar', 'pvar'), ('Psam', 'psam')):
                    source = WDL.Value.File('gs://source/original.' + role)
                    actual = outputs[output].expr.eval(env.bind(role, source), stdlib).coerce(outputs[output].type)
                    self.assertIsInstance(actual, WDL.Value.File)
                    self.assertEqual(actual.value, 'prepared/qc/filtered.' + role if enabled else source.value)

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
        self.assertIsInstance(inputs['MAC'].type, WDL.Type.Int)
        self.assertEqual(inputs['MAC'].expr.eval(WDL.Env.Bindings(), WDL.StdLib.Base('1.0')).value, 0)
        self.assertEqual(inputs['MinCallRate'].expr.eval(WDL.Env.Bindings(), WDL.StdLib.Base('1.0')).value, 0.95)
        preparation = next(node for node in walk(self.doc.workflow)
                           if isinstance(node, WDL.Tree.Call) and node.name == 'PreparePgenBatches')
        batch = next(node for node in walk(self.doc.workflow)
                     if isinstance(node, WDL.Tree.Call) and node.name == 'RunPgenSusieBatch')
        for threshold, rate in ((0, 0.0), (0, .95), (5, 0.0)):
            with self.subTest(mac=threshold, min_call_rate=rate):
                env = WDL.Env.Bindings().bind('MAC', WDL.Value.Int(threshold))
                env = env.bind('MinCallRate', WDL.Value.Float(rate))
                mac = preparation.inputs['mac'].eval(env, WDL.StdLib.Base('1.0'))
                call_rate = preparation.inputs['min_call_rate'].eval(env, WDL.StdLib.Base('1.0'))
                self.assertEqual(mac.value, threshold)
                self.assertEqual(call_rate.value, rate)
                self.assertEqual(batch.inputs['mac'].eval(env, WDL.StdLib.Base('1.0')).value, 0,
                                 'SuSiE tasks must not repeat genotype QC')
                payload = self.execute(self.tasks['PreparePgenBatches'],
                                       overrides={'mac': mac, 'min_call_rate': call_rate})
                self.assertIn('--mac=' + str(threshold), payload['args'])

    def test_genotype_cache_defaults_and_batch_forwarding(self):
        inputs = {decl.name: decl for decl in self.doc.workflow.inputs}
        self.assertIsInstance(inputs['GenotypeCacheMode'].type, WDL.Type.String)
        self.assertEqual(inputs['GenotypeCacheMode'].expr.eval(
            WDL.Env.Bindings(), WDL.StdLib.Base('1.0')).value, 'auto')
        self.assertIsInstance(inputs['GenotypeCacheMemoryGB'].type, WDL.Type.Float)
        self.assertEqual(inputs['GenotypeCacheMemoryGB'].expr.eval(
            WDL.Env.Bindings(), WDL.StdLib.Base('1.0')).value, 0.0)
        batch = next(node for node in walk(self.doc.workflow)
                     if isinstance(node, WDL.Tree.Call) and node.name == 'RunPgenSusieBatch')
        for mode in ('auto', 'whole', 'rolling', 'none'):
            with self.subTest(mode=mode):
                env = WDL.Env.Bindings().bind('GenotypeCacheMode', WDL.Value.String(mode))
                env = env.bind('GenotypeCacheMemoryGB', WDL.Value.Float(3.25))
                self.assertEqual(batch.inputs['genotype_cache_mode'].eval(
                    env, WDL.StdLib.Base('1.0')).value, mode)
                self.assertEqual(batch.inputs['cache_memory_gb'].eval(
                    env, WDL.StdLib.Base('1.0')).value, 3.25)

    def test_genotype_cache_cli_receives_task_memory_and_quoted_mode(self):
        task = self.tasks['RunPgenSusieBatch']
        for mode in ('auto', 'whole', 'rolling', 'none',
                     "auto 'quoted' $(touch PWNED) `touch PWNED`"):
            with self.subTest(mode=mode):
                payload = self.execute(task, overrides={
                    'genotype_cache_mode': WDL.Value.String(mode),
                    'cache_memory_gb': WDL.Value.Float(3.25),
                    'memory_gb': WDL.Value.Int(24),
                })
                args = payload['args']
                self.assertEqual(args[args.index('--genotype-cache-mode') + 1], mode)
                self.assertEqual(float(args[args.index('--genotype-cache-memory-gb') + 1]), 3.25)
                self.assertEqual(float(args[args.index('--job-memory-gb') + 1]), 24)

    def test_genotype_cache_reports_are_collected_as_task_files(self):
        task_outputs = {decl.name: decl for decl in self.tasks['RunPgenSusieBatch'].outputs}
        report = task_outputs['GenotypeCacheReport']
        self.assertIsInstance(report.type, WDL.Type.File)
        self.assertEqual(report.expr.eval(WDL.Env.Bindings(), WDL.StdLib.Base('1.0')).value,
                         'batch/batch.genotype_cache.json')
        workflow_outputs = {decl.name: decl for decl in self.doc.workflow.outputs}
        reports = workflow_outputs['GenotypeCacheReports']
        self.assertEqual(str(reports.type), 'Array[File]')
        uris = ['gs://outputs/batch_1/cache.json', 'gs://outputs/batch_2/cache.json']
        values = WDL.Value.from_json(WDL.Type.Array(WDL.Type.File()), uris)
        env = WDL.Env.Bindings().bind('RunPgenSusieBatch.GenotypeCacheReport', values)
        self.assertEqual(reports.expr.eval(env, WDL.StdLib.Base('1.0')).json, uris)

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

    def test_manifest_paths_become_files_before_preparation_localization(self):
        files = self.manifest_files()
        inputs = {decl.name: decl for decl in self.tasks['PreparePgenBatches'].inputs}
        for name in ('pgen', 'pvar', 'psam'):
            self.assertIsInstance(inputs[name].type, WDL.Type.File)
        for file_bundle in files.value:
            self.execute(self.tasks['PreparePgenBatches'], overrides={
                name: file_bundle.value[name] for name in ('pgen', 'pvar', 'psam')})

    def test_preparation_groups_worklists_without_localizing_other_inputs(self):
        sources = self.manifest_files()
        scatter = next(node for node in walk(self.doc.workflow)
                       if isinstance(node, WDL.Tree.Scatter) and node.variable == 'pgen_index')
        calls = [node for node in walk(scatter) if isinstance(node, WDL.Tree.Call)]
        self.assertEqual([call.name for call in calls], ['PreparePgenBatches'])
        selection = next(node for node in walk(scatter)
                         if isinstance(node, WDL.Tree.Scatter) and node.variable == 'preparation_batch_index')
        self.assertFalse(any(isinstance(node, WDL.Tree.Call) for node in walk(selection)),
                         'Worklist grouping must not create additional compute calls')
        selected = next(node for node in walk(selection)
                        if isinstance(node, WDL.Tree.Decl) and node.name == 'preparation_worklist')
        self.assertIsInstance(selected.type, WDL.Type.File)
        rows = [[value.value['chromosome'].value,
                 *[value.value[role].value for role in ('pgen', 'pvar', 'psam')]]
                for value in sources.value]
        env = WDL.Env.Bindings().bind('pgen_rows', WDL.Value.from_json(
            WDL.Type.Array(WDL.Type.Array(WDL.Type.String())), rows))
        env = env.bind('PlanPgenSusie.PgenBatchIndices', WDL.Value.from_json(
            WDL.Type.Array(WDL.Type.Array(WDL.Type.Int())), [[2], [0, 1]]))
        worklists = ['gs://plan/batch_%06d.tsv' % i for i in (1, 2, 3)]
        env = env.bind('PlanPgenSusie.Worklists', WDL.Value.from_json(
            WDL.Type.Array(WDL.Type.File()), worklists))
        self.assertEqual(scatter.expr.eval(env, WDL.StdLib.Base('1.0')).json, [0, 1])
        for source_index, expected_indices in enumerate(([2], [0, 1])):
            iteration = env.bind(scatter.variable, WDL.Value.Int(source_index))
            for declaration in scatter.body:
                if isinstance(declaration, WDL.Tree.Decl) and declaration.name in ('pgen_row', 'chromosome_pgen'):
                    iteration = iteration.bind(declaration.name,
                        declaration.expr.eval(iteration, WDL.StdLib.Base('1.0')).coerce(declaration.type))
            selected_files = []
            for index in selection.expr.eval(iteration, WDL.StdLib.Base('1.0')).value:
                value = selected.expr.eval(iteration.bind(selection.variable, index), WDL.StdLib.Base('1.0'))
                self.assertIsInstance(value, WDL.Value.File)
                selected_files.append(value)
            iteration = iteration.bind('preparation_worklist', WDL.Value.Array(WDL.Type.File(), selected_files))
            self.assertEqual(calls[0].inputs['worklists'].eval(iteration, WDL.StdLib.Base('1.0')).json,
                             [worklists[index] for index in expected_indices])
            for role in ('pgen', 'pvar', 'psam'):
                actual = calls[0].inputs[role].eval(iteration, WDL.StdLib.Base('1.0'))
                self.assertIsInstance(actual, WDL.Value.File)
                self.assertEqual(actual.value, sources.value[source_index].value[role].value)

    def batch_environment(self, chromosomes=('1', '1', '2'), indices=(1, 1, 0),
                          within=(0, 1, 0), mac=0, rate=.95):
        env = WDL.Env.Bindings().bind('ChromosomePgens', self.manifest_files())
        env = env.bind('MAC', WDL.Value.Int(mac)).bind('MinCallRate', WDL.Value.Float(rate))
        for name, values, typ in (
            ('BatchPgenIndices', indices, WDL.Type.Int()),
            ('BatchWithinPgenIndices', within, WDL.Type.Int()),
            ('Chromosomes', chromosomes, WDL.Type.String()),
            ('Worklists', ['gs://plan/batch_%06d.tsv' % i for i in range(1, len(indices) + 1)], WDL.Type.File()),
            ('PhenotypeBeds', ['gs://plan/batch_%06d.bed' % i for i in range(1, len(indices) + 1)], WDL.Type.File()),
        ):
            env = env.bind('PlanPgenSusie.' + name, WDL.Value.from_json(WDL.Type.Array(typ), list(values)))
        for output, extension in (('BatchPgens', 'pgen'), ('BatchPvars', 'pvar'),
                                  ('BatchPsams', 'psam'), ('Manifests', 'preparation.json')):
            groups = [[] for _ in range(max(indices) + 1)]
            for i, source_index in enumerate(indices, 1):
                groups[source_index].append('gs://prepared/source_%d/batch_%06d.%s' %
                                             (source_index, i, extension))
            env = env.bind('PreparePgenBatches.' + output, WDL.Value.from_json(
                WDL.Type.Array(WDL.Type.Array(WDL.Type.File())), groups))
        reports = [['gs://qc/source_%d/qc.json' % index] if mac > 0 or rate > 0 else []
                   for index in range(max(indices) + 1)]
        env = env.bind('PreparePgenBatches.Reports', WDL.Value.from_json(
            WDL.Type.Array(WDL.Type.Array(WDL.Type.File())), reports))
        return env

    def batch_iteration(self, env, batch):
        iteration = env.bind(self.batch_scatter().variable, batch)
        for declaration in self.batch_scatter().body:
            if isinstance(declaration, WDL.Tree.Decl):
                iteration = iteration.bind(declaration.name,
                    declaration.expr.eval(iteration, WDL.StdLib.Base('1.0')).coerce(declaration.type))
            elif isinstance(declaration, WDL.Tree.Conditional):
                enabled = declaration.expr.eval(iteration, WDL.StdLib.Base('1.0')).value
                for conditional in declaration.body:
                    if isinstance(conditional, WDL.Tree.Decl):
                        value = conditional.expr.eval(iteration, WDL.StdLib.Base('1.0')) if enabled else WDL.Value.Null()
                        iteration = iteration.bind(conditional.name, value)
        return iteration

    def test_shared_qc_reports_and_matching_shards_reach_each_batch(self):
        call = next(node for node in walk(self.doc.workflow)
                    if isinstance(node, WDL.Tree.Call) and node.name == 'RunPgenSusieBatch')
        for mac, rate in ((0, 0.0), (0, .95), (5, 0.0)):
            env = self.batch_environment(mac=mac, rate=rate)
            for number, batch in enumerate(self.batch_scatter().expr.eval(env, WDL.StdLib.Base('1.0')).value, 1):
                iteration = self.batch_iteration(env, batch)
                source_index = (1, 1, 0)[number - 1]
                self.assertEqual(iteration['chromosome_index'].value, source_index)
                self.assertEqual(iteration['within_pgen_index'].value, (0, 1, 0)[number - 1])
                for role, extension in (('pgen', 'pgen'), ('pvar', 'pvar'), ('psam', 'psam'),
                                        ('prepared_batch_manifest', 'preparation.json')):
                    actual = call.inputs[role].eval(iteration, WDL.StdLib.Base('1.0'))
                    self.assertIsInstance(actual, WDL.Value.File)
                    self.assertEqual(actual.value, 'gs://prepared/source_%d/batch_%06d.%s' %
                                     (source_index, number, extension))
                report = call.inputs['genotype_qc_report'].eval(iteration, WDL.StdLib.Base('1.0'))
                self.assertEqual(report.json, 'gs://qc/source_%d/qc.json' % source_index
                                 if mac > 0 or rate > 0 else None)
                self.assertEqual(call.inputs['mac'].eval(iteration, WDL.StdLib.Base('1.0')).value, 0)
                for name, extension in (('worklist', 'tsv'), ('phenotype_bed', 'bed')):
                    actual = call.inputs[name].eval(iteration, WDL.StdLib.Base('1.0'))
                    self.assertIsInstance(actual, WDL.Value.File)
                    self.assertEqual(actual.value, 'gs://plan/batch_%06d.%s' % (number, extension))

    def batch_scatter(self):
        return next(node for node in walk(self.doc.workflow)
                    if isinstance(node, WDL.Tree.Scatter)
                    and any(isinstance(child, WDL.Tree.Call) and child.name == 'RunPgenSusieBatch'
                            for child in walk(node)))

    def test_batch_scatter_rejects_unequal_planner_arrays(self):
        scatter = self.batch_scatter()
        values = {'Worklists': ['gs://bucket/one.tsv', 'gs://bucket/two.tsv'],
                  'PhenotypeBeds': ['gs://bucket/one.bed', 'gs://bucket/two.bed'],
                  'BatchPgenIndices': [1, 0], 'BatchWithinPgenIndices': [0, 0], 'Chromosomes': ['1', '2']}
        for shortened in values:
            env = WDL.Env.Bindings()
            for name, value in values.items():
                item_type = WDL.Type.Int() if name in ('BatchPgenIndices', 'BatchWithinPgenIndices') else (
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
            env = self.batch_environment(chromosomes, indices, list(range(8)) + list(range(4)))
            for name, value in {**cloud_values, 'BatchPgenIndices': evaluated['BatchPgenIndices'],
                                'Chromosomes': evaluated['Chromosomes']}.items():
                env = env.bind('PlanPgenSusie.' + name, value)
            scatter = self.batch_scatter()
            planned = scatter.expr.eval(env, WDL.StdLib.Base('1.0'))
            call = next(node for node in walk(scatter)
                        if isinstance(node, WDL.Tree.Call) and node.name == 'RunPgenSusieBatch')
            self.assertEqual(len(planned.value), 12)
            for number, batch in enumerate(planned.value, 1):
                iteration = self.batch_iteration(env, batch)
                chromosome_index = iteration['chromosome_index']
                self.assertEqual(chromosome_index.value, indices[number - 1])
                for name in ('pgen', 'pvar', 'psam'):
                    genotype = call.inputs[name].eval(iteration, WDL.StdLib.Base('1.0'))
                    self.assertIsInstance(genotype, WDL.Value.File)
                    self.assertEqual(genotype.value, 'gs://prepared/source_%d/batch_%06d.%s' %
                                     (indices[number - 1], number, name))
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
            ('PreparePgenBatches', 'worklists'),
            ('AggregatePgenSusie', 'summary_files'),
            ('AggregatePgenSusie', 'full_susie_files'),
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

    def execute(self, task, optional=True, broken=None, overrides=None, plan_failure=None,
                preparation_failure=None):
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
                elif broken[1] == 'empty':
                    localized[name] = WDL.Value.Array(WDL.Type.File(), [])
                else:
                    value = localized[name]
                    path = value.value[0].value if isinstance(value, WDL.Value.Array) else value.value
                    pathlib.Path(path).unlink()
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
if args[0].endswith('/prepare_pgen_susie.R'):
    prepared = pathlib.Path(args[args.index('--output-dir') + 1])
    batches = prepared / 'batches'
    batches.mkdir(parents=True)
    worklists = pathlib.Path(args[args.index('--worklists-file') + 1]).read_text().splitlines()
    for number in range(1, len(worklists) + 1):
        for extension in ('pgen', 'pvar', 'psam', 'preparation.json'):
            (batches / ('batch_%06d.' % number + extension)).write_text('fixture\\n')
    failure = os.environ.get('SUSIER_TEST_PREPARATION_FAILURE', '')
    if failure == 'missing_pvar':
        (batches / 'batch_000001.pvar').unlink()
    elif failure == 'extra_pgen':
        (batches / 'batch_999999.pgen').write_text('extra\\n')
    elif failure == 'mismatched_stems':
        (batches / 'batch_000001.preparation.json').rename(batches / 'batch_999999.preparation.json')
    elif failure == 'empty_outputs':
        for path in batches.iterdir():
            path.unlink()
    elif failure:
        raise ValueError('Unknown preparation fixture failure: ' + failure)
    (prepared / 'batch_worklists.txt').write_text(''.join(path + '\\n' for path in worklists))
    qc = prepared / 'qc'
    qc.mkdir()
    for extension in ('pgen', 'pvar', 'psam'):
        (qc / ('filtered.' + extension)).write_text('fixture\\n')
    (qc / 'variant_qc.tsv.gz').write_text('fixture\\n')
    (qc / 'qc.json').write_text('{}\\n')
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
    (plan / 'pgen_batch_indices.json').write_text('[[0, 1, 2]]')
    (plan / 'batch_within_pgen_indices.json').write_text('[0, 1, 2]')
''')
            stub.chmod(0o755)
            result = subprocess.run(['bash', 'command.sh'], cwd=cwd,
                                    env={**os.environ, 'PATH': str(bindir) + os.pathsep + os.environ['PATH'],
                                         'SUSIER_TEST_PLAN_FAILURE': plan_failure or '',
                                         'SUSIER_TEST_PREPARATION_FAILURE': preparation_failure or ''},
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
            if preparation_failure:
                self.assertNotEqual(result.returncode, 0, 'Incomplete shard outputs must stop the task')
                self.assertTrue(captured.exists(), 'Preparation must run before output checks')
                self.assertIn('ERROR', result.stderr)
                self.assertIn('CHECK', result.stdout)
                self.assertNotIn('COMPLETE shared PGEN QC and batch preparation', result.stdout)
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
            elif task.name == 'PreparePgenBatches':
                self.assertEqual(args[0], '/opt/r/scripts/prepare_pgen_susie.R')
                self.assertEqual(payload['lists']['--worklists-file'], localized['worklists'].json)
                self.assertIn('--mac=' + str(localized['mac'].value), args)
                rate = next(arg.split('=', 1)[1] for arg in args if arg.startswith('--min-call-rate='))
                self.assertEqual(float(rate), localized['min_call_rate'].value)
            else:
                self.assertEqual(args[0], '/opt/r/scripts/aggregate_pgen_susie.R')
                self.assertIn('--full-susie-files', args)
                for flag, paths in payload['lists'].items():
                    if flag in GENE_URI_FLAGS or flag in PGEN_METADATA_FLAGS:
                        self.assertEqual(paths, raw[{**GENE_URI_FLAGS, **PGEN_METADATA_FLAGS}[flag]].json)
                    else:
                        self.assertEqual(paths, localized[flag[2:].replace('-', '_')].json)
                        for path in paths:
                            self.assertTrue(pathlib.Path(path).is_file())
                for flag, name in GENE_URI_FLAGS.items():
                    if name in raw:
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
            name = {'AggregatePgenSusie': 'expected_outcomes', 'PreparePgenBatches': 'pgen'}.get(task.name, 'phenotype_bed')
            for failure in ('uri', 'missing'):
                with self.subTest(task=task.name, failure=failure):
                    self.execute(task, broken=(name, failure))

    def test_aggregate_file_lists_reject_unresolved_cloud_paths(self):
        for name in ('summary_files', 'full_susie_files', 'credible_set_files', 'index_files', 'diagnostic_files'):
            for failure in ('uri', 'missing'):
                with self.subTest(file_list=name, failure=failure):
                    self.execute(self.tasks['AggregatePgenSusie'], broken=(name, failure))

    def test_preparation_worklists_reject_unresolved_or_missing_files(self):
        for failure in ('uri', 'missing', 'empty'):
            with self.subTest(failure=failure):
                self.execute(self.tasks['PreparePgenBatches'], broken=('worklists', failure))

    def test_preparation_rejects_incomplete_or_mismatched_shard_files(self):
        for failure in ('missing_pvar', 'extra_pgen', 'mismatched_stems', 'empty_outputs'):
            with self.subTest(failure=failure):
                self.execute(self.tasks['PreparePgenBatches'], preparation_failure=failure)

    def test_optional_sample_list_rejects_unresolved_or_missing_file(self):
        for failure in ('uri', 'missing'):
            with self.subTest(failure=failure):
                self.execute(self.tasks['RunPgenSusieBatch'], broken=('sample_list', failure))

    def test_genotype_files_reject_unresolved_or_missing_file(self):
        for name in ('pgen', 'pvar', 'psam', 'prepared_batch_manifest'):
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
