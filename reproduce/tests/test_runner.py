import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import sys
import tempfile
import unittest

RUNNER = Path(__file__).resolve().parents[1] / 'run_all.py'
spec = importlib.util.spec_from_file_location('repro_runner', RUNNER)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class RunnerTests(unittest.TestCase):
    def test_paths_cannot_escape(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                module.safe_relative(Path(tmp), '../outside.txt')

    def test_metric_extraction_and_comparison(self):
        with tempfile.TemporaryDirectory() as tmp:
            project = Path(tmp)
            p = project / 'results' / 'step36_final_validation' / 'final_validation.json'
            p.parent.mkdir(parents=True)
            p.write_text(json.dumps({'experiments': [
                {'experiment': '256_SNP', 'actual_auc': 0.6772, 'shadow_auc': 0.5036,
                 'auc_gap': 0.1736, 'calibrated_auc': 0.6885,
                 'bootstrap': {'gap_ci': [0.0976, 0.2509]}}], 'relatedness_sensitivity': []}))
            r = module.load_metrics(project)
            self.assertTrue(r['available'])
            self.assertEqual(len(r['experiments']), 1)
            c = module.compare(r, {'reference_results': {'256_SNP': {'actual_auc': 0.6772}}}, False)
            self.assertTrue(c[0]['within_tolerance'])
            self.assertFalse(c[0]['fresh_rerun'])

    def test_run_stage_output_refresh(self):
        with tempfile.TemporaryDirectory() as tmp:
            project = Path(tmp)
            (project / 'src').mkdir()
            script = project / 'src' / 'make.py'
            script.write_text('from pathlib import Path\np=Path("results/output.txt")\np.parent.mkdir(exist_ok=True)\np.write_text("new")\n')
            out = project / 'runs' / 'demo'
            out.mkdir(parents=True)
            stage = {'id': 'test', 'label': 'test', 'script': 'src/make.py',
                     'executor': 'windows_python', 'expected_outputs': ['results/output.txt'], 'timeout_seconds': 10}
            args = SimpleNamespace(python_win=Path(sys.executable), python_wsl='', distro='Ubuntu')
            r = module.run_stage(stage, project, out, args)
            self.assertEqual(r['status'], 'PASS', r)
            self.assertEqual((project / 'results/output.txt').read_text(), 'new')


if __name__ == '__main__':
    unittest.main()
