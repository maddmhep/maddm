"""Plot inputs: consistent fit blocks and residuals at actual AMS nodes."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('pbar_plot', ROOT/'plotting/plot_ams02_pbar_fluxes.py')
plot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plot)
FIXTURE = ROOT/'tests/fixtures/pbar/AMS02_pbar_fluxes.dat'


class PlotDataTests(unittest.TestCase):
    def test_all_fit_blocks_and_residual_denominator_at_data_points(self):
        raw = np.loadtxt(FIXTURE)
        for prediction, offset in [('model', 3), ('thermal', 6), ('best', 9), ('upper', 12)]:
            with self.subTest(prediction=prediction):
                flux = plot.read_fluxes(FIXTURE, prediction)
                self.assertEqual(len(flux['R']), 58)
                np.testing.assert_array_equal(flux['R'], raw[:, 0])
                np.testing.assert_array_equal(flux['secondary'], raw[:, offset])
                np.testing.assert_array_equal(flux['dm'], raw[:, offset+1])
                np.testing.assert_array_equal(flux['total'], raw[:, offset+2])
                np.testing.assert_allclose(flux['residual']*flux['total'], raw[:, offset+2]-raw[:, 1])
                np.testing.assert_allclose(flux['residual_error']*flux['total'], raw[:, 2])

    def test_invalid_components_cannot_make_a_plot(self):
        for failure in ['nan', 'mismatched_sum', 'missing_header']:
            with self.subTest(failure=failure), tempfile.TemporaryDirectory() as directory:
                data = np.loadtxt(FIXTURE)
                header = FIXTURE.read_text().splitlines()[0].lstrip('# ')
                if failure == 'nan': data[0, 4] = np.nan
                if failure == 'mismatched_sum': data[0, 5] *= 1.1
                path = Path(directory)/'fluxes.dat'
                np.savetxt(path, data, header='' if failure == 'missing_header' else header)
                with self.assertRaises(ValueError): plot.read_fluxes(path)


if __name__ == '__main__':
    unittest.main()
