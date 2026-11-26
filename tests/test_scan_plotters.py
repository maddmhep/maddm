"""Regression coverage for optional analyses and scans over several inputs."""
import contextlib
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT/'plotting'/(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SUMMARY = load('plot_sigmav_limits')
COUPLING = load('plot_coupling_mass')


def scan_file(path, columns):
    with path.open('w') as stream:
        for i, name in enumerate(columns, 1):
            stream.write('# [%02d] : %s\n' % (i, name))
        for row in zip(*columns.values()):
            stream.write(' '.join(map(str, row))+'\n')


class ScanPlotsTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        output = contextlib.redirect_stdout(io.StringIO())
        output.__enter__()
        self.addCleanup(output.__exit__, None, None, None)

    def base(self):
        return dict(run=np.array(['run_01','run_02','run_03']),
                    **{'mass#32':np.array([.45,4.5,45.]),
                       'mass#9000005':np.array([1.,10.,100.])})

    def test_two_mass_direct_only_includes_electrons_and_all_nuclear_limits(self):
        cols = self.base()
        cols.update(sigmaN_SD_p=np.array([1e-40,2e-40,3e-40]),
                    lim_sigmaN_SD_p=np.array([2e-41,3e-41,4e-41]),
                    lim_sigmaN_SD_p_LZ=np.array([-1.,1e-41,2e-41]),
                    sigma_e=np.array([1e-38,2e-38,3e-38]),
                    lim_sigma_e=np.array([-1.,3e-39,5e-39]),
                    lim_sigma_e_rescaled=np.array([-1.,2e-39,4e-39]),
                    lim_sigma_e_Xenon10=np.array([1e-38,2e-38,-1.]))
        manifest = SUMMARY.generate(cols, self.root, 'dd')
        self.assertEqual(manifest['x_column'], 'mass#9000005')
        self.assertIn('electron_recoil', manifest['plots'])
        for result in manifest['plots'].values():
            colors=[style['color'] for style in result['styles']]
            self.assertEqual(len(colors),len(set(colors)))
        self.assertIn('lim_sigmaN_SD_p_LZ', manifest['plots']['dd_limits']['columns'])
        self.assertIn('relic_density', manifest['skipped'])
        path = self.root/'direct.txt'
        scan_file(path, cols)
        curves, meta = COUPLING.generate(path, coupling=.2)
        expected = .2*(cols['lim_sigmaN_SD_p']/cols['sigmaN_SD_p'])**.25
        np.testing.assert_allclose(curves['SD_p'], expected)
        self.assertIn('electron_Xenon10', meta['envelope_inputs'])
        self.assertNotIn('electron_rescaled', meta['envelope_inputs'])
        self.assertEqual(meta['other_mass_columns'], ['mass#32'])

    def test_damicm_published_and_rescaled_are_alternatives(self):
        """only one convention of each published/rescaled pair (XENONnT, DAMIC-M) enters the envelope, following
        --electron-limit-column"""
        cols = self.base()
        cols.update(sigma_e=np.array([1e-38,2e-38,3e-38]),
                    lim_sigma_e=np.array([-1.,3e-39,5e-39]),
                    lim_sigma_e_rescaled=np.array([-1.,2e-39,4e-39]),
                    lim_sigma_e_DAMICM=np.array([1e-38,2e-38,3e-38]),
                    lim_sigma_e_DAMICM_rescaled=np.array([1.3e-38,2.6e-38,3.9e-38]),
                    lim_sigma_e_DS50=np.array([-1.,5e-38,6e-38]))
        path = self.root/'damicm.txt'
        scan_file(path, cols)
        _, meta = COUPLING.generate(path, coupling=.2)
        self.assertIn('electron_DAMICM', meta['envelope_inputs'])
        self.assertNotIn('electron_DAMICM_rescaled', meta['envelope_inputs'])
        self.assertIn('electron_DS50', meta['envelope_inputs'])
        _, meta = COUPLING.generate(path, coupling=.2, electron_limit_column='lim_sigma_e_rescaled')
        self.assertIn('electron_DAMICM_rescaled', meta['envelope_inputs'])
        self.assertNotIn('electron_DAMICM', meta['envelope_inputs'])
        self.assertIn('electron_rescaled', meta['envelope_inputs'])
        self.assertNotIn('electron', meta['envelope_inputs'])
        self.assertEqual(COUPLING.electron_label('lim_sigma_e_DAMICM_rescaled'), 'DAMIC-M; RAPIDD (QCDark2)')
        self.assertEqual(SUMMARY.electron_label('lim_sigma_e_DAMICM_rescaled'), 'DAMIC-M; RAPIDD (QCDark2)')

    def test_direct_nucleon_uses_equivalent_cross_sections(self):
        """LZ / XLZD limits and the neutrino fog against sigmaN_*_eq, PICO-60 against sigmaN_SD_p; the sigma_p,
        sigma_n limits of the same scan are not plotted"""
        cols = self.base()
        n = lambda *v: np.array(v, float)
        cols.update(sigmaN_SI_p=n(1e-46, 2e-46, 3e-46), sigmaN_SI_n=n(1e-46, 2e-46, 3e-46),
                    sigmaN_SD_p=n(1e-40, 2e-40, 3e-40), sigmaN_SD_n=n(1e-41, 2e-41, 3e-41),
                    lim_sigmaN_SI_p=n(1e-47, 1e-47, 1e-47), lim_sigmaN_SI_n=n(1e-47, 1e-47, 1e-47),
                    lim_sigmaN_SD_p=n(2e-41, 3e-41, 4e-41), lim_sigmaN_SD_p_LZ=n(1e-41, 1e-41, 2e-41),
                    lim_sigmaN_SD_n=n(1e-42, 1e-42, 1e-42),
                    sigmaN_SI_eq=n(5e-47, 6e-47, 7e-47), lim_sigmaN_SI_eq=n(1e-47, 1e-47, 1e-47),
                    sigmaN_SD_n_eq=n(4e-41, 5e-41, 6e-41), lim_sigmaN_SD_n_eq=n(1e-42, 1e-42, 1e-42),
                    sigmaN_SD_p_eq=n(9e-40, 9e-40, 9e-40), lim_sigmaN_SD_p_eq=n(1e-41, 1e-41, 2e-41),
                    lim_sigmaN_SI_eq_XLZD=n(1e-49, 1e-49, 1e-49), neutrino_fog_SI=n(1e-49, 2e-49, 3e-49))
        manifest = SUMMARY.generate(cols, self.root, 'eq')
        plotted = manifest['plots']['dd_comparison']['columns']
        self.assertEqual(sorted(manifest['plots']['dd_limits']['columns']),
                         sorted(['lim_sigmaN_SI_eq', 'lim_sigmaN_SI_eq_XLZD', 'lim_sigmaN_SD_n_eq',
                                 'lim_sigmaN_SD_p_eq', 'lim_sigmaN_SD_p']))
        for key in ('sigmaN_SI_eq', 'sigmaN_SD_n_eq', 'sigmaN_SD_p_eq', 'sigmaN_SD_p'):
            self.assertIn(key, plotted)
        for key in ('sigmaN_SI_p', 'sigmaN_SD_n', 'lim_sigmaN_SI_p', 'lim_sigmaN_SD_p_LZ'):
            self.assertNotIn(key, plotted)
        path = self.root/'direct_eq.txt'
        scan_file(path, cols)
        curves, meta = COUPLING.generate(path, coupling=.2)
        np.testing.assert_allclose(curves['SI_eq'], .2*(cols['lim_sigmaN_SI_eq']/cols['sigmaN_SI_eq'])**.25)
        np.testing.assert_allclose(curves['SD_p'], .2*(cols['lim_sigmaN_SD_p']/cols['sigmaN_SD_p'])**.25)
        np.testing.assert_allclose(curves['neutrino_fog_SI'], .2*(cols['neutrino_fog_SI']/cols['sigmaN_SI_eq'])**.25)
        self.assertNotIn('SI_p', curves)
        self.assertNotIn('SD_p_LZ', curves)

    def test_relic_only_and_indirect_only(self):
        for kind in ('relic','indirect'):
            cols = self.base()
            if kind == 'relic':
                cols['Omegah^2'] = np.array([.06,.12,.24])
                expected_plot = 'relic_density'
            else:
                cols['tot_SM_xsec'] = np.array([1e-26,2e-26,3e-26])
                cols['Fermi_sigmav'] = np.array([2e-26,4e-26,6e-26])
                expected_plot = 'sigmav_limits'
            with self.subTest(kind=kind):
                manifest = SUMMARY.generate(cols,self.root,kind)
                self.assertIn(expected_plot, manifest['plots'])
                path = self.root/(kind+'.txt')
                scan_file(path,cols)
                curves,meta = COUPLING.generate(path,coupling=.1)
                key = 'relic_target' if kind=='relic' else 'dwarfs'
                expected = .1*(cols['Omegah^2']/.12)**.25 if kind=='relic' else np.full(3,.1*2**.25)
                np.testing.assert_allclose(curves[key],expected)

    def test_parameter_grid_is_split_and_selection_accepts_a_coupling_axis(self):
        cols = {'run':np.array(['a','b','c','d']), 'mass#1':np.array([1.,2.,1.,2.]),
                'DMINPUTS#1':np.array([.1,.1,.2,.2]), 'Omegah2':np.array([.12,.12,.06,.06])}
        for module in (SUMMARY, COUPLING):
            axis,slices = module.scan_slices(cols,'mass#1')
            self.assertEqual(axis,'mass#1')
            self.assertEqual(len(slices),2)
            self.assertEqual(slices[0][1],'DMINPUTS#1=0.1')
            axis,selected = module.scan_slices(cols,'DMINPUTS#1',select=['mass#1=2'])
            self.assertEqual(axis,'DMINPUTS#1')
            self.assertEqual(list(selected[0][2]['run']),['b','d'])
        path=self.root/'grid.txt'
        scan_file(path,cols)
        result=COUPLING.generate(path,coupling_column='DMINPUTS#1')
        self.assertEqual(len(result),2)
        self.assertEqual(len(list(self.root.glob('*coupling_mass.pdf'))),2)
        manifest=SUMMARY.generate(cols,self.root,'grid')
        self.assertEqual(len(manifest['slices']),2)
        with self.assertRaisesRegex(ValueError,'No rows'):
            SUMMARY.scan_slices(cols,select=['mass#1=100'])

    def test_combined_analyses_and_coupling_column_use_the_same_selected_rows(self):
        cols=self.base()
        cols.update(**{'DMINPUTS#1':np.array([.1,.2,.3]),
                       'xsi':np.ones(3), 'Omegah2':np.array([.12,.1,.08]),
                       'sigmav(xf)':np.array([3e-26,3e-26,3e-26]),
                       'sigmaN_SI_p':np.array([1e-46,1e-45,1e-44]),
                       'lim_sigmaN_SI_p':np.array([1e-47,1e-46,1e-45]),
                       'tot_SM_xsec':np.array([1e-26,2e-26,3e-26]),
                       'Fermi_sigmav':np.array([2e-26,4e-26,6e-26])})
        manifest=SUMMARY.generate(cols,self.root,'combined',mass_column='DMINPUTS#1')
        for name in ('relic_density','sigmav_limits','dd_comparison'):
            self.assertIn(name,manifest['plots'])
        path=self.root/'scan.v1.txt'
        scan_file(path,cols)
        curves,meta=COUPLING.generate(path,coupling_column='DMINPUTS#1')
        np.testing.assert_allclose(curves['SI_p'],cols['DMINPUTS#1']*.1**.25)
        self.assertEqual(meta['relic_target_method'],'freezeout_rate')
        self.assertTrue((self.root/'scan.v1_coupling_mass.pdf').exists())
        self.assertEqual(meta['coupling_column'],'DMINPUTS#1')
        colors=[style['color'] for style in meta['plot_styles']]
        self.assertEqual(len(colors),len(set(colors)))

    def test_zero_or_unavailable_limits_are_not_drawn_as_constraints(self):
        cols=self.base()
        cols.update(sigma_e=np.zeros(3),lim_sigma_e=np.full(3,-1.),Omegah2=np.full(3,-1.))
        manifest=SUMMARY.generate(cols,self.root,'empty')
        self.assertEqual(manifest['plots'],{})
        path=self.root/'empty.txt'
        scan_file(path,cols)
        curves,meta=COUPLING.generate(path)
        self.assertFalse(meta['plot_created'])
        self.assertTrue(np.isnan(curves['upper_envelope']).all())
        self.assertFalse((self.root/'empty_coupling_mass.pdf').exists())


if __name__=='__main__':
    unittest.main()
