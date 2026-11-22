"""Reference coverage, first-use download, and projected-versus-current semantics."""
import io
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import data_install
import direct_references as refs
from indirect import pbar_interface, gce_interface
from test_pbar_integration import record
from test_scan_plotters import SUMMARY, COUPLING, scan_file


class ReferenceTests(unittest.TestCase):
    def test_log_interpolation_and_no_extrapolation(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'table.txt'
            path.write_text('# mass sigma\n1 1e-40\n100 1e-44\n')
            self.assertAlmostEqual(refs.interpolate(path, 10)/1.e-42, 1.)
            for mass in (0., .99, 101., float('nan')):
                self.assertEqual(refs.interpolate(path, mass), -1.)
            self.assertEqual(refs.interpolate(path, 100), 1.e-44)

    def test_xlzd_is_opt_in_and_recovers_supplied_nodes(self):
        default = refs.evaluate(50.)
        self.assertFalse(any(k in default for k in refs.XLZD_KEYS))
        for exposure, key in zip((200,1000), refs.XLZD_KEYS):
            data = np.loadtxt(refs.ROOT/'ExpData/XLZD'/('XLZD_SI_90CL_%sty_digitized_smoothed.txt' % exposure))
            i = len(data)//2
            self.assertAlmostEqual(refs.evaluate(data[i,0], True)[key]/data[i,1], 1., places=12)
        self.assertIn('NOT a current exclusion', '\n'.join(refs.report(refs.evaluate(50., True))))

    def test_launch_reference_check_with_actual_madgraph_card(self):
        import ast
        import types
        try:
            from PLUGIN.maddm.maddm_run_interface import MadDMCard
        except ImportError:
            self.skipTest('requires MadGraph plugin imports on PYTHONPATH')
        tree = ast.parse((refs.ROOT/'maddm_run_interface.py').read_text())
        launch = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == 'do_launch'
                      and any(isinstance(x, ast.Name) and x.id == 'direct_references' for x in ast.walk(n)))
        statement = next(n for n in ast.walk(launch) if isinstance(n, ast.Expr)
                         and any(isinstance(x, ast.Attribute) and x.attr == 'evaluate'
                                 and isinstance(x.value, ast.Name) and x.value.id == 'direct_references'
                                 for x in ast.walk(n)))
        code = compile(ast.Module(body=[statement], type_ignores=[]), '<actual DD reference launch>', 'exec')
        for flag in (None, False, True):
            with self.subTest(xlzd=flag), tempfile.TemporaryDirectory() as folder:
                path = Path(folder)/'maddm_card.dat'
                # A card from before XLZD was added must also load the default.
                path.write_text('True = relic_canonical\n' if flag is None else str(flag)+' = xlzd_limits\n')
                card = MadDMCard(str(path))
                values = {}
                exec(code, {'self':types.SimpleNamespace(maddm_card=card),
                            'direct_references':refs, 'mdm':100., 'result':values})
                self.assertEqual(any(key in values for key in refs.XLZD_KEYS), flag is True)
                if flag is True:
                    self.assertTrue(all(values[key] > 0 for key in refs.XLZD_KEYS))

    def test_first_use_download_is_atomic_and_cached(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)/'plugin'; root.mkdir()
            archive = Path(folder)/'upstream.tar.gz'
            prefix = 'NeutrinoFog-'+data_install.NEUTRINO_FOG_REVISION
            with tarfile.open(archive, 'w:gz') as tar:
                for name in ('data/floors/Xe_SI.txt','data/floors/Xe_SDp.txt','data/floors/Xe_SDn.txt','LICENSE'):
                    raw = b'1 1e-40\n100 1e-44\n' if name.endswith('txt') else b'MIT'
                    member = tarfile.TarInfo(prefix+'/'+name); member.size=len(raw)
                    tar.addfile(member, io.BytesIO(raw))
            with patch.object(data_install, 'verified_archive', return_value=archive) as download:
                destination=data_install.install_neutrino_fog(root)
                data_install.install_neutrino_fog(root)
                self.assertEqual(download.call_count,1)
                self.assertTrue((destination/'LICENSE').is_file())
                self.assertAlmostEqual(refs.evaluate(10., root=root)['neutrino_fog_SI']/1e-42,1.)
                self.assertFalse(list(root.glob('.neutrino-fog-*')))

    def test_missing_mediator_contribution_preserves_gce_signal_rate(self):
        engine=gce_interface.GCEFit()
        known=[record([5,-5],3e-26)]
        complete=engine.evaluate(50.,known,1000022,4,{})
        partial=engine.evaluate(50.,known+[record([9000001,9000001],9e-26)],1000022,4,{})
        for key in ('sigmav_best_cm3s','chi2','model_sigmav_cm3s'):
            self.assertEqual(partial[key],complete[key])
        self.assertAlmostEqual(partial['omitted_fraction'],.75)
        self.assertTrue(partial['warnings'])
        output={};pbar_interface.record_omissions(output,'GCE',partial)
        self.assertAlmostEqual(output['GCE_omitted_sigmav']/1e-26,9.)
        self.assertIn('unknown source spectra', output['ID_incomplete_warnings']['GCE'][0])

    def test_new_references_are_plotted_but_not_in_current_exclusion_envelope(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            cols={'run':np.array(['a','b','c']),'mass#1':np.array([10.,50.,100.]),
                  'sigmaN_SI_p':np.full(3,1e-45),'lim_sigmaN_SI_p':np.full(3,1e-46),
                  'lim_sigmaN_SI_XLZD_200ty':np.full(3,1e-49),
                  'lim_sigmaN_SI_XLZD_1000ty':np.full(3,5e-50),
                  'neutrino_fog_SI':np.full(3,1e-50)}
            result=SUMMARY.generate(cols,root,'reference')
            self.assertIn('neutrino_fog_SI',result['plots']['dd_comparison']['columns'])
            self.assertIn('lim_sigmaN_SI_XLZD_1000ty',result['plots']['dd_limits']['columns'])
            path=root/'scan.txt';scan_file(path,cols)
            curves,meta=COUPLING.generate(path,coupling=.1)
            self.assertIn('neutrino_fog_SI',curves)
            self.assertEqual(meta['envelope_inputs'],['SI_p'])
            np.testing.assert_allclose(curves['upper_envelope'],.1*(.1**.25))


class PartialTemplateTests(unittest.TestCase):
    def test_antiproton_antideuteron_and_earth_templates_keep_known_amplitude(self):
        from types import SimpleNamespace
        from indirect import cr_config, dbar_interface, cr_earth
        from indirect.pbar.core import CHANNELS, FluxTable
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); table=root/'synthetic.dat'
            energy=np.geomspace(.001,1e5,100)
            # Controlled template fixture, not a physical transport prediction.
            rows=[[mass,np.log10(e),*([1e-5*np.exp(-e/100.)]*len(CHANNELS))]
                  for mass in (10.,100.,1000.) for e in energy]
            np.savetxt(table, rows, header='mDM log10K '+' '.join(CHANNELS))
            secondary=root/'secondary.dat'
            np.savetxt(secondary,np.column_stack((energy,1e-7/(1+energy)**2)))
            known=[record([5,-5],3e-26)]; extra=known+[record([9999,9999],9e-26)]
            a=pbar_interface.channel_mixture(known,1000022,4)
            b=pbar_interface.channel_mixture(extra,1000022,4)
            with patch.object(cr_config,'table_path',return_value=table), \
                 patch.object(cr_config,'secondary_path',return_value=secondary):
                d=dbar_interface.AntideuteronEvents()
                x=d.evaluate(100.,known,1000022,4,1.,{})
                y=d.evaluate(100.,extra,1000022,4,1.,{})
                self.assertEqual(x['cases'],y['cases'])
                self.assertTrue(y['warnings'])
                for species in ('pbar','dbar'):
                    earth=cr_earth.EarthFluxes()
                    x=earth.evaluate(species,100.,a,4,1.,{})
                    y=earth.evaluate(species,100.,b,4,1.,{})
                    np.testing.assert_array_equal(x['_columns'],y['_columns'])
                    self.assertEqual(y['prediction_scope'],'supported_channels_only')
            class Likelihood:
                secondary_path=secondary
                config=SimpleNamespace(delta_chi2=2.71)
                def upper_limit(self, template):
                    return dict(status='ok',sigmav95_cm3s=1e-25,
                                best_fit=dict(chi2=0.,sigmav_cm3s=0.),warnings=[])
                def profile(self,template,rate):
                    value=(rate/1e-26)**2
                    return SimpleNamespace(chi2=value,to_dict=lambda:dict(chi2=value))
            engine=pbar_interface.AntiprotonLimits();engine.table=FluxTable(table)
            with patch.object(engine,'configure',return_value=(Likelihood(),('fixture',))):
                x=engine.evaluate(100.,known,1000022,4,1.,{})
                y=engine.evaluate(100.,extra,1000022,4,1.,{})
                self.assertEqual(x['model_points'],y['model_points'])
                self.assertEqual(y['total_sigmav_cm3s'],3e-26)
                self.assertTrue(y['warnings'][0].startswith('Incomplete indirect'))

    def test_prompt_source_and_dwarf_likelihood_use_supported_rate(self):
        import ast,types
        from test_sigmav_channel_limits import RUN, ENV
        import auxiliary
        scope=dict(ENV,HAS_SCIPY=True,aux=auxiliary)
        names=('read_tabulatedspectra','cont_spectra_available_channels','calculate_fermi_limits',
               'dNdE_dPhidE','dNdlog10x_int','Phi')
        for node in RUN.body:
            if isinstance(node,ast.FunctionDef) and node.name in names:
                exec(compile(ast.Module(body=[node],type_ignores=[]),'<actual run methods>','exec'),scope)
        h=types.SimpleNamespace()
        for name in names:
            setattr(h,name,types.MethodType(scope[name],h))
        h.maddm_card={'indirect_flux_source_method':'cosmixs','sigmav_method':'reshuffling',
                      'vave_indirect_cont':2e-5,'fermi_fast_dwarfs':'all'}
        h.mode={'indirect':'flux_source'};h.vave_indirect_cont_range=(3e-6,1.5e-4)
        h.options={'pppc4dmid_path':'fixture'};h.proc_characteristics={'dm_candidate':[1000022]}
        h.param_card=types.SimpleNamespace(get_value=lambda *a:100.)
        h.is_spectral_finalstate=lambda p:False;h.norm_Majorana_Dirac=lambda:4
        h.str_processes={'known':types.SimpleNamespace(final_states=['5','-5']),
                         'dark':types.SimpleNamespace(final_states=['9999','9999'])}
        h.last_results={'taacsID#known':3e-26,'taacsID#dark':9e-26,'taacsID':12e-26,'xsi':1.}
        spectra={'x':[],'gammas':[]}
        h.Spectra=types.SimpleNamespace(spectra=spectra,spectra_id={'gx':'gammas'},
            map_allowed_final_state_PPPC={'5.-5':'bb'},initialize_spectra=lambda:spectra.update(x=[],gammas=[]),
            check_mass=lambda m:True,load_tabulated_sourcespectra=lambda *a,**k:{'x':[.1,.3,.5]},
            interpolate_spectra=lambda *a,**k:np.array([1.,2.,3.]))
        self.assertTrue(h.read_tabulatedspectra())
        np.testing.assert_array_equal(spectra['gammas'],[1.,2.,3.])
        self.assertEqual(h.last_results['tot_SM_xsec'],3e-26)
        self.assertAlmostEqual(h.last_results['Fermi_omitted_fraction'],.75)
        evaluated=[]
        def limit(*args,**kw):
            if 'sigmav_th' in kw:
                evaluated.append(kw['sigmav_th']);return .5,1.
            return 1e-25
        h.Fermi=types.SimpleNamespace(Fermi_sigmav_lim=limit,last_limit={'n_dwarfs':1})
        with tempfile.TemporaryDirectory() as folder:
            h.dir_path=folder;h.run_name='point';(Path(folder)/'output/point').mkdir(parents=True)
            h.calculate_fermi_limits(100.)
        self.assertEqual(evaluated,[3e-26])
        self.assertEqual(h.last_results['Fermi_model_sigmav'],3e-26)

        h.maddm_card.update(halo_profile='NFW', jfactors={'NFW':1.}, npts_for_flux=50)
        h.Spectra.flux_source={'e':100.*np.asarray(spectra['x']), 'gammas':{}}
        h.Spectra.flux_earth={'gammas':{}}
        h.dNdE_dPhidE('gammas')
        differential=h.Spectra.flux_earth['gammas']['dPhidE'].copy()
        integrated=h.Phi('gammas')
        # Removing the unmodeled contribution from the full total must not
        # change either flux: the original supported rate was already used.
        h.last_results['taacsID']=3e-26
        h.dNdE_dPhidE('gammas')
        np.testing.assert_array_equal(differential,h.Spectra.flux_earth['gammas']['dPhidE'])
        self.assertEqual(integrated,h.Phi('gammas'))
