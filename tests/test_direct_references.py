"""Neutrino fog and XLZD references: rhoDM rescaling, download at generation, projected-versus-current semantics."""
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

    @staticmethod
    def fog_root(folder):
        root = Path(folder)
        floors = root/'NeutrinoFog'/'data'/'floors'; floors.mkdir(parents=True)
        for name in refs.FOG_FILES.values():
            (floors/name).write_text('# mass sigma\n1 1e-40\n100 1e-44\n')
        return root

    def test_fog_is_rescaled_to_rhoDM(self):
        with tempfile.TemporaryDirectory() as folder:
            root = self.fog_root(folder)
            native = refs.evaluate(10., rhoDM=.3, root=root)
            double = refs.evaluate(10., rhoDM=.6, root=root)
            for key in refs.FOG_KEYS:
                self.assertAlmostEqual(native[key]/1e-42, 1.)
                self.assertAlmostEqual(double[key]/native[key], .5)
            self.assertEqual(refs.evaluate(1e3, rhoDM=.6, root=root)['neutrino_fog_SI'], -1.)
            mass, sigma = refs.curves(.6, root=root)['neutrino_fog_SD_n']
            np.testing.assert_allclose(sigma, [.5e-40, .5e-44])
            lines = refs.report(double, .6)
            self.assertIn('NOT an exclusion', lines[0])
            self.assertIn('rhoDM = 0.6', lines[1])
        with tempfile.TemporaryDirectory() as folder:
            missing = refs.evaluate(10., rhoDM=.3, root=folder)
            self.assertTrue(all(missing[key] == -1. for key in refs.FOG_KEYS))
            self.assertEqual(len(missing['DD_reference_warnings']), 1)
            self.assertIn('install neutrino_fog', missing['DD_reference_warnings'][0])

    def test_fog_is_downloaded_by_generate_and_evaluated_only_for_direct_nucleon(self):
        import ast
        run = ast.parse((refs.ROOT/'maddm_run_interface.py').read_text())
        calls = lambda tree, name: [n for n in ast.walk(tree) if isinstance(n, ast.Call)
                                    and isinstance(n.func, ast.Attribute) and n.func.attr == name]
        # launch never downloads
        self.assertEqual(calls(run, 'install_neutrino_fog'), [])
        # the fog is computed in the `if self.mode['direct_nucleon']:` block only
        evaluate = [n for n in calls(run, 'evaluate') if isinstance(n.func.value, ast.Name)
                    and n.func.value.id == 'direct_references']
        self.assertEqual(len(evaluate), 1)
        blocks = [n for n in ast.walk(run) if isinstance(n, ast.If) and evaluate[0] in ast.walk(n)]
        self.assertEqual(ast.unparse(blocks[-1].test), "self.mode['direct_nucleon']")
        interface = ast.parse((refs.ROOT/'maddm_interface.py').read_text())
        generate = next(n for n in ast.walk(interface) if isinstance(n, ast.FunctionDef) and n.name == 'generate_direct')
        self.assertTrue(calls(generate, 'install_neutrino_fog'))
        launch = next(n for n in ast.walk(interface) if isinstance(n, ast.FunctionDef) and n.name == 'do_launch')
        self.assertEqual(calls(launch, 'install_neutrino_fog'), [])

    def test_xlzd_limits_come_from_rapidd_and_are_off_by_default(self):
        import ast
        from test_sigmav_channel_limits import ENV
        tree = ast.parse((refs.ROOT/'maddm_run_interface.py').read_text())
        param = next(n for n in ast.walk(tree) if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
                     and n.func.attr == 'add_param' and n.args and getattr(n.args[0], 'value', None) == 'xlzd_limits')
        self.assertEqual(param.args[1].value, 'off')
        self.assertEqual(ast.literal_eval(next(k.value for k in param.keywords if k.arg == 'allowed')),
                         ['off', '200ty', '1000ty'])
        limits = ENV['ExpConstraints']()
        limits.load_constraints()
        xlzd = refs.ROOT/'vendor'/'RAPIDD_for_DM'/'lib'/'experiments'/'XLZD'
        if not (xlzd/'XLZD_SI_200ty.dat').is_file():
            self.skipTest('RAPIDD submodule (vendor/RAPIDD_for_DM) not checked out')
        for exposure in ('200ty', '1000ty'):
            for name, value in (('SI', lambda m: limits.SI_max(m, 'XLZD_'+exposure)),
                                ('SDp', lambda m: limits.SD_max(m, 'p', 'XLZD_'+exposure)),
                                ('SDn', lambda m: limits.SD_max(m, 'n', 'XLZD_'+exposure))):
                data = np.loadtxt(xlzd/('XLZD_%s_%s.dat' % (name, exposure)))
                i = len(data)//2
                limits.set_rhoDM(.3)
                self.assertAlmostEqual(value(data[i, 0])/data[i, 1], 1., places=10)
                limits.set_rhoDM(.6)
                self.assertAlmostEqual(value(data[i, 0])/data[i, 1], .5, places=10)
                # linear extrapolation above the table, nothing above 1e15 GeV
                top = data[np.argmax(data[:, 0])]
                self.assertAlmostEqual(value(10.*top[0])/(5.*top[1]), 1., places=10)
                self.assertEqual(value(2e15), -1)
        limits.set_rhoDM(.3)
        mass, sigma, table_max = limits.dd_limit_curve('si_xlzd_200ty', 1e9)
        self.assertEqual(mass[-1], 1e9)
        self.assertEqual(table_max, np.loadtxt(xlzd/'XLZD_SI_200ty.dat')[:, 0].max())
        self.assertIsNone(limits.dd_limit_curve('missing', 1e3))

    def test_direct_detection_rows_use_equivalent_cross_sections(self):
        """LZ verdicts on the equivalent cross sections (no verdict when sigma_eq = -1), PICO-60 on sigma_SD_p,
        sigma_p and sigma_n as information only; the XLZD rows use the equivalent cross sections"""
        import ast
        import types
        from test_sigmav_channel_limits import ENV
        source = (refs.ROOT/'maddm_run_interface.py').read_text()
        tree = ast.parse(source)
        scope = dict(ENV)
        for node in tree.body:
            if isinstance(node, ast.Assign) and any(getattr(t, 'id', '').startswith('DD_') for t in node.targets):
                exec(compile(ast.Module(body=[node], type_ignores=[]), 'maddm_run_interface.py', 'exec'), scope)
        run = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'MADDMRunCmd')
        for node in run.body:
            if isinstance(node, ast.FunctionDef) and node.name in ('direct_detection_rows', 'xlzd_exposure'):
                exec(compile(ast.Module(body=[node], type_ignores=[]), 'maddm_run_interface.py', 'exec'), scope)
        res = {'sigmaN_SI_p': 1e-45, 'sigmaN_SI_n': 1e-45, 'sigmaN_SD_p': 1e-40, 'sigmaN_SD_n': 1e-41,
               'sigmaN_SI_eq': 2e-46, 'sigmaN_SD_n_eq': 3e-41, 'sigmaN_SD_p_eq': -1}
        res.update({'lim_sigmaN_' + k: 1e-44 for k in ('SI_p', 'SI_n', 'SD_p', 'SD_n', 'SI_eq', 'SD_n_eq', 'SD_p_eq',
                                                        'SI_eq_XLZD', 'SD_n_eq_XLZD', 'SD_p_eq_XLZD')})
        h = types.SimpleNamespace(maddm_card={'xlzd_limits': '200ty'}, last_results=res,
                                  limits=types.SimpleNamespace())
        scope['ExpConstraints'] = types.SimpleNamespace(XLZD_EXPOSURES=('200ty', '1000ty'))
        for name in ('direct_detection_rows', 'xlzd_exposure'):
            setattr(h, name, types.MethodType(scope[name], h))
        rows = {D['n']: D for D in h.direct_detection_rows()}
        self.assertEqual([n for n, D in rows.items() if D['verdict']],
                         ['SigmaN_SI_eq', 'SigmaN_SD_n_eq', 'SigmaN_SD_p_eq', 'SigmaN_SD_p', 'SigmaN_SI_eq_XLZD',
                          'SigmaN_SD_n_eq_XLZD', 'SigmaN_SD_p_eq_XLZD'])
        self.assertEqual(rows['SigmaN_SI_eq']['sig'], 2e-46)
        self.assertEqual(rows['SigmaN_SD_p']['exp'], 'Pico60 (2019)')
        self.assertEqual(rows['SigmaN_SD_p_eq']['lim'], -1)
        self.assertEqual(rows['SigmaN_SD_p_eq_XLZD']['lim'], -1)
        self.assertEqual(rows['SigmaN_SD_n_eq_XLZD']['sig'], 3e-41)
        self.assertFalse(rows['SigmaN_SI_p']['verdict'])

    def test_multinest_direct_nucleon_likelihood_uses_equivalent_cross_sections(self):
        """half_gauss: penalty only above the limit, on sigma_eq for LZ; SD-p also on sigma_SD_p for PICO-60 (the
        most constraining); widths from the multinest card names; gaussian with the dd_si_limits values"""
        import ast
        import types
        from test_sigmav_channel_limits import ENV
        tree = ast.parse((refs.ROOT/'maddm_run_interface.py').read_text())
        scope = dict(ENV)
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name == 'secure_float_f77' or \
                    isinstance(node, ast.ClassDef) and node.name == 'Multinest':
                exec(compile(ast.Module(body=[node], type_ignores=[]), 'maddm_run_interface.py', 'exec'), scope)
        mn = object.__new__(scope['Multinest'])
        mn.options = {'half_gauss_width': {'spinSI': 0.01, 'spinSD_p': 0.01, 'spinSD_n': 0.01, 'directSI': 0.5}}
        mn.maddm_run = types.SimpleNamespace(limits=types.SimpleNamespace(
            _sigma_SI=2e-47, _sigma_SI_width=1e-47, _sigma_SDp=-1, _sigma_SDp_width=-1, _sigma_SDn=-1,
            _sigma_SDn_width=-1))
        res = {'sigmaN_SI_eq': 1e-46, 'lim_sigmaN_SI_eq': 1e-47, 'sigmaN_SI_p': 1e-50, 'lim_sigmaN_SI_p': 1e-47,
               'sigmaN_SD_n_eq': 1e-43, 'lim_sigmaN_SD_n_eq': 1e-42,
               'sigmaN_SD_p_eq': 1e-41, 'lim_sigmaN_SD_p_eq': 1e-40, 'sigmaN_SD_p': 1e-39, 'lim_sigmaN_SD_p': 1e-41}
        # SI above the LZ limit by a factor 10, width 0.5 from the card (key directSI)
        self.assertAlmostEqual(mn.direct_nucleon_loglike('directSI', 'half_gauss', res), -0.5*(1/0.5)**2)
        # SD-n below the limit: no penalty
        self.assertEqual(mn.direct_nucleon_loglike('directSD_n', 'half_gauss', res), 0.)
        # SD-p: sigma_eq below the LZ limit, sigma_SD_p 100 times the PICO-60 limit (legacy width spinSD_p)
        self.assertAlmostEqual(mn.direct_nucleon_loglike('directSD_p', 'half_gauss', res), -0.5*(2/0.01)**2)
        # no prediction (-1) or no limit: no contribution
        self.assertEqual(mn.direct_nucleon_loglike('directSI', 'half_gauss', dict(res, sigmaN_SI_eq=-1)), 0.)
        self.assertEqual(mn.direct_nucleon_loglike('directSI', 'half_gauss', dict(res, lim_sigmaN_SI_eq=-1)), 0.)
        self.assertAlmostEqual(mn.direct_nucleon_loglike('directSI', 'gaussian', res), -0.5*(8e-47/1e-47)**2)
        self.assertEqual(mn.direct_nucleon_loglike('directSD_n', 'user', res), 0.)

    def test_run_folder_dd_limits_plot(self):
        import ast, json, types
        from test_sigmav_channel_limits import ENV
        import run_products
        tree = ast.parse((refs.ROOT/'maddm_run_interface.py').read_text())
        run = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'MADDMRunCmd')
        with tempfile.TemporaryDirectory() as folder:
            root = self.fog_root(folder)
            scope = dict(ENV, json=json, DD_EQUIVALENT_KEYS=('SI_eq', 'SD_n_eq', 'SD_p_eq'),
                         direct_references=types.SimpleNamespace(
                curves=lambda rho: refs.curves(rho, root=root), FOG_SOURCE=refs.FOG_SOURCE))
            for node in run.body:
                if isinstance(node, ast.FunctionDef) and node.name in ('write_dd_limits', 'xlzd_exposure'):
                    exec(compile(ast.Module(body=[node], type_ignores=[]), 'maddm_run_interface.py', 'exec'), scope)
            h = types.SimpleNamespace(maddm_card={'rhoDM': .6, 'xlzd_limits': '1000ty'},
                                      proc_characteristics={'dm_candidate': [1000022]},
                                      param_card=types.SimpleNamespace(get_value=lambda *a: 50.),
                                      last_results={'xsi': .5, 'sigmaN_SI_p': 1e-47, 'sigmaN_SI_n': 1.1e-47,
                                                    'sigmaN_SD_p': 1e-42, 'sigmaN_SD_n': -1.,
                                                    'sigmaN_SI_eq': 9e-48, 'sigmaN_SD_n_eq': 2e-42,
                                                    'sigmaN_SD_p_eq': 5e-41},
                                      limits=scope['ExpConstraints']())
            h.limits.load_constraints()
            for name in ('write_dd_limits', 'xlzd_exposure'):
                setattr(h, name, types.MethodType(scope[name], h))
            out = Path(folder)/'run_01'; out.mkdir()
            h.write_dd_limits(str(out))
            data = json.loads((out/'DD_limits.json').read_text())
            self.assertEqual(data['rhoDM'], .6)
            self.assertEqual(data['predictions_eq'], {'SI': 9e-48, 'SD_n': 2e-42, 'SD_p': 5e-41})
            labels = {(c['kind'], c['nucleons'], c['label']) for c in data['curves']}
            self.assertIn(('SI', 'pn', 'LZ (2024)'), labels)
            self.assertIn(('SI', 'pn', 'Xe neutrino fog'), labels)
            if (refs.ROOT/'vendor/RAPIDD_for_DM/lib/experiments/XLZD').is_dir():
                self.assertIn(('SD', 'n', 'XLZD 1000 t yr (projected)'), labels)
            fog = next(c for c in data['curves'] if c['style'] == 'fog' and c['kind'] == 'SI')
            np.testing.assert_allclose(fog['sigma'], [.5e-40, .5e-44])
            run_products.write_run_products(out, recoil=True)
            self.assertTrue((out/'DD_limits.pdf').is_file())
            self.assertIn('DD_limits.pdf', (out/'README.md').read_text())

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
                  'lim_sigmaN_SI_p_XLZD':np.full(3,1e-49),
                  'neutrino_fog_SI':np.full(3,1e-50)}
            result=SUMMARY.generate(cols,root,'reference')
            self.assertIn('neutrino_fog_SI',result['plots']['dd_comparison']['columns'])
            self.assertIn('lim_sigmaN_SI_p_XLZD',result['plots']['dd_limits']['columns'])
            self.assertEqual(SUMMARY.dd_label('SI','p','_XLZD'),'XLZD (projected); SI proton')
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
