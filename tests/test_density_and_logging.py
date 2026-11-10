"""Cross-analysis density scaling, user data and plot logging regressions."""
import io
import json
import logging
from pathlib import Path
import tempfile
import unittest
import numpy as np
from indirect import cr_config, dbar_interface, gce_interface, pbar_interface
from plotting import plot_ams02_pbar_fluxes as pbar_plot
from plotting import signal_plots
from test_pbar_integration import record,Harness

from flux_test_data import requires_dbar


class ConsistencyTests(unittest.TestCase):
    @requires_dbar
    def test_common_density_controls_gce_dbar_and_channel_bound(self):
        rows=[record([5,-5],3.e-26)]
        engine=gce_interface.GCEFit();db=dbar_interface.AntideuteronEvents()
        r1=engine.evaluate(50.,rows,1000022,4,{'rholocal':.4})
        r2=engine.evaluate(50.,rows,1000022,4,{'rholocal':.8,'gce_rho_sun':.345})
        self.assertAlmostEqual(r2['sigmav_best_cm3s']/r1['sigmav_best_cm3s'],.25,places=9)
        self.assertAlmostEqual(r1['chi2'],r2['chi2'],places=9)
        np.testing.assert_allclose(r2['_spectrum'][:,1],4*r1['_spectrum'][:,1])
        e1=db.evaluate(50.,rows,1000022,4,1.,{'rholocal':.4})
        e2=db.evaluate(50.,rows,1000022,4,1.,{'rholocal':.8,'pbar_rho_sun':.1})
        self.assertAlmostEqual(e2['cases']['all_dm']['AMS02_total_events']/e1['cases']['all_dm']['AMS02_total_events'],4.,places=9)
        table=pbar_interface.SingleChannelLimits()
        a=pbar_interface.evaluate_channel_limits(50.,rows,1000022,4,1.,{'rholocal':.4},table.upper_limit)
        b=pbar_interface.evaluate_channel_limits(50.,rows,1000022,4,1.,{'rholocal':.8},table.upper_limit)
        self.assertAlmostEqual(b['channels'][0]['sigmav95_cm3s']/a['channels'][0]['sigmav95_cm3s'],.25,places=9)

    def test_density_validation_and_supplied_secondary_invariance(self):
        for v in [0.,-1.,float('nan'),float('inf')]:
            with self.assertRaises(ValueError):cr_config.local_density({'rholocal':v})
        for species in ['pbar','dbar']:
            for scheme in ['BIG','QUAINT','SLIM']:
                paths=[cr_config.secondary_path(dict(cr_propagation_scheme=scheme,cr_propagation_set=s),species) for s in ['MIN','MED','MAX']]
                self.assertEqual(len(set(paths)),1)
                a=np.loadtxt(paths[0]);self.assertTrue(np.all(a[:,1]>0))

    def test_summary_names_and_one_sided_statistic_explained(self):
        fixture=Path(__file__).resolve().parents[1]/'tests/fixtures/pbar/AMS02_pbar_limits.json'
        result=json.loads(fixture.read_text())
        values=pbar_interface.observables(result,.00075)
        for new,old in [('AMS02_pbar_sigmav_UL95','AMS02_pbar_sigmav'),('AMS02_pbar_sigmav_best_fit','AMS02_pbar_sigmav_best'),('AMS02_pbar_q_upper_model','AMS02_pbar_q')]:
            self.assertEqual(values[new],values[old])
        with tempfile.TemporaryDirectory() as tmp:
            h=Harness(tmp);h.last_results.update(values);h._pbar_result=result
            out=io.StringIO();h.save_pbar_summary(out);text=out.getvalue()
        self.assertIn('AMS02_pbar_sigmav_UL95 =',text)
        self.assertNotIn('\nAMS02_pbar_sigmav =',text)
        self.assertIn('not a p-value',text)
        self.assertIn('xsi=1',text)

    def test_plot_logging_preserves_warnings_and_other_debug(self):
        stream=io.StringIO();handler=logging.StreamHandler(stream)
        root=logging.getLogger();old=root.level;root.setLevel(logging.DEBUG);root.addHandler(handler)
        try:
            for fn in [signal_plots.quiet_plot_logging,pbar_plot.quiet_plot_logging]:
                logging.getLogger('fontTools.subset').setLevel(logging.DEBUG)
                logging.getLogger('matplotlib.ticker').setLevel(logging.DEBUG)
                fn()
                logging.getLogger('fontTools.subset').debug('debug noise')
                logging.getLogger('matplotlib.ticker').debug('debug noise')
            logging.getLogger('fontTools.subset').warning('real warning')
            logging.getLogger('maddm_regression').debug('physics debug')
        finally:root.removeHandler(handler);root.setLevel(old)
        self.assertNotIn('debug noise',stream.getvalue())
        self.assertIn('real warning',stream.getvalue());self.assertIn('physics debug',stream.getvalue())

    def test_gce_saved_statistical_errors_and_argonne_defaults(self):
        result=gce_interface.GCEFit().evaluate(50.,[record([5,-5],3.e-26)],1000022,4,{})
        self.assertEqual(result['_bins'].shape,(26,11))
        _,f,used=gce_interface.read_data();np.testing.assert_array_equal(result['_bins'][:,10],f[:,1])
        p=Path(__file__).resolve().parents[1]
        for n in ['pythia8_card.dat','pythia8_card_default.dat']:
            self.assertIn('Main:methodDbar=4.0000000000e+00',(p/'Templates/Cards'/n).read_text())

if __name__=='__main__':unittest.main()
