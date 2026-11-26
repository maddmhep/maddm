"""Heavy-mediator validity warnings of the direct detection results: q_ref estimates, thresholds and flags
(maddm_run_interface.DDMediatorCheck) and the classification of the exchanged particles of the DD diagrams
(maddm_interface.dd_exchanged_particles). The class and the function are exec'd alone."""
import ast
import math
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(source, names, scope=None):
    tree = ast.parse((ROOT/source).read_text())
    nodes = [x for x in tree.body if isinstance(x, (ast.ClassDef, ast.FunctionDef, ast.Assign))
             and (getattr(x, 'name', None) in names
                  or (isinstance(x, ast.Assign) and getattr(x.targets[0], 'id', None) in names))]
    scope = dict(scope or {})
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(ROOT/source), 'exec'), scope)
    return scope


CHECK = load('maddm_run_interface.py', ('DDMediatorCheck',), {'math': math})['DDMediatorCheck']
GEN = load('maddm_interface.py', ('DD_MEDIATOR_KINDS', 'dd_exchanged_particles'))
CARD = {'vescape': 544., 'vMP': 238., 'vSun_r': 11.1, 'vSun_phi': 12.2, 'vSun_theta': 7.3, 'vEarth_mod': 29.8}


class QRefTests(unittest.TestCase):
    def test_v_max(self):
        self.assertAlmostEqual(CHECK.v_max(CARD) * CHECK.C_KMS, 544. + math.sqrt(11.1**2 + 250.2**2 + 7.3**2) + 29.8)

    def test_q_ref_nucleus(self):
        v = CHECK.v_max(CARD)
        m_xe = 131.293 * CHECK.AMU_GEV
        # heavy DM: recoil-energy scale sqrt(2 m_Xe 10 keV) ~ 49 MeV
        self.assertAlmostEqual(CHECK.q_ref_nucleus(1000., v), math.sqrt(2 * m_xe * 1e-5))
        self.assertAlmostEqual(CHECK.q_ref_nucleus(1000., v), 0.0494, places=3)
        # light DM: kinematic limit 2 mu v_max (mu ~ mdm)
        mu = 1. * m_xe / (1. + m_xe)
        self.assertAlmostEqual(CHECK.q_ref_nucleus(1., v), 2 * mu * v)
        # xenon is always included (LZ), and it is the largest for heavy DM
        self.assertEqual(CHECK.q_ref_nucleus(1000., v, ['f']), CHECK.q_ref_nucleus(1000., v))
        self.assertLess(math.sqrt(2 * 18.998403 * CHECK.AMU_GEV * 1e-5), CHECK.q_ref_nucleus(1000., v, ['f']))

    def test_q_ref_electron(self):
        v = CHECK.v_max(CARD)
        self.assertAlmostEqual(CHECK.q_ref_electron(1., v), 10 * 7.2973525643e-3 * 0.51099895e-3)
        self.assertAlmostEqual(CHECK.q_ref_electron(1e-3, v), 2e-3 * v)


class ThresholdTests(unittest.TestCase):
    q = 0.05

    def test_single_heavy(self):
        self.assertEqual(CHECK.check('nucleon', 100., self.q, [('Y1 (32)', 1000.)], []), (0, False, []))

    def test_single_light(self):
        flag, invalid, msg = CHECK.check('nucleon', 100., self.q, [('Y1 (32)', 0.03)], [])
        self.assertEqual((flag, invalid, len(msg)), (CHECK.FLAG_LIGHT, False, 1))
        self.assertIn('%.3g' % (1 + (0.05/0.03)**2)**2, msg[0])
        self.assertIn('light mediators are not implemented', msg[0])
        # just above / below 5 q_ref
        self.assertEqual(CHECK.check('nucleon', 100., self.q, [('Y1', 5.01 * self.q)], [])[0], 0)
        self.assertEqual(CHECK.check('nucleon', 100., self.q, [('Y1', 4.99 * self.q)], [])[0], CHECK.FLAG_LIGHT)

    def test_massless(self):
        flag, invalid, msg = CHECK.check('electron', 0.1, 3.7e-5, [('a (22)', 0.)], [])
        self.assertEqual((flag, invalid), (CHECK.FLAG_LIGHT, True))
        self.assertIn('massless', msg[0])
        self.assertIn('set to -1', msg[0])

    def test_several(self):
        flag, invalid, msg = CHECK.check('nucleon', 100., self.q, [('Z (23)', 91.19), ('Y1 (32)', 1000.)], [])
        self.assertEqual((flag, invalid, len(msg)), (CHECK.FLAG_SEVERAL, False, 1))
        self.assertIn('Z (23) M = 91.2', msg[0])
        flag, _, msg = CHECK.check('nucleon', 100., self.q, [('Z (23)', 91.19), ('Y1 (32)', 0.03)], [])
        self.assertEqual(flag, CHECK.FLAG_SEVERAL | CHECK.FLAG_LIGHT)
        self.assertEqual(len(msg), 1)
        self.assertIn('Below the threshold', msg[0])
        self.assertIn('Y1 (32)', msg[0].split('Below the threshold')[1])

    def test_su_channel(self):
        # coloured mediator coupling DM to the u quark: compressed by 5 GeV (coannihilation) is fine, 0.5 GeV is not
        self.assertEqual(CHECK.check('nucleon', 500., self.q, [], [('su (1000002)', 505., 2e-3)]), (0, False, []))
        flag, invalid, msg = CHECK.check('nucleon', 500., self.q, [], [('su (1000002)', 500.5, 2e-3)])
        self.assertEqual((flag, invalid), (CHECK.FLAG_SU_DEGENERATE, False))
        self.assertIn('M - m_DM - m_f = 0.498', msg[0])
        # heavy quark: the mass of the quark counts
        self.assertEqual(CHECK.check('nucleon', 500., self.q, [], [('st (1000006)', 673., 172.7)])[0],
                         CHECK.FLAG_SU_DEGENERATE)
        # s/u-channel mediators do not count as several t-channel mediators and have no q warning
        self.assertEqual(CHECK.check('nucleon', 0.5, self.q, [('Y1', 1000.)], [('su', 2., 2e-3)])[0], 0)
        # electrons: atomic scale
        self.assertEqual(CHECK.check('electron', 0.1, 3.7e-5, [], [('E (9000010)', 0.1 + 0.000511 + 5e-3, 0.000511)])[0], 0)
        self.assertEqual(CHECK.check('electron', 0.1, 3.7e-5, [], [('E (9000010)', 0.1 + 0.000511 + 5e-4, 0.000511)])[0],
                         CHECK.FLAG_SU_DEGENERATE)


class Leg(dict):
    pass


def leg(number, pdg):
    return Leg(number=number, id=pdg)


class Obj(dict):
    pass


def amplitude(fermion, diagrams):
    """DM(1) f(2) > DM(3) f(4): diagrams as lists of vertices, each a list of (number, pdg) legs"""
    process = Obj(legs=[leg(1, 52), leg(2, fermion), leg(3, 52), leg(4, fermion)])
    return Obj(process=process,
               diagrams=[Obj(vertices=[Obj(legs=[leg(n, p) for n, p in v]) for v in d]) for d in diagrams])


class ExchangedParticlesTests(unittest.TestCase):
    def test_channels(self):
        f = GEN['dd_exchanged_particles']
        t = [[(1, 52), (3, 52), (1, 32)], [(1, 32), (2, 2), (4, 2)]]           # DM DM Y1, Y1 u u
        t2 = [[(2, 2), (4, 2), (2, 23)], [(2, 23), (1, 52), (3, 52)]]           # t channel joined on the quark line
        s = [[(1, 52), (2, 2), (1, 1000002)], [(1, 1000002), (3, 52), (4, 2)]]  # s channel
        u = [[(1, 52), (4, 2), (1, -1000002)], [(1, -1000002), (2, 2), (3, 52)]]  # u channel
        contact = [[(1, 52), (2, 2), (3, 52), (4, 2)]]
        out = f([amplitude(2, [t, t2, s, u, contact]), amplitude(11, [t])], GEN['DD_MEDIATOR_KINDS']['nucleon'])
        self.assertEqual(out, {'t': {(32, 2), (23, 2)}, 'su': {(1000002, 2)}})
        out = f([amplitude(2, [t]), amplitude(-11, [t, s])], GEN['DD_MEDIATOR_KINDS']['electron'])
        self.assertEqual(out, {'t': {(32, 11)}, 'su': {(1000002, 11)}})
        self.assertEqual(f([], (11,)), {'t': set(), 'su': set()})


if __name__ == '__main__':
    unittest.main()
