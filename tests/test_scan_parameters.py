"""Linked SLHA inputs using the real MadGraph iterator when available."""
import ast
from pathlib import Path
import sys
import tempfile
import unittest

import scan_parameters as links

try:
    from madgraph.models import check_param_card as cards
except ImportError:
    try:
        from models import check_param_card as cards
    except ImportError:
        cards = None


@unittest.skipIf(cards is None, 'requires MadGraph on PYTHONPATH')
class LinkedScanTests(unittest.TestCase):
    def iterator(self, expression, extra='', values='[1,3,4,5,7,8,10,100]'):
        text=('BLOCK MASS\n 54 '+expression+' # my0\n 5000001 scan1:'+values+' # mxr\n'+extra)
        iterator=links.make_scan_iterator(cards,text)
        iterator.logging=False
        return iterator

    def test_requested_rules_and_restoration(self):
        for expression, calculation in [('2*mxr',lambda x:2*x),('mxr/2',lambda x:x/2),
                                        ('mxr**2',lambda x:x**2),('2 * (mxr + 1)',lambda x:2*(x+1))]:
            with self.subTest(expression=expression):
                iterator=self.iterator(expression)
                self.assertEqual(iterator.next(autostart=True).get_value('mass',54),calculation(1))
                self.assertEqual(iterator.param_order,['mass#5000001','mass#54'])
                rows=[(1.,calculation(1))]
                for card in iterator:
                    source=card.get_value('mass',5000001)
                    rows.append((source,card.get_value('mass',54)))
                    self.assertEqual(card.get_value('mass',54),calculation(source))
                    self.assertEqual(iterator.itertag,[source,calculation(source)])
                    self.assertNotIn('mxr',card['mass'].get(54).__str__().split('#')[0])
                self.assertEqual(len(rows),8)
                original=iterator.write()
                self.assertIn('scan1:',original)
                self.assertIn(expression.replace(' ',''),original)

    def test_chained_links_cartesian_and_correlated_scans(self):
        text='''BLOCK MASS
 54 2*mxr # my0
 5000001 scan1:[1,3] # mxr
 55 my0+g # second
BLOCK DMINPUTS
 1 scan2:[0.1,0.2] # g
 2 scan1:[4,5] # correlated
 3 sqrt(second**2) # derived
'''
        iterator=links.make_scan_iterator(cards,text)
        iterator.logging=False
        seen=[]
        for card in iterator:
            mass=card.get_value('mass',5000001)
            g=card.get_value('dminputs',1)
            self.assertAlmostEqual(card.get_value('mass',55),2*mass+g)
            self.assertAlmostEqual(card.get_value('dminputs',3),2*mass+g)
            seen.append((mass,g,card.get_value('dminputs',2)))
        self.assertEqual(seen,[(1,.1,4),(1,.2,4),(3,.1,5),(3,.2,5)])

    def test_ident_card_matrix_indices_and_fortran_numbers(self):
        with tempfile.TemporaryDirectory() as directory:
            ident=Path(directory)/'ident_card.dat'
            ident.write_text('mass 5000001 physicalmass\nyukawa 1 2 mix\n')
            text='''BLOCK MASS
 5000001 scan1:[2,4] # label
 54 physicalmass / 2 # my0
BLOCK YUKAWA
 1 2 2 * physicalmass # mix
 2 2 1.0D-2 # fixed
'''
            iterator=links.make_scan_iterator(cards,text,ident)
            iterator.logging=False
            card=iterator.next(autostart=True)
            self.assertEqual(card.get_value('mass',54),1.)
            self.assertEqual(card.get_value('yukawa',(1,2)),4.)
            self.assertEqual(card.get_value('yukawa',(2,2)),.01)
            self.assertIn('yukawa#1_2',iterator.param_order)

    def test_unknown_cycles_ambiguous_and_invalid_values_fail_early(self):
        for expr,extra,pattern in [('unknown/2','','Unknown'),('my0*2','','Cyclic'),
                                   ('other','','Unknown'),('mxr/0','','Cannot evaluate')]:
            with self.subTest(expr=expr),self.assertRaisesRegex(links.ParameterExpressionError,pattern):
                self.iterator(expr,extra).next(autostart=True)
        text='BLOCK MASS\n 1 scan1:[1] # x\n 2 2 # x\n 3 x/2 # y\n'
        with self.assertRaisesRegex(links.ParameterExpressionError,'Ambiguous'):
            links.make_scan_iterator(cards,text)
        for expression in ['mxr.__class__',"__import__('os')",'[x for x in []]','mxr[0]','True', 'lambda: 1']:
            with self.subTest(expression=expression),self.assertRaises(links.ParameterExpressionError):
                self.iterator(expression)

    def test_real_maddm_card_check_resolves_before_auto_widths(self):
        root=Path(__file__).resolve().parents[1]
        tree=ast.parse((root/'maddm_run_interface.py').read_text())
        cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='MADDMRunCmd')
        node=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='check_param_card')
        import os, re, types
        env=dict(re=re,param_card_mod=cards,scan_parameters=links,pjoin=os.path.join,
                 logger=types.SimpleNamespace(info=lambda *a:None),
                 cmd=types.SimpleNamespace(CmdShell=object))
        exec(compile(ast.Module(body=[node],type_ignores=[]),'<card check>','exec'),env)
        class Run:
            check_param_card=env['check_param_card']
            in_scan_mode=False
            def do_compute_widths(self,line):
                current=cards.ParamCard(str(self.path))
                self.width_mass=current.get_value('mass',54)
        with tempfile.TemporaryDirectory() as directory:
            run=Run();run.dir_path=directory
            path=Path(directory)/'Cards';path.mkdir()
            run.path=path/'param_card.dat'
            run.path.write_text('BLOCK MASS\n 54 2*mxr # my0\n 5000001 scan1:[1,3] # mxr\nDECAY 54 Auto # wy0\n')
            run.check_param_card(str(run.path))
            self.assertEqual(run.width_mass,2.)
            self.assertEqual(cards.ParamCard(str(run.path)).get_value('mass',54),2.)
            self.assertEqual(run.param_card_iterator.param_order,['mass#5000001','mass#54'])


if __name__=='__main__':
    unittest.main()
