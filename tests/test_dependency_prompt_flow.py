"""Installation must be offered before Pythia8, on the configured interface."""
import ast
from pathlib import Path
import types
import unittest
from unittest.mock import patch

import data_install

ROOT=Path(__file__).resolve().parents[1]
TREE=ast.parse((ROOT/'maddm_interface.py').read_text())
CLASS=next(n for n in TREE.body if isinstance(n,ast.ClassDef) and n.name=='MadDM_interface')
NAMES=('prepare_python_requirements','install_indirect')
ENV=dict(data_install=data_install,MG5DIR=str(ROOT),
         logger=types.SimpleNamespace(info=lambda *a:None,warning=lambda *a:None),
         AskMadDMInstaller=object,MadGraph5Error=RuntimeError)
for node in CLASS.body:
    if isinstance(node,ast.FunctionDef) and node.name in NAMES:
        exec(compile(ast.Module(body=[node],type_ignores=[]),'<main interface>','exec'),ENV)


class PromptTests(unittest.TestCase):
    def harness(self, answer):
        class Interface:
            prepare_python_requirements=ENV['prepare_python_requirements']
            install_indirect=ENV['install_indirect']
            InvalidCmd=RuntimeError
            options=dict(maddm_first_indirect=True)
            def ask(self, text, default, choices=None, **kwargs):
                events.append(('ask',text,kwargs))
                if 'required packages through' in text:
                    return answer
                return {'pythia8':'install','cosmixs':'install'}
            def do_install(self,name):events.append(('install',name))
            def exec_cmd(self,line):events.append(('command',line))
        events=[]
        return Interface(),events

    def test_acceptance_checks_before_source_menu_and_pythia(self):
        interface,events=self.harness('y')
        with patch.object(data_install,'missing_python_requirements',side_effect=[['matplotlib'],[]]), \
             patch.object(data_install,'python_install_command',return_value=('conda',['conda'])), \
             patch.object(data_install,'install_python_requirements',side_effect=lambda *a,**kw:events.append(('python','install'))):
            interface.install_indirect()
        self.assertEqual(events[0][0],'ask')
        self.assertIn('through conda [y]',events[0][1])
        self.assertIs(events[0][2]['timeout'],False)
        self.assertLess(events.index(('python','install')),events.index(('install','pythia8')))

    def test_manual_choice_skips_pythia_and_preserves_first_use(self):
        interface,events=self.harness('n')
        with patch.object(data_install,'missing_python_requirements',return_value=['matplotlib']), \
             patch.object(data_install,'python_install_command',return_value=('pip',['pip'])):
            interface.install_indirect()
        self.assertFalse(any(event[0] in ('install','command') for event in events))
        self.assertTrue(interface.options['maddm_first_indirect'])

    def test_explicit_pythia_target_has_a_preflight_and_constructor_never_asks(self):
        source=(ROOT/'maddm_interface.py').read_text()
        block=source[source.index("if package == 'pythia8':"):source.index("if not shutil.which('rsync')")]
        self.assertIn('prepare_python_requirements(build_tools=True)',block)
        run=ast.parse((ROOT/'maddm_run_interface.py').read_text())
        cls=next(n for n in run.body if isinstance(n,ast.ClassDef) and n.name=='MADDMRunCmd')
        init=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='__init__')
        ensure=next(n for n in ast.walk(init) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='ensure_python_requirements')
        ask=next(k.value for k in ensure.keywords if k.arg=='ask')
        self.assertIsNone(ask.value)


if __name__=='__main__':unittest.main()
