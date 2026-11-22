"""Resolve arithmetic links between SLHA inputs at every MadDM scan point.

The ordinary MadGraph scan iterator retains its scan/scanN grouping rules.
Linked values use a restricted arithmetic interpreter, never Python eval.
"""
import ast
import math
import operator
from pathlib import Path
import re


class ParameterExpressionError(ValueError):
    """An invalid, ambiguous or nonnumeric parameter relation."""


_BINARY = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
           ast.Div: operator.truediv, ast.FloorDiv: operator.floordiv,
           ast.Mod: operator.mod, ast.Pow: operator.pow}
_UNARY = {ast.UAdd: operator.pos, ast.USub: operator.neg}
_FUNCTIONS = {name: getattr(math, name) for name in
              ('sqrt', 'exp', 'log', 'log10', 'sin', 'cos', 'tan')}
_FUNCTIONS.update(abs=abs, min=min, max=max)
_CONSTANTS = {'pi': math.pi, 'e': math.e}
_IDENTIFIER = re.compile(r'^[a-zA-Z_][a-zA-Z0-9_]*$')
_SCAN = re.compile(r'^scan\s*\d*\s*:', re.I)
_NUMBER = re.compile(r'^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?$')


def _tree(expression):
    if len(expression) > 1024:
        raise ParameterExpressionError('Parameter expressions must be shorter than 1024 characters.')
    try:
        tree = ast.parse(expression.strip(), mode='eval').body
    except (SyntaxError, ValueError) as error:
        raise ParameterExpressionError('Invalid parameter expression %r: %s' % (expression, error)) from error
    if sum(1 for _ in ast.walk(tree)) > 128:
        raise ParameterExpressionError('Parameter expression is too complex: '+expression)
    def validate(node):
        if isinstance(node, ast.Constant) and type(node.value) in (int, float):
            return
        if isinstance(node, ast.Name) and not node.id.startswith('_'):
            return
        if isinstance(node, ast.BinOp) and type(node.op) in _BINARY:
            validate(node.left); validate(node.right)
            return
        if isinstance(node, ast.UnaryOp) and type(node.op) in _UNARY:
            validate(node.operand)
            return
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and
                node.func.id.lower() in _FUNCTIONS and not node.keywords and 1 <= len(node.args) <= 8):
            for arg in node.args:
                validate(arg)
            return
        raise ParameterExpressionError('Unsupported syntax in %r. Use parameter names, numbers, '
                                       'arithmetic and the documented math functions.' % expression)
    validate(tree)
    return tree


def _names(tree):
    if isinstance(tree, ast.Call):
        return set().union(*(_names(arg) for arg in tree.args))
    if isinstance(tree, ast.Name):
        return {tree.id.lower()}
    return set().union(*(_names(child) for child in ast.iter_child_nodes(tree)))


def _identifiers(path):
    result = []
    if path is None or not Path(path).is_file():
        return result
    for line in Path(path).read_text().splitlines():
        fields = line.split('#', 1)[0].split()
        if len(fields) < 3 or not _IDENTIFIER.fullmatch(fields[-1]):
            continue
        try:
            code = tuple(int(value) for value in fields[1:-1])
        except ValueError:
            continue
        result.append(((fields[0].lower(), code), fields[-1].lower()))
    return result


def normalize_card_text(text, ident_path=None):
    """Make spaced arithmetic a single SLHA value without changing scan lists."""
    identifiers = _identifiers(ident_path)
    dimensions = {}
    for (block, code), _ in identifiers:
        dimensions.setdefault(block, set()).add(len(code))
    block = None
    lines = text.splitlines(keepends=True)
    # Infer matrix index counts from ordinary numeric rows if no ident card exists.
    for line in lines:
        fields = line.split('#', 1)[0].split()
        if not fields:
            continue
        if fields[0].lower() == 'block':
            block = fields[1].lower()
        elif fields[0].lower() in ('decay', 'xsection'):
            block = None
        elif block and len(fields)>1 and _NUMBER.fullmatch(fields[-1]):
            try:
                tuple(int(value) for value in fields[:-1])
            except ValueError:
                continue
            dimensions.setdefault(block, set()).add(len(fields)-1)
    block = None
    output = []
    for line in lines:
        data, marker, comment = line.partition('#')
        fields = data.split()
        if not fields:
            output.append(line); continue
        keyword = fields[0].lower()
        if keyword == 'block':
            block = fields[1].lower()
            output.append(line); continue
        if keyword == 'xsection':
            block = None
            output.append(line); continue
        decay = keyword == 'decay'
        if decay:
            block = None
        current = 'decay' if decay else block
        if current is None or current in ('spinfo', 'dcinfo') or len(fields) < 2:
            output.append(line); continue
        nindex = 1 if decay or current == 'mass' else min(dimensions.get(current, {1}))
        # Ident-card entries disambiguate matrix indices and integer operands.
        if not decay:
            for (name, code), _ in identifiers:
                if name == current and fields[:len(code)] == [str(v) for v in code]:
                    nindex = len(code)
                    break
        start = nindex+1 if decay else nindex
        value = ' '.join(fields[start:])
        if not value or _SCAN.match(value) or value.lower().startswith('auto'):
            output.append(line); continue
        if _NUMBER.fullmatch(value):
            normalized = value.replace('D','e').replace('d','e')
        else:
            _tree(value)  # Reject invalid syntax before removing whitespace.
            normalized = re.sub(r'\s+', '', value)
        ending = '\n' if line.endswith('\n') else ''
        if normalized != value or len(fields[start:]) > 1:
            line = ' '.join(fields[:start])+' '+normalized+(' #'+comment.rstrip('\r\n') if marker else '')+ending
        output.append(line)
    return ''.join(output)


def normalize_card_file(path, ident_path=None):
    path = Path(path)
    text = path.read_text()
    normalized = normalize_card_text(text, ident_path)
    if text != normalized:
        path.write_text(normalized)


class ParameterLinks:
    def __init__(self, card, ident_path=None):
        self.aliases, self.rules = {}, {}
        keys = set()
        def alias(name, key):
            if _IDENTIFIER.fullmatch(name):
                self.aliases.setdefault(name.lower(), set()).add(key)
        for block in card.order:
            if block.name in ('spinfo', 'dcinfo'):
                continue
            for param in block:
                key = (block.name, tuple(param.lhacode))
                keys.add(key)
                alias(block.name+'_'+'_'.join(str(i).replace('-','m') for i in key[1]), key)
                # Generated cards put the UFO name first in the comment.
                comment = (param.comment or '').strip().split()
                if comment:
                    alias(comment[0], key)
                value = param.value
                if isinstance(value, str) and not (_SCAN.match(value) or value.lower().startswith('auto')):
                    if _NUMBER.fullmatch(value):
                        continue
                    self.rules[key] = (value, _tree(value))
        for key, name in _identifiers(ident_path):
            if key in keys:
                alias(name, key)
        self.order, visiting, done = [], set(), set()
        def visit(key):
            if key in done:
                return
            if key in visiting:
                raise ParameterExpressionError('Cyclic parameter relation involving '+self.label(key))
            visiting.add(key)
            for name in _names(self.rules[key][1]):
                dependency = self.resolve(name)
                if dependency in self.rules:
                    visit(dependency)
            visiting.remove(key)
            done.add(key)
            self.order.append(key)
        for key in self.rules:
            visit(key)

    @staticmethod
    def label(key):
        return key[0]+'#'+'_'.join(map(str,key[1]))

    def resolve(self, name):
        matches = self.aliases.get(name, ())
        if len(matches)>1:
            raise ParameterExpressionError('Ambiguous parameter name %r: %s. Use a BLOCK_CODE name.' %
                                           (name, ', '.join(sorted(self.label(key) for key in matches))))
        if matches:
            return next(iter(matches))
        if name in _CONSTANTS:
            return None
        raise ParameterExpressionError('Unknown parameter name %r. Use its UFO name from ident_card.dat, '
                                       'its # comment name, or BLOCK_CODE.' % name)

    def evaluate(self, card):
        def value(node):
            if isinstance(node, ast.Constant):
                return float(node.value)
            if isinstance(node, ast.Name):
                key = self.resolve(node.id.lower())
                return _CONSTANTS[node.id.lower()] if key is None else float(card[key[0]].get(key[1]).value)
            if isinstance(node, ast.UnaryOp):
                return _UNARY[type(node.op)](value(node.operand))
            if isinstance(node, ast.BinOp):
                left, right = value(node.left), value(node.right)
                if isinstance(node.op, ast.Pow) and abs(right)>1024:
                    raise ValueError('Exponent exceeds 1024')
                return _BINARY[type(node.op)](left, right)
            return _FUNCTIONS[node.func.id.lower()](*(value(arg) for arg in node.args))
        for key in self.order:
            expression, tree = self.rules[key]
            try:
                number = float(value(tree))
                if not math.isfinite(number):
                    raise ValueError('result is not finite')
            except (ArithmeticError, TypeError, ValueError) as error:
                raise ParameterExpressionError('Cannot evaluate %s = %s: %s' %
                                               (self.label(key), expression, error)) from error
            param = card[key[0]].get(key[1])
            param.value, param.format = number, 'float'
        return [float(card[block].get(code).value) for block,code in self.order]


def make_scan_iterator(card_module, text, ident_path=None):
    """Retain MadGraph's iterator and append linked inputs to its summary columns."""
    normalized = normalize_card_text(text, ident_path)
    class LinkedParamCardIterator(card_module.ParamCardIterator):
        def __init__(self):
            super().__init__(normalized)
            self.links = ParameterLinks(self, ident_path)
        def iterate(self):
            configured = False
            for card in super().iterate():
                derived = self.links.evaluate(card)
                if not configured:
                    self.param_order.extend(self.links.label(key) for key in self.links.order)
                    configured = True
                self.itertag.extend(derived)
                yield card
    return LinkedParamCardIterator()
