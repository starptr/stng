from stencil_ir import *
import asp.codegen.ast_tools as ast_tools
import os, binascii
import logging
from backend_halide import IndexResugar

class ToZ3(ast_tools.NodeVisitor):
    """
    Converts a set of assertions into Z3 syntax.
    """
    def __init__(self, tree, loopvars=None, additional_conds=None, convert_floats=False, funcdict={}, inputs=[]):
        self.tree = tree
        self.convert_floats = convert_floats
        if loopvars:
            self.loopvars = loopvars
        else:
            self.loopvars = ["i"]
        if additional_conds:
            self.additional_conds = additional_conds
        else:
            self.additional_conds = []

        self.funcdict = funcdict
        logging.debug("Func dict is %s", funcdict)
        self.inputs = inputs

    def to_str(self):
        return self.visit(self.tree) + '\n'.join([self.visit(x) for x in self.additional_conds])

    def visit_NumNode(self, node):
        return str(node.val)

    def visit_VarNode(self, node):
        if self.convert_floats and "__float__" in node.name:
            import re
            num = re.sub("__float__", "", node.name)
            num = re.sub("_", ".", num)
            logging.debug("Converting %s to %s", node.name, float(num))
            return str(float(num))
        return node.name

    def visit_BinExp(self, node):
        translate_op = {"&&":"and", "||":"or"}
        if node.op in translate_op.keys():
            op = translate_op[node.op]
        else:
            op = node.op
        return "(%s %s %s)" % (op, self.visit(node.left), self.visit(node.right))

    def visit_CallExp(self, node):
        return "(%s %s)" % (self.visit(node.fname), ' '.join([self.visit(x) for x in node.params]))

    def visit_NotExp(self, node):
        return "(not %s)" % (self.visit(node.ex))

    def visit_ArrExp(self, node):
        def in2pre(m):
            """Only handles t +|-|*|/ n"""
            import re
            splitted = re.split('(\+|\-|\*|\/)', m)
            if len(splitted) > 1:
                return "(%s %s %s)" % (splitted[1].strip(), splitted[0].strip(), splitted[2].strip())
            else:
                return m
        def construct_select(terms):
            if len(terms) == 2:
                return "(select %s %s)" % (terms[0], terms[1])
            else:
                return "(select %s %s)" % (construct_select(terms[:-1]), terms[-1])

        logging.debug("inputs: %s", self.inputs)
        logging.debug("funcdict %s", self.funcdict)
        logging.debug("loopvars %s", self.loopvars)
        logging.debug("node: %s", tree_to_str(node))
        idx_expression,self.should_reverse = IndexResugar(self.funcdict, self.loopvars, self.inputs).resugar(node.name.name, node.loc)
        logging.debug("idx expression is: %s", idx_expression)
        return construct_select([self.visit(node.name)] + map(in2pre, idx_expression.split(',')))

    def visit_AugArrayVarNode(self, node):
        # helper to turn into nested store/select statements
        def construct_select(terms):
            if len(terms) == 2:
                return "(select %s %s)" % (terms[0], terms[1])
            else:
                return "(select %s %s)" % (construct_select(terms[:-1]), terms[-1])


        def nester(m, aug):
            if len(m)==2:
                return "(store %s %s %s)" % (m[0], m[1], aug)
            else:
                logging.debug("Stuff to concat: %s and %s", m[:2], m[2:])
                return "(store %s %s %s)" % (m[0], m[1], nester([construct_select(m[:2])]+m[2:], aug))

        # helper to recursively add augmentations
        def helper(node, aug_keys, cur_idx):
            if cur_idx == len(aug_keys)-1:
                name = self.visit(node.name)
                idx_expression,should_reverse = IndexResugar(self.funcdict, self.loopvars, self.inputs).resugar(name, aug_keys[cur_idx])
                return nester([name]+idx_expression.split(','), self.visit(node.augmentation[aug_keys[cur_idx]]))
            else:
                name = helper(node, aug_keys, cur_idx+1)
            a = self.visit(aug_keys[cur_idx])
            b = self.visit(node.augmentation[aug_keys[cur_idx]])
            return "(store %s %s %s)" % (name, a, b)
        return helper(node, node.augmentation.keys(), 0)

    def visit_ImplicationExp(self, node):
        # we need to change loopvar to loopvar_valp for the postcondition clause,
        # due to the way the Z3 script is being structured
        class PostconFinder(ast_tools.NodeVisitor):
            def __init__(self):
                self.found = False
            def find(self, n):
                self.visit(n)
                return self.found
            def visit_CallExp(self, n):
                if n.fname.name == "postcondition":
                    self.found = True

        ifx = node.ifx

        return "(=> %s %s)\n" % (self.visit(ifx), self.visit(node.then))
