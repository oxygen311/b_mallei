from ete3 import Tree, TreeStyle

from glob import glob

tree_file = glob('*_tree.treefile')[0]
t = Tree(tree_file)


circular_style = TreeStyle()
circular_style.mode = "c" # draw tree in circular mode
circular_style.show_branch_support = True
# circular_style.scale = 20

# t.set_outgroup(t.get_common_ancestor(t & 'LEPN.0521.00010', t & 'LEPN.0521.00006'))
# t.set_outgroup(t.get_common_ancestor(t & 'PSAE.0521.00253', t & 'PSAE.0521.00085'))
# t.set_outgroup(t.get_common_ancestor(t & 'PSAE.0521.00253', t & 'PSAE.0521.00148'))
# t.set_outgroup(t.get_common_ancestor(t & 'BUMA.1121.00035', t & 'BUMA.1121.00003'))
t.set_outgroup(t & 'BUMA.1121.00068')

t.render("test.pdf", tree_style=circular_style)
## RIGHT ROOT
t.write(outfile=tree_file.replace('_tree', '_tree_right_root'))
