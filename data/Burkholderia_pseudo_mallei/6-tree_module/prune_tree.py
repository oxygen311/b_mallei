from ete3 import Tree, TreeStyle

from glob import glob
import pandas as pd

tree_file = glob('*_tree_right_root.treefile')[0]
t = Tree(tree_file)

circular_style = TreeStyle()
circular_style.mode = "c" # draw tree in circular mode
circular_style.show_branch_support = True
# circular_style.scale = 20


labels = pd.read_csv('../labels.csv')
mallei = [strain for strain, label in zip(labels['strain'], labels['label']) if 'pseudomallei' not in label]

t.prune(mallei)

t.render("test.pdf", tree_style=circular_style)
## RIGHT ROOT
# t.write(outfile=tree_file.replace('_tree', '_tree_right_root'))
