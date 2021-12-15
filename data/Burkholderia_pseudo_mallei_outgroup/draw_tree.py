from parebrick.tree.tree_holder import TreeHolder
from parebrick.utils.data.parsers import make_labels_dict

from collections import defaultdict

import logging

labels_file = 'labels.csv'
tree_file = '6-tree_module/BUMA_113.nucl.grp.aln.iqtree_tree.treefile'

tree_holder = TreeHolder(tree_file, logging.getLogger(), labels_dict=make_labels_dict(labels_file))
tree_holder.count_innovations_fitch(defaultdict(int))

tree_holder.draw(tree_file.replace('.treefile', '.pdf'), colors=['white'])