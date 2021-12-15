from ete3 import Tree, TreeStyle, TextFace, RectFace, NodeStyle, CircleFace, faces, PieChartFace

from collections import defaultdict, Counter
from itertools import combinations
from matplotlib.colors import LinearSegmentedColormap

import math
import matplotlib.cm

import numpy as np

from parebrick.tree.neighbours_utils import generate_neighbour_face, align_neighbours, get_offsets


class TreeHolder:
    def __init__(self, tree, logger, scale=None, labels_dict=None, node_colors=defaultdict(lambda: 'black')):
        self.tree = Tree(tree)
        self.scale = scale

        for node in self.tree.traverse():
            if len(node.children) == 3:
                logger.info("Trying to root tree by first child of root")
                logger.info(f'Children of root: {node.children}')
                self.tree.set_outgroup(node.children[0])
            break

        for node in self.tree.traverse():
            # Hide node circles
            node.img_style['size'] = 0

            if node.is_leaf():
                try:
                    name_face = TextFace(labels_dict[node.name] if labels_dict else node.name,
                                         fgcolor=node_colors[node.name])
                except KeyError:
                    msg = f'There is not label for leaf {node.name} in labels file'
                    logger.error(msg)
                    raise KeyError(msg)
                node.add_face(name_face, column=1)

    def draw(self, file, block_to_length, show_branch_support=True, show_scale=True, legend_scale=1,
             mode="c"):
        def get_genome_sizes(node):
            return sum(block_to_length[int(block)] * c for block, c in node.blocks1.items()),\
                   sum(block_to_length[int(block)] * c for block, c in node.blocks2.items())

        # def get_genome_size(node):
        #     return sum(block_to_length[int(block)] * c for block, c in node.blocks.items())

        def calculate_diffs(v, parent):
            l_before1, l_before2 = get_genome_sizes(v)
            l_after1, l_after2 = get_genome_sizes(parent)
            diff1 = l_after1 - l_before1
            diff2 = l_after2 - l_before2

            # return math.sqrt(diff) * scale
            return diff1, diff2

        # def calculate_diffs(v, parent):
        #     l_before = get_genome_size(v)
        #     l_after = get_genome_size(parent)
        #     diff = l_after - l_before
        #
        #     # return math.sqrt(diff) * scale
        #     return diff

        def color_nodes(v, parent):
            nonlocal max_stat
            nstyle = NodeStyle()
            nstyle["size"] = 0
            v.set_style(nstyle)

            diffs = np.array(calculate_diffs(v, parent))
            diffs = np.max([diffs, [0, 0]], axis=0)
            # diff = calculate_diffs(v, parent)
            # diff = max(calculate_diffs(v, parent), 0)

            branch_length = self.tree.get_distance(v, parent)
            scale = 1e-4

            stat = math.sqrt(sum(diffs) * scale)
            # stat = math.sqrt(diff * scale)

            max_stat = max(stat, max_stat)

            # v.add_face(CircleFace(stat, 'blue'), column=0)
            v.add_face(PieChartFace(diffs / sum(diffs) * 100 if sum(diffs) > 0 else [50, 50], stat, stat), column=0)

            for child in v.children:
                color_nodes(child, v)

        # assign colors
        root = self.tree.get_tree_root()
        max_stat = 0

        some_mallei_node = self.tree & 'BUMA.1121.00036'
        some_pmallei_node = self.tree & 'BUMA.1121.00111'

        for child in root.get_children():
            color_nodes(child, root)
            print('hey!', len(child.get_leaves()))

        print('Root length:', get_genome_sizes(root))
        print('Some mallei length:', get_genome_sizes(some_mallei_node))
        print('Some pseudimallei length:', get_genome_sizes(some_pmallei_node))

        print('Max stat:', max_stat)

        ts = TreeStyle()
        ts.mode = mode
        ts.scale = self.scale
        # Disable the default tip names config
        ts.show_leaf_name = False
        ts.show_branch_support = show_branch_support
        # ts.branch_vertical_margin = 20
        ts.show_scale = show_scale

        self.tree.render(file, w=1000, tree_style=ts)

    def get_all_leafs(self):
        return {node.name for node in self.tree.get_leaves()}

    def recover_internal_states(self, genomes_to_blocks):
        def assign_colorset_feature(v):
            if v.is_leaf():
                bs1 = genomes_to_blocks[(v.name, '1')]
                bs2 = genomes_to_blocks[(v.name, '2')]
                v.add_features(blocks1=(bs1), blocks2=bs2)
                # v.add_features(blocks=genomes_to_blocks[v.name])
            else:
                try:
                    child1, child2 = v.children
                except ValueError:
                    print(v.children)
                    raise ValueError('Tree must me binary')
                c1_blocks1, c1_blocks2 = assign_colorset_feature(child1)
                c2_blocks1, c2_blocks2 = assign_colorset_feature(child2)
                # c1_blocks = assign_colorset_feature(child1)
                # c2_blocks = assign_colorset_feature(child2)

                v.add_features(blocks1=c1_blocks1 | c2_blocks1, blocks2=c1_blocks2 | c2_blocks2)
                # v.add_features(blocks=c1_blocks | c2_blocks)

            return v.blocks1, v.blocks2
            # return v.blocks

        # get colorsets for internal nodes
        root = self.tree.get_tree_root()

        assign_colorset_feature(root)

    def prune(self, ls):
        self.tree.prune(list(ls))
