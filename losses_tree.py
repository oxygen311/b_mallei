from src.utils.parsers import parse_infercars_to_df, make_labels_dict, get_genomes_contain_blocks_grimm
from src.utils.converters import block_coords_to_infercars
from src.tree.tree_holder import TreeHolder

from itertools import repeat

import logging

logger = logging.getLogger()

res = 5000
mode = 'c'
data_folder = 'data/Burkholderia_pseudo_mallei/7-sibeliaz/sibeliaz_out/fine/%d/' % res
GRIMM_FILENAME = 'genomes_permutations.txt'
BLOCKS_COORD_FILENAME = 'blocks_coords.txt'
INFERCARS_FILENAME = 'blocks_coords.infercars'

tree_file = 'data/Burkholderia_pseudo_mallei/6-tree_module/BUMA_112.nucl.grp.aln.iqtree_tree_right_root.treefile'
labels_file = 'data/Burkholderia_pseudo_mallei/labels.csv'

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # parsing
    tree_holder = TreeHolder(tree_file, logger, labels_dict=make_labels_dict(labels_file))

    block_coords_to_infercars(data_folder + BLOCKS_COORD_FILENAME, data_folder + INFERCARS_FILENAME)
    blocks_df = parse_infercars_to_df(data_folder + INFERCARS_FILENAME)
    blocks_df['length'] = blocks_df['chr_end'] - blocks_df['chr_beg']

    genomes, blocks, block_genome_count, genomes_to_blocks = \
        get_genomes_contain_blocks_grimm(data_folder + GRIMM_FILENAME)

    block_to_length = blocks_df.groupby('block')['length'].mean().to_dict()

    # # processing?
    tree_holder.recover_internal_states(genomes_to_blocks)
    tree_holder.draw('tree_losses_2chr_%d_%s.pdf' % (res, mode), block_to_length, mode=mode, show_branch_support=False)
