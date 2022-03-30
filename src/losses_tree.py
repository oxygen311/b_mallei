from src.utils.parsers import parse_infercars_to_df, make_labels_dict, get_genomes_contain_blocks_grimm
from src.utils.converters import block_coords_to_infercars
from src.tree.tree_holder import TreeHolder

import logging
import json
import pandas as pd

logger = logging.getLogger()

res = 3000
mode = 'r'
data_folder = 'data/Burkholderia_pseudo_mallei/7-sibeliaz/sibeliaz_out_all/fine/%d/' % res
GRIMM_FILENAME = 'genomes_permutations.txt'
BLOCKS_COORD_FILENAME = 'blocks_coords.txt'
INFERCARS_FILENAME = 'blocks_coords.infercars_old'

ALLOWED_BLOCKS_FILE = 'data/Burkholderia_pseudo_mallei/7-sibeliaz/infercars_old/%d/allowed_blocks_20.json' % res

tree_file = 'data/Burkholderia_pseudo_mallei/6-tree_module/BUMA_112.nucl.grp.aln.iqtree_tree_right_root.treefile'
labels_file = 'data/Burkholderia_pseudo_mallei/labels.csv'

def prune_mallei(tree_holder):
    labels = pd.read_csv('data/Burkholderia_pseudo_mallei/labels.csv')
    mallei = [strain for strain, label in zip(labels['strain'], labels['label']) if 'pseudomallei' not in label
                                                                                    or 'VB29710' in label]
    tree_holder.tree.prune(mallei)

if __name__ == '__main__':
    # parsing
    tree_holder = TreeHolder(tree_file, logger, labels_dict=make_labels_dict(labels_file), scale=300000)

    with open(ALLOWED_BLOCKS_FILE) as json_file:
        allowed_blocks = set(json.load(json_file))

    block_coords_to_infercars(data_folder + BLOCKS_COORD_FILENAME, data_folder + INFERCARS_FILENAME)
    blocks_df = parse_infercars_to_df(data_folder + INFERCARS_FILENAME)

    # filter for allowed blocks
    blocks_df = blocks_df[blocks_df.apply(lambda x: x['block'] in allowed_blocks, axis=1)]

    blocks_df['length'] = blocks_df['chr_end'] - blocks_df['chr_beg']

    genomes, blocks, block_genome_count, genomes_to_blocks = \
        get_genomes_contain_blocks_grimm(data_folder + GRIMM_FILENAME, allowed_blocks)

    block_to_length = blocks_df.groupby('block')['length'].mean().to_dict()

    # # processing?
    tree_holder.recover_internal_states(genomes_to_blocks)
    prune_mallei(tree_holder)
    tree_holder.draw('mallei_tree_losses_2chr_%d_%s.pdf' % (res, mode), block_to_length, mode=mode, show_branch_support=False)
