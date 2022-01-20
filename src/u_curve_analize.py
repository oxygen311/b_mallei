from src.utils.parsers import parse_infercars_to_df, make_labels_dict, get_genomes_contain_blocks_grimm

import json

res = 5000
species = 'pseudomallei'
data_folder = 'data/Burkholderia_pseudo_mallei/7-sibeliaz/infercars/%d/' % res
INFERCARS_FILENAME = 'blocks_coords_%s.infercars' % species
THRESHHOLD = 20

labels_file = 'data/Burkholderia_pseudo_mallei/labels.csv'

df_blocks = parse_infercars_to_df(data_folder + INFERCARS_FILENAME)

block_count = df_blocks.groupby('block')['species'].nunique().sort_values()

allowed_blocks = [b for b, copies in block_count.iteritems() if copies >= THRESHHOLD]

with open(data_folder + 'allowed_blocks_%d.json' % THRESHHOLD, 'w') as outfile:
    json.dump(allowed_blocks, outfile)