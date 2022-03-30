from parebrick.utils.data.parsers import parse_infercars_to_df, export_df_to_infercars

import pandas as pd

file = '7-sibeliaz/infercars_old/1000/blocks_coords.infercars'
file_pseudomallei = file.replace('.infercars_old', '_pseudomallei.infercars_old')
file_mallei = file.replace('.infercars_old', '_mallei.infercars_old')

df = parse_infercars_to_df(file)

labels = pd.read_csv('labels.csv')

pseudo_mallei = set(strain for strain, label
                    in zip(labels['strain'], labels['label'])
                    if 'pseudomallei' in label)

pseudo_mallei_df = df[df.apply(lambda x: x['species'] in pseudo_mallei, axis=1)]
mallei_df = df[df.apply(lambda x: x['species'] not in pseudo_mallei, axis=1)]

export_df_to_infercars(pseudo_mallei_df, file_pseudomallei)
export_df_to_infercars(mallei_df, file_mallei)