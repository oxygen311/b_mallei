import os.path

import pandas as pd
import scipy.stats

import numpy as np

from utils.parsers import export_df_to_infercars, parse_infercars_to_df
# from parebrick.utils.data.unique_gene_filters import filter_dataframe_unique

prefix = 'BUMA'

folder_genes = 'data/Burkholderia_pseudo_mallei/2-annotate_module/LSTINFO/'
lst_info_columns = ['start', 'end', 'orientation', 'type', 'gene_id', 'gene', 'desc']
# species = 'pseudomallei'
species = 'mallei'

labels_file = 'data/Burkholderia_pseudo_mallei/labels.csv'
og_file = f'data/Burkholderia_pseudo_mallei/{species}.proteinortho.tsv'

res = 3000
blocks_file = f'data/Burkholderia_pseudo_mallei/7-sibeliaz/infercars/{res}/blocks_coords_{species}.infercars'
genes_file = f'data/Burkholderia_pseudo_mallei/genes_{species}.infercars'


def diff_len_chr1_chr2():
    df = parse_infercars_to_df(genes_file)

    vs1 = df[df.chr == '1']['chr_end'] - df[df.chr == '1']['chr_beg']
    vs2 = df[df.chr == '2']['chr_end'] - df[df.chr == '2']['chr_beg']

    print(len(vs1), np.mean(vs1), np.percentile(vs1, 50))
    print(len(vs2), np.mean(vs2), np.percentile(vs2, 50))

    print(scipy.stats.mannwhitneyu(vs1, vs2))
    print(scipy.stats.ttest_ind(vs1, vs2, equal_var=False))


def is_there_common_genes_on_both_chr():
    df = parse_infercars_to_df(genes_file)
    all_sp = len(df['species'].unique())

    cnt = 0
    for block, df_block in df.groupby('block'):
        genomes_count = len(df_block['species'].unique())
        chrs_count = len(df_block['chr'].unique())

        if genomes_count == all_sp and chrs_count > 1:
            cnt += 1
            print(genomes_count, chrs_count)

    print(cnt)


def sd_length_for_every_block():
    df = parse_infercars_to_df(blocks_file)
    df['length'] = df['chr_end'] - df['chr_beg']

    for block, df_block in df.groupby('block'):
        sd = np.std(df_block.length.values)
        if sd > 50:
            print(block, sd)



def test_mallei_core():
    df = parse_infercars_to_df(blocks_file)

    block_to_len = {b: np.mean(df_b.chr_end - df_b.chr_beg) for b, df_b in df.groupby('block')}
    all_sp = len(df['species'].unique())

    for chr, df_chr in df.groupby('chr'):
        len_acc = 0
        for block, df_block in df_chr.groupby('block'):
            if len(df_block['species'].unique()) == all_sp:
                len_acc += block_to_len[block]

        print(chr, len_acc)


def filter_dataframe_unique(df):
    allowed_blocks = set()
    all_sp = len(df['species'].unique())

    for block, df_block in df.groupby('block'):
        if len(df_block['species'].unique()) == all_sp:
            allowed_blocks.add(block)

    return df.loc[df['block'].isin(allowed_blocks)].copy()


def test_mallei_core_2():
    df = parse_infercars_to_df(blocks_file)
    df = filter_dataframe_unique(df)

    ln_2d = []
    for (genome, chr), df_tmp in df.groupby(['species', 'chr']):
        ln = sum(df_tmp['chr_end'] - df_tmp['chr_beg'])
        ln_2d.append([genome, chr, ln])

    ln_df = pd.DataFrame(ln_2d, columns=['strain', 'chr', 'length'])

    for chr, df_chr in ln_df.groupby('chr'):
        print(chr, np.mean(df_chr.length))


if __name__ == "__main__":
    test_mallei_core_2()
