import os.path

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from itertools import chain
from glob import glob
from bisect import bisect
from collections import defaultdict, Counter
from utils.parsers import export_df_to_infercars, parse_infercars_to_df

prefix = 'BUMA'

folder_genes = 'data/Burkholderia_pseudo_mallei/2-annotate_module/LSTINFO/'
lst_info_columns = ['start', 'end', 'orientation', 'type', 'gene_id', 'gene', 'desc']
# species = 'pseudomallei'
species = 'mallei'

labels_file = 'data/Burkholderia_pseudo_mallei/labels.csv'
og_file = f'data/Burkholderia_pseudo_mallei/{species}.proteinortho.tsv'

res = 3000
blocks_file = f'data/Burkholderia_pseudo_mallei/7-sibeliaz/infercars/{res}/blocks_coords_{species}.infercars'
annotated_blocks_file = f'data/Burkholderia_pseudo_mallei/7-sibeliaz/annotated/{res}/blocks_coords_{species}.csv'

chromo_order_file = 'data/Burkholderia_pseudo_mallei/chromo_order.csv'
infercars_file = f'data/Burkholderia_pseudo_mallei/genes_{species}.infercars'

genome_lengths_file = 'data/Burkholderia_pseudo_mallei/8-parebrick/3000/preprocessed_data/genomes_lengths.csv'


def get_pseudo_labels(file):
    labels = pd.read_csv(file)

    return set(strain for strain, label
               in zip(labels['strain'], labels['label'])
               # if 'pseudomallei' not in label)
               if 'pseudomallei' in label)


def parse_genes_df(allowed_strains):
    genes_dfs = []
    for lst_file in glob(folder_genes + '*.lst'):
        strain = os.path.basename(lst_file).replace('.lst', '')
        if strain not in allowed_strains: continue

        genes_df = pd.read_csv(lst_file, sep='\t', header=None, names=lst_info_columns)
        genes_dfs.append(genes_df)

    return pd.concat(genes_dfs)


def annotate_genes_og(genes_df):
    ogs = pd.read_csv(og_file, sep='\t')

    gene_to_group = {}
    for i, row in ogs.iterrows():
        for gene_v in row.values:
            if type(gene_v) == str and gene_v[:4] == prefix:
                genes = gene_v.split(',')
                for gene in genes:
                    gene_to_group[gene] = i + 1

    genes_df['og'] = [gene_to_group.get(id, '') for id in genes_df['gene_id']]


def order_chromos(genes_df, order_file):
    genes_df['strain'] = [id[:15] for id in genes_df.gene_id]

    order_df = pd.read_csv(order_file)
    chromo_order_dict = {(strain, old): new for strain, old, new in zip(order_df.strain, order_df.old, order_df.new)}

    # FILTER CHROMO 3
    genes_df = genes_df[
        (genes_df.gene_id.str.slice(19, 20) == '1') | (genes_df.gene_id.str.slice(19, 20) == '2')].copy()

    genes_df['chromo'] = [chromo_order_dict[(strain, int(id[19]))]
                          for id, strain in zip(genes_df.gene_id, genes_df.strain)]

    return genes_df


def export_infercars(genes_df, infercars_file):
    genes_df = genes_df.drop(columns=['type', 'gene_id', 'gene', 'desc'])

    genes_df.rename(columns={'start': 'chr_beg', 'end': 'chr_end', 'strain': 'species', 'og': 'block', 'chromo': 'chr'},
                    inplace=True)
    genes_df['orientation'] = ['+' if o == 'D' else '-' for o in genes_df['orientation']]
    genes_df = genes_df[genes_df.block != ''].copy()

    export_df_to_infercars(genes_df, infercars_file)


def match_blocks_genes(blocks_df, genes_df):
    def get_ogs(genes_df, inverted=False):
        ogs = [0 if og == '' else (og if or_ == 'D' else -og) for og, or_ in zip(genes_df.og, genes_df.orientation)]
        return [-og for og in ogs[::-1]] if inverted else ogs

    block_to_left_border_gene = {}
    block_to_right_border_gene = {}
    block_to_genes = {}
    gene_to_blocks = defaultdict(list)

    for strain, genes_strain_df in genes_df.groupby('strain'):
        for chr, genes_strain_chr_df in genes_strain_df.groupby('chromo'):
            blocks_strain_df = blocks_df[(blocks_df['species'] == strain) & (blocks_df['chr'] == str(chr))] \
                .sort_values(by=['chr_beg'])
            genes_strain_chr_df = genes_strain_chr_df.sort_values('start')

            for b_id, b_start, b_end, b_or in zip(blocks_strain_df.block, blocks_strain_df.chr_beg,
                                                  blocks_strain_df.chr_end, blocks_strain_df.orientation):

                genes_start_index = bisect(genes_strain_chr_df.start.values, b_start)
                genes_start_index = genes_start_index - 1 if genes_start_index > 0 else 0
                genes_end_index = bisect(genes_strain_chr_df.start.values, b_end)

                current_genes = genes_strain_chr_df[genes_start_index:genes_end_index]

                genes_start_index, genes_end_index = 0, len(current_genes)

                start_is_on_border = current_genes.start.values[0] < b_start
                end_is_on_border = current_genes.end.values[-1] > b_end

                if start_is_on_border: genes_start_index += 1
                if end_is_on_border: genes_end_index -= 1

                inverted = b_or == '-'
                left_ogs = get_ogs(current_genes[0:genes_start_index], inverted)
                center_ogs = get_ogs(current_genes[genes_start_index:genes_end_index], inverted)
                right_ogs = get_ogs(current_genes[genes_end_index:], inverted)

                [gene_to_blocks[abs(og)].append(b_id) for og in center_ogs]

                block_to_genes[(strain, b_id)] = center_ogs
                block_to_left_border_gene[(strain, b_id)], block_to_right_border_gene[(strain, b_id)] = \
                    (right_ogs, left_ogs) if inverted else (left_ogs, right_ogs)

    blocks_df['genes_tail'] = [', '.join(map(str, block_to_left_border_gene[(b_strain, b_id)]))
                               for b_id, b_strain in zip(blocks_df.block, blocks_df.species)]
    blocks_df['genes_inside'] = [', '.join(map(str, block_to_genes[(b_strain, b_id)]))
                                 for b_id, b_strain in zip(blocks_df.block, blocks_df.species)]
    blocks_df['genes_head'] = [', '.join(map(str, block_to_right_border_gene[(b_strain, b_id)]))
                               for b_id, b_strain in zip(blocks_df.block, blocks_df.species)]

    genes_df.rename(columns={'chr_beg': 'start', 'chr_end': 'end', 'species': 'strain'}, inplace=True)
    # blocks_df.to_csv(annotated_blocks_file, index=False)

    return gene_to_blocks


def visualize_gene_to_blocks(og_to_blocks):
    for og, blocks in og_to_blocks.items():
        cnt = Counter(blocks)
        if len(cnt) > 1:
            print(og, Counter(blocks))

    blocks_count = [len(set(bs)) for g, bs in og_to_blocks.items() if g != 0]

    sns.set_style('whitegrid')
    sns.histplot(blocks_count, binwidth=1, log_scale=(False, True))

    plt.tight_layout()
    plt.savefig('genes_in_block_occ.pdf')
    plt.show()


def construct_genome_length(file):
    genome_lengths_df = pd.read_csv(file)
    return {tuple(g.rsplit('.', 1)): l for g, l in zip(genome_lengths_df.Genome, genome_lengths_df.Length)}


def calculate_coverages(blocks_df, genes_df, genome_lengths, filter_singletons=False):
    def construct_locus_labels(df, c_start, c_end, label):
        return [(p, 's', label) for p in df[c_start]], [(p, 'e', label) for p in df[c_end]]

    type_coverages_2d = []
    type_coverages_flat = []
    columns = []

    if filter_singletons:
        genes_df = genes_df[genes_df.og != ''].copy()

    chr_keep = 1
    for strain, genes_strain_df in genes_df.groupby('strain'):
        for chr, genes_strain_chr_df in genes_strain_df.groupby('chromo'):
            if chr != chr_keep: continue

            blocks_strain_df = blocks_df[(blocks_df['species'] == strain) & (blocks_df['chr'] == str(chr))] \
                .sort_values(by=['chr_beg'])
            genes_strain_chr_df = genes_strain_chr_df.sort_values('start')

            ls1, le1 = construct_locus_labels(blocks_strain_df, 'chr_beg', 'chr_end', 'block')
            ls2, le2 = construct_locus_labels(genes_strain_chr_df, 'start', 'end', 'gene')

            covered = defaultdict(int)
            events = list(sorted(chain(ls1, le1, ls2, le2)))

            prev = 0
            state = {'block': 0, 'gene': 0}

            for cur_pos, event, type in events:
                cov_str = ', '.join(t for t, c in state.items() if c > 0)
                covered[cov_str] += cur_pos - prev

                if event == 's':
                    state[type] += 1
                elif event == 'e':
                    state[type] -= 1
                prev = cur_pos

            chr_len = genome_lengths[(strain, str(chr))]

            type_coverages_flat.append([covered[t] / chr_len * 100 for t in sorted(covered.keys()) if t != ''])
            columns = [t for t in sorted(covered.keys()) if t != '']
            for t, c in covered.items():
                if t == '': continue
                type_coverages_2d.append([strain, chr, t, c / chr_len * 100])

    cov_df = pd.DataFrame(type_coverages_2d, columns=['strain', 'chr', 'type', 'covered'])
    cov_df_flat = pd.DataFrame(type_coverages_flat, columns=columns).sort_values('block, gene')

    sns.set_style('whitegrid')

    plt.stackplot(range(len(cov_df_flat)),
                  cov_df_flat['block, gene'],
                  cov_df_flat['block'],
                  cov_df_flat['gene'],
                  colors=['#7dc0a6', '#ed936b', '#919fc7'],
                  labels=['block, gene','block','gene'])

    # sns.boxplot(data=cov_df, x='type', y='covered', hue='chr', order=['block', 'gene', 'block, gene'])

    plt.legend(loc='best', ncol=3)
    plt.ylim(ymax=107)
    plt.ylabel('covered')
    plt.xlabel('starins')
    # plt.xlabel('type')

    plt.tight_layout()
    # plt.savefig(f'coverages_perc_boxplot_{species}_no_s.pdf')
    plt.savefig(f'coverages_perc_area_{species}_chr_{chr_keep}_no_s.pdf')
    plt.show()


def main():
    genome_lengths = construct_genome_length(genome_lengths_file)
    allowed_strains = get_pseudo_labels(labels_file)

    genes_df = parse_genes_df(allowed_strains)
    annotate_genes_og(genes_df)

    genes_df = order_chromos(genes_df, chromo_order_file)
    # export_infercars(genes_df, infercars_file=infercars_file)

    blocks_df = parse_infercars_to_df(blocks_file)
    gene_to_blocks = match_blocks_genes(blocks_df, genes_df)
    visualize_gene_to_blocks(gene_to_blocks)

    calculate_coverages(blocks_df, genes_df, genome_lengths, filter_singletons=True)


if __name__ == "__main__":
    main()
