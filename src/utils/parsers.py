import os
import re

import numpy as np
import pandas as pd

from itertools import takewhile
from io import StringIO
from bg.grimm import GRIMMReader
from collections import defaultdict, Counter

PATTERN = re.compile("([A-Za-z0-9_\(\)\/\s\.-]+)\.([A-Za-z0-9_]+):(\d+)-(\d+) ([+|-]).*")
COLUMNS = ["block", "species", "chr", "chr_beg", "chr_end", "orientation"]
BLOCKS_SEPARATOR = '-' * 80


def parse_infercars_to_df(file_name):
    def find_indices(lst, condition):
        return [i for i, elem in enumerate(lst) if condition(elem)]

    with open(file_name) as f:
        lines = f.readlines()

    last_line = len(lines) - 1
    while lines[last_line] == '\n': last_line -= 1

    n_at_end = len(lines) - 1 - last_line
    for _ in range(1 - n_at_end): lines.append('\n')

    bs = np.split(lines, find_indices(lines, lambda x: x[0] == ">"))
    temp = []

    for i, b in enumerate(bs):
        if len(b) == 0: continue
        b_i = int(b[0][1:])

        for oc in b[1:-1]:
            m = PATTERN.match(oc)
            temp.append([b_i, m.group(1), m.group(2), int(m.group(3)), int(m.group(4)), m.group(5)])

    return pd.DataFrame(temp, columns=COLUMNS)


def export_df_to_infercars(df, file_name):
    with open(file_name, 'w') as f:
        for block, block_df in df.groupby('block'):
            print(f'>{block}', file=f)
            for i, row in block_df.iterrows():
                print(f'{row["species"]}.{row["chr"]}:{row["chr_beg"]}-{row["chr_end"]} {row["orientation"]}', file=f)
            print(file=f)


def genome_lengths_from_block_coords(in_file):
    with open(in_file) as f:
        head_lines = list(takewhile(lambda line: (line != BLOCKS_SEPARATOR + os.linesep) and
                                                 (line != BLOCKS_SEPARATOR + '\n'), f))

    # names of chromosomes
    df_head = pd.read_csv(StringIO(''.join(head_lines)), sep='\t')
    return {row['Description']: row['Size'] for index, row in df_head.iterrows()}

def genome_genome_lengths_from_chromosomes_lengths(chr_lengths):
    lengths = defaultdict(int)

    for chr, len in chr_lengths.items():
        strain, _ = chr.rsplit('.', 1)
        lengths[strain] += len

    return lengths


def get_genomes_contain_blocks_grimm(grimm_file, allowed_blocks):
    genomes, blocks = set(), set()

    with open(grimm_file) as f:
        ls = f.readlines()
    block_genome_count, genomes_to_blocks = defaultdict(Counter), {}
    # block_genome_count, genomes_to_blocks = defaultdict(Counter), defaultdict(Counter)

    for i in range(0, len(ls), 2):
        name, chromo = GRIMMReader.parse_genome_declaration_string(ls[i]).name.rsplit('.', 1)
        # print(name, chromo)
        data = GRIMMReader.parse_data_string(ls[i + 1])[1]
        genomes.add(name)
        for _, block in data:
            if int(block) in allowed_blocks:
                blocks.add(int(block))
                block_genome_count[int(block)][name] += 1
        genomes_to_blocks[(name, chromo)] = Counter([block for _, block in data if int(block) in allowed_blocks])
        # genomes_to_blocks[name] += Counter([block for _, block in data])

    return list(sorted(genomes)), list(sorted(blocks)), block_genome_count, genomes_to_blocks


def make_labels_dict(file, row_from='strain', row_to='label'):
    try:
        df = pd.read_csv(file)
        return {row[row_from]: row[row_to] for i, row in df.iterrows()}
    except FileNotFoundError:
        return None