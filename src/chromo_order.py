import os

import pandas as pd
import numpy as np

from Bio import SeqIO
from glob import glob
from collections import defaultdict

import argparse

folder = 'Burkholderia_pseudo_mallei/'
gembase_file = glob(folder + '2-annotate_module/LSTINFO-LSTINFO*')[0]

gembase_df = pd.read_csv(gembase_file, sep='\t')
all_contigs = []
from collections import Counter
cnt = Counter()
contig_lengths = defaultdict(list)

data_2d = []

for _, row in gembase_df.iterrows():
    # print(row['gembase_name'], row['orig_name'])

    contigs = [contig for contig in SeqIO.parse(open(row['orig_name']), 'fasta')][0:2]

    contig_to_old_num = {contig.id : i + 1 for i, contig in enumerate(contigs)}
    contigs.sort(key=len, reverse=True)
    contig_to_new_num = {contig.id : i + 1 for i, contig in enumerate(contigs)}

    for i, contig in enumerate(contigs):
        id = contig.id
        data_2d.append([row['gembase_name'], id, contig_to_old_num[id], contig_to_new_num[id]])


print(data_2d)
data = pd.DataFrame(data=data_2d, columns=['strain', 'contig_id', 'old', 'new'])

data.to_csv('Burkholderia_pseudo_mallei/chromo_order.csv', index=False)
print(data)
print(sum(data.old != data.new))
# print(cnt)
# for chr, lengths in contig_lengths.items():
#     print('Chr:', chr, '   mean length:', np.mean(lengths))