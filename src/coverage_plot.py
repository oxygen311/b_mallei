from parebrick.utils.data.parsers import parse_infercars_to_df, export_df_to_infercars

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from collections import Counter


file = 'data/Burkholderia_pseudo_mallei/7-sibeliaz/sibeliaz_out_all/fine/3000/coverage_report.txt'
labels = pd.read_csv('data/Burkholderia_pseudo_mallei/labels.csv')

df = pd.read_csv(file, sep='\t', skiprows=67, header=None)

pseudo_mallei = set(strain for strain, label
                    in zip(labels['strain'], labels['label'])
                    if 'pseudomallei' in label)

cov_data_2d = []

print(pseudo_mallei)

for i, row in df.iterrows():
    genome, chr = row[0].rsplit('.', 1)
    species = 'pseudomallei' if genome in pseudo_mallei else 'mallei'
    cov_data_2d.append([genome, species, chr, row[1]])

cov_df = pd.DataFrame(cov_data_2d, columns=['strain', 'species', 'chr', 'coverage'])
print(Counter(cov_df.species.values))

sns.set(style='whitegrid', font_scale=1.2)
sns.boxplot(data=cov_df, x='species', y='coverage', hue='chr')

plt.tight_layout()
plt.savefig('coverage.pdf')
plt.show()