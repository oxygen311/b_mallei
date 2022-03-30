import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('losses_tree_3000.csv')

df['losses_chr1'] = [l if l > 10 else 0 for l in df['losses_chr1']]
df['losses_chr2'] = [l if l > 10 else 0 for l in df['losses_chr2']]

df['losses'] = df['losses_chr1'] + df['losses_chr2']

# name,root_length,branch_length,losses_chr1,losses_chr2
x = 'losses'
y = 'root_length'

sns.set_style('whitegrid')
# sns.scatterplot(x='losses', y=y, data=df)
sns.scatterplot(x='losses_chr1', y=y, data=df, label='chr1')
sns.scatterplot(x='losses_chr2', y=y, data=df, label='chr2')

plt.tight_layout()

# plt.xscale('log')
# plt.yscale('log')

plt.savefig(f'scatter_{x}_{y}_2chr.pdf')
plt.show()