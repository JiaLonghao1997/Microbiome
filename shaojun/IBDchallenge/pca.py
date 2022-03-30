from medic import load

import numpy as np
import pandas as pd
from sklearn import decomposition
from matplotlib import pyplot as plt
from os import makedirs

makedirs('plots', exist_ok=True)


tax_He = load.load_dataset('He', 'Taxonomic')
tax_Schirmer = load.load_dataset('Schirmer', 'Taxonomic')
tax_PMI = load.load_dataset('PMI', 'Taxonomic')
tax = pd.concat((tax_He, tax_Schirmer, tax_PMI))
dataset = np.concatenate([np.repeat('H', len(tax_He)), np.repeat('S', len(tax_Schirmer)), np.repeat('P', len(tax_PMI))])

# Filter low abundance variables
tax = tax.T[tax.max() > 1e-3].T

# PCA on log-transformed data
p = decomposition.PCA(2)
p.fit(np.log10(1e-6 + tax).T)

fig,ax = plt.subplots()
for d in 'HSP':
    label = {'H':'He', 'S':'Schirmer', 'P':'PMI'}[d]
    ax.scatter(p.components_[0][dataset == d], p.components_[1][dataset ==d], s=10, label=label)
    
ax.set_xlabel('PC1 ({:.1%} explained variance)'.format(p.explained_variance_ratio_[0]))
ax.set_ylabel('PC2 ({:.1%} explained variance)'.format(p.explained_variance_ratio_[1]))

ax.legend(loc='best')
fig.savefig('plots/PCA.log.l2.svg')

