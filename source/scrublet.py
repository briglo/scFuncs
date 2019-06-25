module load briglo/miniconda/3
source activate magic
python

import sys ### for later, parses arguments to script (sys.arg)
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

#filtered
#input_dir = '/share/ScratchGeneral/briglo/scRNA/venchi/data/hg19/VENCHI_SampleBCITE/outs/filtered_feature_bc_matrix/'
input_dir = '/share/ScratchGeneral/briglo/scRNA/venchi/outputs/seurat/'
counts_matrix = scipy.io.mmread(input_dir + 'combined.human.mtx').T.tocsc()
genes = np.array(scr.load_genes(input_dir + 'genes.tsv', delimiter='\t', column=0))

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))

#Counts matrix shape: 12865 rows, 32738 columns
#Number of genes in gene list: 32738

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                         n_prin_comps=30,
                                                         get_doublet_neighbor_parents=True)

scrub.plot_histogram();
plt.savefig('/share/ScratchGeneral/briglo/scRNA/venchi/plt.png')
scrub.call_doublets(threshold=0.24)
scrub.plot_histogram();
plt.savefig('/share/ScratchGeneral/briglo/scRNA/venchi/plt.png')

print('Running UMAP...')
# scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
# scrub.plot_embedding('UMAP', order_points=True);
# plt.savefig('/share/ScratchGeneral/briglo/scRNA/venchi/plt.png')

np.savetxt(fname="/share/ScratchGeneral/briglo/scRNA/venchi/outputs/scrublet/Sample_BCITE_doubletScores_seuratCombined_HUMAN.tsv",X=doublet_scores)