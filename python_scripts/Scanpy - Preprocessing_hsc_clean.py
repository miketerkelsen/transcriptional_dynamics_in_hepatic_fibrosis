import numpy as np
import pandas as pd
import scanpy.api as sc
import matplotlib as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

sc.settings.set_figure_params(dpi=80, dpi_save=600, format='pdf',
                              color_map='Reds')
sc.settings.verbosity = 3
sc.settings.autoshow = False
sc.logging.print_versions()
results_file = '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_anal2.h5ad'

path = '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/'
aggr_all = sc.read(path + 'matrix_HSC_norep1.mtx', cache=True).T
aggr_all.var_names = pd.read_csv(path + 'genes_HSC_norep1.tsv', header=None, sep='\t')[0]
aggr_all.obs_names = pd.read_csv(path + 'barcodes_HSC_norep1.tsv', header=None, sep='\t')[0]

aggr_all.var_names_make_unique()

LibraryID = pd.read_csv('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/HSC_norep1.csv', sep=',')
aggr_all.obs['LibraryID'] = LibraryID['LibraryID'].values

sc.pl.highest_expr_genes(aggr_all, n_top=30)
plt.pyplot.savefig('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/highest_expr_genes.pdf', bbox_inches='tight')

sc.pp.filter_genes(aggr_all, min_cells=50)
sc.pp.filter_cells(aggr_all, min_genes=200)

mito_genes = [name for name in aggr_all.var_names if name.startswith('mt.')]
aggr_all.obs['percent_mito'] = np.sum(
    aggr_all[:, mito_genes].X, axis=1).A1 / np.sum(aggr_all.X, axis=1).A1
aggr_all.obs['n_counts'] = aggr_all.X.sum(axis=1).A1

sc.pl.violin(aggr_all, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4,
             multi_panel=True)
plt.pyplot.savefig('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/violin.pdf', bbox_inches='tight')
sc.pl.scatter(aggr_all, x='n_counts', y='n_genes')
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/scatter_counts_genes.pdf', bbox_inches='tight')
sc.pl.scatter(aggr_all, x='n_counts', y='percent_mito')
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/scatter_counts_mito.pdf', bbox_inches='tight')

aggr_all.raw = sc.pp.log1p(aggr_all, copy=True)

aggr_all = aggr_all[aggr_all.obs['n_genes'] < 2500, :]
aggr_all = aggr_all[aggr_all.obs['percent_mito'] < 0.05, :]

sc.pp.normalize_per_cell(aggr_all)

sc.pp.log1p(aggr_all)

# dt = pd.DataFrame(data=aggr_all.var)
# index = dt.index.values
# listOfRowIndexLabels = list(index)
# np.savetxt("HSC_background_genes.csv", listOfRowIndexLabels, delimiter=',', fmt='%s')

sc.pp.highly_variable_genes(aggr_all, min_mean=0.0125, max_mean=3, min_disp=0.5, max_disp=20)
sc.pl.highly_variable_genes(aggr_all)
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/filter_genes_dispersion.pdf', bbox_inches='tight')
aggr_all = aggr_all[:, aggr_all.var['highly_variable']]

sc.pp.regress_out(aggr_all, ['n_counts', 'percent_mito'], n_jobs=1)

sc.pp.scale(aggr_all)

# PCA
sc.tl.pca(aggr_all, svd_solver='arpack')
sc.pl.pca_scatter(aggr_all)
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/pca_scatter.pdf', bbox_inches='tight')

sc.pl.pca_variance_ratio(aggr_all, log=True, n_pcs=40)
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/pca_variance_ratio.pdf', bbox_inches='tight')

# tSNE
sc.tl.tsne(aggr_all, n_pcs=25)
sc.pl.tsne(aggr_all)
plt.pyplot.savefig('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/tsne.pdf', bbox_inches='tight')
sc.pl.tsne(aggr_all, color=['Lrat', 'Ptprb', 'Adgre1'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/tsne_lrat_ptprb_adgre1.pdf', bbox_inches='tight')

# Compute neightbours
sc.pp.neighbors(aggr_all, n_neighbors=10, n_pcs=25)
sc.tl.umap(aggr_all)
sc.pl.umap(aggr_all)
plt.pyplot.savefig('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap.pdf', bbox_inches='tight')
sc.pl.umap(aggr_all, color=['Lrat', 'Ptprb', 'Adgre1'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_lrat_ptprb_adgre1.pdf', bbox_inches='tight')

aggr_all.write(results_file)

# Clustering the graph
sc.tl.louvain(aggr_all, resolution=0.1)

colors = ['tomato', 'coral', 'darkorange', 'goldenrod', 'olive', 'limegreen', 'mediumseagreen', 'lightseagreen',
          'deepskyblue', 'dodgerblue', 'mediumpurple', 'orchid', 'hotpink', 'fuchsia', 'crimson']

sc.pl.tsne(aggr_all, color='louvain', palette=colors)
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/tsne_louvain.pdf', bbox_inches='tight')
sc.pl.umap(aggr_all, color='louvain', palette=colors)
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_louvain.pdf', bbox_inches='tight')

color = ['maroon', 'firebrick', 'red', 'darkorange', 'orange',
         'darkolivegreen', 'olivedrab', 'olive']

sc.pl.tsne(aggr_all, color='LibraryID', palette=color, groups=['Ctrl - 1', 'Ctrl - 2', 'Ctrl - 3', 'CCl4 2w - 1', 'CCl4 2w - 2',
                                                               'CCl4 2w - 3', 'CCl4 4w - 2', 'CCl4 4w - 3'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/tsne_library.pdf', bbox_inches='tight')
sc.pl.umap(aggr_all, color='LibraryID', palette=color, groups=['Ctrl - 1', 'Ctrl - 2', 'Ctrl - 3', 'CCl4 2w - 1', 'CCl4 2w - 2',
                                                               'CCl4 2w - 3', 'CCl4 4w - 2', 'CCl4 4w - 3'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_library.pdf', bbox_inches='tight')

aggr_all.write(results_file)



# Analysis of states in HSC pseudotime
results_file = '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_anal_DE.h5ad'

aggr_all = sc.read(results_file)
ini_sub = pd.read_csv('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/cell_state_analysis/initial_state_bars.csv', sep=',')
mid_2w = pd.read_csv('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/cell_state_analysis/middle_2w_state.csv', sep=',')
mid_sub = pd.read_csv('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/cell_state_analysis/middle_state.csv', sep=',')
late_sub = pd.read_csv('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/cell_state_analysis/late_state.csv', sep=',')

HSC_ini = pd.DataFrame(aggr_all.obs.index)
HSC_ini.to_csv('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/cell_state_analysis/HSC_barcodes.csv')
HSC_comp_list = pd.read_csv('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/cell_state_analysis/HSC_comp_list.csv', header=None, sep=',')[0]

aggr_all.obs['Ctrl state'] = pd.Categorical(HSC_ini['index'].isin(ini_sub['barcodes']))
sc.pl.umap(aggr_all, color='Ctrl state', palette=['grey', 'red'], alpha=0.6)
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/cell_state_analysis/umap_ctrl_state.pdf', bbox_inches='tight')

aggr_all.obs['Mid 2w state'] = pd.Categorical(HSC_ini['index'].isin(mid_2w['barcodes']))
sc.pl.umap(aggr_all, color='Mid 2w state', palette=['grey', 'red'], alpha=0.6)
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/cell_state_analysis/umap_mid_2w_state.pdf', bbox_inches='tight')

aggr_all.obs['Mid state'] = pd.Categorical(HSC_ini['index'].isin(mid_sub['barcodes']))
sc.pl.umap(aggr_all, color='Mid state', palette=['grey', 'red'], alpha=0.6)
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/cell_state_analysis/umap_mid_state.pdf', bbox_inches='tight')

aggr_all.obs['Late state'] = pd.Categorical(HSC_ini['index'].isin(late_sub['barcodes']))
sc.pl.umap(aggr_all, color='Late state', palette=['grey', 'red'], alpha=0.6)
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/cell_state_analysis/umap_late_state.pdf', bbox_inches='tight')

sc.tl.rank_genes_groups(aggr_all, groupby='Ctrl state', method='wilcoxon')
sc.pl.rank_genes_groups(aggr_all, n_genes=50, sharey=False, fontsize=6)
plt.pyplot.savefig('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/cell_state_analysis/Ctrl_state_genes.pdf', bbox_inches='tight')
np.savetxt('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/cell_state_analysis/Ctrl_state_markers.csv', [aggr_all.uns['rank_genes_groups']], delimiter=',', fmt='%s')

sc.tl.rank_genes_groups(aggr_all, groupby='Mid 2w state', method='wilcoxon')
sc.pl.rank_genes_groups(aggr_all, n_genes=50, sharey=False, fontsize=6)
plt.pyplot.savefig('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/cell_state_analysis/Middle_2w_genes.pdf', bbox_inches='tight')
np.savetxt('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/cell_state_analysis/Middle_2w_markers.csv', [aggr_all.uns['rank_genes_groups']], delimiter=',', fmt='%s')

sc.tl.rank_genes_groups(aggr_all, groupby='Mid state', method='wilcoxon')
sc.pl.rank_genes_groups(aggr_all, n_genes=50, sharey=False, fontsize=6)
plt.pyplot.savefig('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/cell_state_analysis/Middle_state_genes.pdf', bbox_inches='tight')
np.savetxt('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/cell_state_analysis/Middle_state_markers.csv', [aggr_all.uns['rank_genes_groups']], delimiter=',', fmt='%s')

sc.tl.rank_genes_groups(aggr_all, groupby='Late state', method='wilcoxon')
sc.pl.rank_genes_groups(aggr_all, n_genes=50, sharey=False, fontsize=6)
plt.pyplot.savefig('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/cell_state_analysis/Late_state_genes.pdf', bbox_inches='tight')
np.savetxt('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/cell_state_analysis/Late_state_markers.csv', [aggr_all.uns['rank_genes_groups']], delimiter=',', fmt='%s')
