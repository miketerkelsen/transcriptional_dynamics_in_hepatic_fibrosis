import scanpy.api as sc
import matplotlib as plt
import numpy as np
from numpy import genfromtxt
import pandas as pd
import scipy.io as si
import random
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

sc.settings.set_figure_params(dpi=80, dpi_save=1200, format='pdf',
                              color_map='Reds')
sc.settings.verbosity = 3
sc.settings.autoshow = False
sc.logging.print_versions()
results_file = '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_anal_DE.h5ad'
result_file = '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_anal_filt.h5ad'


aggr_all = sc.read(results_file)
del aggr_all.uns['louvain_colors']
aggr_all = aggr_all[aggr_all.obs['louvain'].isin(['0', '1', '2', '3'])]
len(aggr_all.obs)
sc.pl.umap(aggr_all, color='louvain')

aggr_all.write(result_file)

# sc.tl.diffmap(aggr_all, n_comps=30)

# sc.pl.diffmap(aggr_all, color='LibraryID')
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/diffmap.pdf', bbox_inches='tight')

df = pd.DataFrame(aggr_all.obs)
df.to_csv('HSC_louvain_new.csv')

sc.pl.umap(aggr_all, color='LibraryID', palette=['#D0D0D0', '#D0D0D0', '#D0D0D0', '#D0D0D0', '#D0D0D0', '#2D2D2D', '#D0D0D0', '#D0D0D0'], groups=[
           'Ctrl - 1', 'Ctrl - 2', 'Ctrl - 3', 'CCl4 2w - 1', 'CCl4 2w - 2', 'CCl4 2w - 3', 'CCl4 4w - 2', 'CCl4 4w - 3'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_library_rep_1.pdf', bbox_inches='tight')
sc.pl.umap(aggr_all, color='LibraryID', palette=['#D0D0D0', '#D0D0D0', '#D0D0D0', '#D0D0D0', '#D0D0D0', '#D0D0D0', '#2D2D2D', '#D0D0D0'], groups=[
           'Ctrl - 1', 'Ctrl - 2', 'Ctrl - 3', 'CCl4 2w - 1', 'CCl4 2w - 2', 'CCl4 2w - 3', 'CCl4 4w - 2', 'CCl4 4w - 3'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_library_rep_2.pdf', bbox_inches='tight')
sc.pl.umap(aggr_all, color='LibraryID', palette=['#D0D0D0', '#D0D0D0', '#D0D0D0', '#D0D0D0', '#D0D0D0', '#D0D0D0', '#D0D0D0', '#2D2D2D'], groups=[
           'Ctrl - 1', 'Ctrl - 2', 'Ctrl - 3', 'CCl4 2w - 1', 'CCl4 2w - 2', 'CCl4 2w - 3', 'CCl4 4w - 2', 'CCl4 4w - 3'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_library_rep_3.pdf', bbox_inches='tight')
sc.pl.umap(aggr_all, color='LibraryID', palette=['#2D2D2D', '#D0D0D0', '#D0D0D0', '#D0D0D0', '#D0D0D0', '#D0D0D0', '#D0D0D0', '#D0D0D0'], groups=[
           'Ctrl - 1', 'Ctrl - 2', 'Ctrl - 3', 'CCl4 2w - 1', 'CCl4 2w - 2', 'CCl4 2w - 3', 'CCl4 4w - 2', 'CCl4 4w - 3'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_library_rep_4.pdf', bbox_inches='tight')
sc.pl.umap(aggr_all, color='LibraryID', palette=['#D0D0D0', '#2D2D2D', '#D0D0D0', '#D0D0D0', '#D0D0D0', '#D0D0D0', '#D0D0D0', '#D0D0D0'], groups=[
           'Ctrl - 1', 'Ctrl - 2', 'Ctrl - 3', 'CCl4 2w - 1', 'CCl4 2w - 2', 'CCl4 2w - 3', 'CCl4 4w - 2', 'CCl4 4w - 3'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_library_rep_5.pdf', bbox_inches='tight')
sc.pl.umap(aggr_all, color='LibraryID', palette=['#D0D0D0', '#D0D0D0', '#2D2D2D', '#D0D0D0', '#D0D0D0', '#D0D0D0', '#D0D0D0', '#D0D0D0'], groups=[
           'Ctrl - 1', 'Ctrl - 2', 'Ctrl - 3', 'CCl4 2w - 1', 'CCl4 2w - 2', 'CCl4 2w - 3', 'CCl4 4w - 2', 'CCl4 4w - 3'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_library_rep_6.pdf', bbox_inches='tight')
sc.pl.umap(aggr_all, color='LibraryID', palette=['#D0D0D0', '#D0D0D0', '#D0D0D0', '#2D2D2D', '#D0D0D0', '#D0D0D0', '#D0D0D0', '#D0D0D0'], groups=[
           'Ctrl - 1', 'Ctrl - 2', 'Ctrl - 3', 'CCl4 2w - 1', 'CCl4 2w - 2', 'CCl4 2w - 3', 'CCl4 4w - 2', 'CCl4 4w - 3'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_library_rep_7.pdf', bbox_inches='tight')
sc.pl.umap(aggr_all, color='LibraryID', palette=['#D0D0D0', '#D0D0D0', '#D0D0D0', '#D0D0D0', '#2D2D2D', '#D0D0D0', '#D0D0D0', '#D0D0D0'], groups=[
           'Ctrl - 1', 'Ctrl - 2', 'Ctrl - 3', 'CCl4 2w - 1', 'CCl4 2w - 2', 'CCl4 2w - 3', 'CCl4 4w - 2', 'CCl4 4w - 3'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_library_rep_8.pdf', bbox_inches='tight')


# color = ['#A37C27', '#A7414A', '#6A8A82']
aggr_all.obs['LibraryID'].replace({'Ctrl - 1': 'Ctrl', 'Ctrl - 2': 'Ctrl', 'Ctrl - 3': 'Ctrl',
                                   'CCl4 2w - 1': 'CCl4 2w', 'CCl4 2w - 2': 'CCl4 2w', 'CCl4 2w - 3': 'CCl4 2w',
                                   'CCl4 4w - 2': 'CCl4 4w', 'CCl4 4w - 3': 'CCl4 4w'}, inplace=True)

new_treatment_color = ['#4F94CD', '#274763', '#BCE1FF']

sc.pl.umap(aggr_all, color='LibraryID', palette=new_treatment_color,
           groups=['Ctrl', 'CCl4 2w', 'CCl4 4w'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_library_new_color.pdf', bbox_inches='tight')

sc.pl.tsne(aggr_all, color='LibraryID', palette=new_treatment_color,
           groups=['Ctrl', 'CCl4 2w', 'CCl4 4w'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/tsne_library.pdf', bbox_inches='tight')

# np.savetxt("HSC_louvain.csv", aggr_all.obs, delimiter=',', fmt='%s')

# sc.tl.louvain(aggr_all, resolution=0.6)

colors = ['#B22222', '#B47746', '#50A85E', '#33689E', '#50A4A8', '#674CC7', '#C64CC7', '#C74C83']

# colors = ['#DD4949', '#DD7D49', '#DDB249', '#40DD4C', '#49D2DD', '#498CDD', '#7A49DD', '#DD49B0', '#DD4978']

# sc.pl.tsne(aggr_all, color='louvain', palette=colors)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/tsne_louvain.pdf', bbox_inches='tight')
sc.pl.umap(aggr_all, color='louvain', palette=colors)
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_louvain_new_color.pdf', bbox_inches='tight')

# sc.tl.rank_genes_groups(aggr_all, groupby='louvain', method='wilcoxon')
# sc.pl.rank_genes_groups(aggr_all, n_genes=50, sharey=False)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/cluster_marker_list.pdf', bbox_inches='tight')

# color = ['maroon', 'firebrick', 'red', 'darkorange', 'orange', 'goldenrod',
#          'darkolivegreen', 'olivedrab', 'olive']

# sc.pl.tsne(aggr_all, color='LibraryID', palette=color, groups=['Ctrl - 1', 'Ctrl - 2', 'Ctrl - 3', 'CCl4 2w - 1', 'CCl4 2w - 2',
#                                                                'CCl4 2w - 3', 'CCl4 4w - 1', 'CCl4 4w - 2', 'CCl4 4w - 3'])
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_clean_clean_clean/scanpy_figures/tsne_library_clear.png')
# sc.pl.umap(aggr_all, color='LibraryID', palette=color, groups=['Ctrl - 1', 'Ctrl - 2', 'Ctrl - 3', 'CCl4 2w - 1', 'CCl4 2w - 2',
#                                                                'CCl4 2w - 3', 'CCl4 4w - 1', 'CCl4 4w - 2', 'CCl4 4w - 3'])
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_clean_clean_clean/scanpy_figures/umap_library_clear.png')

plt.pyplot.rcParams['image.cmap'] = 'jet'
sc.pl.umap(aggr_all, color='n_counts', palette='jet')
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_n_counts.pdf', bbox_inches='tight')

# cluster_markers_0 = genfromtxt('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/markers_cluster_0.csv', delimiter=',', dtype=str)
# sc.pl.umap(aggr_all, color=cluster_markers_0)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_cluster_markers_plot_0.png', bbox_inches='tight')

# cluster_markers_1 = genfromtxt('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/markers_cluster_1.csv', delimiter=',', dtype=str)
# sc.pl.umap(aggr_all, color=cluster_markers_1)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_cluster_markers_plot_1.png', bbox_inches='tight')

# cluster_markers_2 = genfromtxt('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/markers_cluster_2.csv', delimiter=',', dtype=str)
# sc.pl.umap(aggr_all, color=cluster_markers_2)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_cluster_markers_plot_2.png', bbox_inches='tight')

# cluster_markers_3 = genfromtxt('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/markers_cluster_3.csv', delimiter=',', dtype=str)
# sc.pl.umap(aggr_all, color=cluster_markers_3)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_cluster_markers_plot_3.png', bbox_inches='tight')

# cluster_markers_4 = genfromtxt('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/markers_cluster_4.csv', delimiter=',', dtype=str)
# sc.pl.umap(aggr_all, color=cluster_markers_4)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_cluster_markers_plot_4.png', bbox_inches='tight')

# Export random HSC subset
row_index = pd.read_csv(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/genes_HSC_norep1.tsv', header=None, sep='\t')[0]
col_index = pd.read_csv(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/barcodes_HSC_norep1.tsv', header=None, sep='\t')[0]
ccl4 = si.mmread(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/matrix_HSC_norep1.mtx')
B = ccl4.todense()
df = pd.DataFrame(B, index=row_index.values, columns=col_index.values)

df = df.T
len(df)
rand = random.sample(range(12392), 6196)
df = df[col_index.values[rand]]
np.savetxt('/Users/mikekrogh/Desktop/subset_HSC.csv', rand, delimiter=',')

np.savetxt('/Users/mikekrogh/Desktop/barcodes_HSC.tsv',
           df.columns.values.tolist(), delimiter='\t', fmt='%s')
np.savetxt('/Users/mikekrogh/Desktop/genes_HSC.tsv',
           df.index.values.tolist(), delimiter='\t', fmt='%s')

df.to_csv('/Users/mikekrogh/Desktop/genes_HSC.csv', sep=',', index=False)
si.mmwrite('/Users/mikekrogh/Desktop/matrix_HSC', df)

lrat_scale = {'red': ((0.0, 0.89, 0.89),
                      (0.005, 0.80, 0.80),
                      (0.25, 0.70, 0.70),
                      (0.5, 0.60, 0.60),
                      (0.75, 0.41, 0.41),
                      (1.0, 0.31, 0.31)),

              'green': ((0.0, 0.89, 0.89),
                        (0.005, 0.84, 0.84),
                        (0.25, 0.79, 0.79),
                        (0.5, 0.74, 0.74),
                        (0.75, 0.63, 0.63),
                        (1.0, 0.58, 0.58)),

              'blue': ((0.0, 0.89, 0.89),
                       (0.005, 0.88, 0.88),
                       (0.25, 0.86, 0.86),
                       (0.5, 0.85, 0.85),
                       (0.75, 0.82, 0.82),
                       (1.0, 0.80, 0.80))
              }

plt.pyplot.register_cmap(name='Grey_Blue', data=lrat_scale)
plt.pyplot.rcParams['image.cmap'] = 'Grey_Blue'
sc.pl.umap(aggr_all, color='Col1a1')
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_marker_col1a1.pdf', bbox_inches='tight')

# UMAP stage assignments
stage_id = pd.read_csv(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/HSC_stages_500.csv', sep=',')
aggr_all.obs['Stage_1'] = stage_id['Stage_1'].values
aggr_all.obs['Stage_2'] = stage_id['Stage_2'].values
aggr_all.obs['Stage_3'] = stage_id['Stage_3'].values
aggr_all.obs['Stage_4'] = stage_id['Stage_4'].values

sc.pl.umap(aggr_all, color='Stage_1', palette=['#D0D0D0', '#E72525'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_stage_1_new.pdf', bbox_inches='tight')

sc.pl.umap(aggr_all, color='Stage_2', palette=['#D0D0D0', '#E72525'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_stage_2.pdf', bbox_inches='tight')

sc.pl.umap(aggr_all, color='Stage_3', palette=['#D0D0D0', '#E72525'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_stage_3.pdf', bbox_inches='tight')

sc.pl.umap(aggr_all, color='Stage_4', palette=['#D0D0D0', '#E72525'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_stage_4.pdf', bbox_inches='tight')

sc.tl.rank_genes_groups(aggr_all, groupby='Stage_4', method='wilcoxon')
sc.pl.rank_genes_groups(aggr_all, n_genes=50, sharey=False)
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/cluster_marker_list.pdf', bbox_inches='tight')
aggr_all.uns['rank_genes_groups']


# Output genes for manuscript
lrat_scale = {'red': ((0.0, 0.89, 0.89),
                      (0.005, 0.80, 0.80),
                      (0.10, 0.70, 0.70),
                      (0.30, 0.60, 0.60),
                      (0.50, 0.41, 0.41),
                      (0.70, 0.31, 0.31),
                      (1.00, 0.13, 0.13)),

              'green': ((0.0, 0.89, 0.89),
                        (0.005, 0.84, 0.84),
                        (0.10, 0.79, 0.79),
                        (0.30, 0.74, 0.74),
                        (0.50, 0.63, 0.63),
                        (0.70, 0.58, 0.58),
                        (1.00, 0.32, 0.32)),

              'blue': ((0.0, 0.89, 0.89),
                       (0.005, 0.88, 0.88),
                       (0.10, 0.86, 0.86),
                       (0.30, 0.85, 0.85),
                       (0.50, 0.82, 0.82),
                       (0.70, 0.80, 0.80),
                       (1.00, 0.55, 0.55))
              }

lrat_scale2 = {'red': ((0.0, 0.89, 0.89),
                      (0.005, 0.80, 0.80),
                      (0.10, 0.70, 0.70),
                      (0.20, 0.60, 0.60),
                      (0.30, 0.41, 0.41),
                      (0.40, 0.31, 0.31),
                      (1.00, 0.09, 0.09)),

              'green': ((0.0, 0.89, 0.89),
                        (0.005, 0.84, 0.84),
                        (0.10, 0.79, 0.79),
                        (0.20, 0.74, 0.74),
                        (0.30, 0.63, 0.63),
                        (0.40, 0.58, 0.58),
                        (1.00, 0.26, 0.26)),

              'blue': ((0.0, 0.89, 0.89),
                       (0.005, 0.88, 0.88),
                       (0.10, 0.86, 0.86),
                       (0.20, 0.85, 0.85),
                       (0.30, 0.82, 0.82),
                       (0.40, 0.80, 0.80),
                       (1.00, 0.41, 0.41))
              }

plt.pyplot.register_cmap(name='Grey_Blue', data=lrat_scale2)
plt.pyplot.rcParams['image.cmap'] = 'Grey_Blue'

sc.pl.umap(aggr_all, color=['Fcna', 'Fn1', 'Mmp2', 'Dpep1'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_genes.pdf', bbox_inches='tight')

df = pd.DataFrame(aggr_all.raw.var)
df.to_csv('HSC_genes.csv')
from scipy import sparse

len(aggr_all.raw.var)
np.amax(aggr_all.raw.X.todense()[:, 6795])
aggr_all.raw.X.shape[1]
np.savetxt('HSC_plot_val.csv', aggr_all.raw.X.todense(), delimiter=',')
aggr_all.raw.X.todense()

sc.pl.umap(aggr_all, color=['Plvap', 'Col1a1'], vmin=0, vmax=3.85)
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_plvap2.pdf', bbox_inches='tight')

sc.pl.tsne(aggr_all, color=['Col1a1'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/tsne_plvap.pdf', bbox_inches='tight')


eigengenes = pd.read_csv(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/HSC_eigengenes.csv', sep=',')
eigengenes.head()
aggr_all.obs['Module 1 Eigengene'] = eigengenes['Module 1 Eigengene'].values
aggr_all.obs['Module 2 Eigengene'] = eigengenes['Module 2 Eigengene'].values
aggr_all.obs['Module 3 Eigengene'] = eigengenes['Module 3 Eigengene'].values

sc.pl.umap(aggr_all, color=['Module 1 Eigengene', 'Module 2 Eigengene', 'Module 3 Eigengene'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/HSC_norep1/scanpy_figures/umap_modules_eigengenes.pdf', bbox_inches='tight')
