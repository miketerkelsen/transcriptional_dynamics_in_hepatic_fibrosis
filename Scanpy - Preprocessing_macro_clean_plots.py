import scanpy.api as sc
import matplotlib as plt
import numpy as np
from numpy import genfromtxt
import pandas as pd
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

sc.settings.set_figure_params(dpi=80, dpi_save=1200, format='pdf',
                              color_map='Reds')
sc.settings.verbosity = 2
sc.settings.autoshow = False
sc.logging.print_versions()
results_file = '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/scanpy_anal_DE.h5ad'

aggr_all_im = sc.read(results_file)
sc.pl.umap(aggr_all, color='louvain')
# print(aggr_all_im.obs)
# np.savetxt("Macrophage_louvain.csv", aggr_all.obs, delimiter=',', fmt='%s')
del aggr_all_im.uns['louvain_colors']
aggr_all = aggr_all_im[aggr_all_im.obs['louvain'].isin(['0', '1', '2', '3'])]
len(aggr_all.obs)

# sc.tl.louvain(aggr_all, resolution=0.4)
df = pd.DataFrame(aggr_all.obs)
df.to_csv('macro_louvain_new.csv')

# color = ['#A37C27', '#A7414A', '#6A8A82']
aggr_all.obs['LibraryID'].replace({'Ctrl - 1': 'Ctrl', 'Ctrl - 2': 'Ctrl', 'Ctrl - 3': 'Ctrl',
                                   'CCl4 2w - 1': 'CCl4 2w', 'CCl4 2w - 2': 'CCl4 2w', 'CCl4 2w - 3': 'CCl4 2w',
                                   'CCl4 4w - 2': 'CCl4 4w', 'CCl4 4w - 3': 'CCl4 4w'}, inplace=True)

new_treatment_color = ['#7CCD7C', '#3C633C', '#C4FDC4']
sc.pl.umap(aggr_all, color='LibraryID', palette=new_treatment_color,
           groups=['Ctrl', 'CCl4 2w', 'CCl4 4w'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/scanpy_figures/umap_library_new_color.pdf', bbox_inches='tight')

colors = ['#B22222', '#B47746', '#50A85E', '#33689E', '#50A4A8', '#674CC7', '#C64CC7', '#C74C83']

# colors = ['#DD4949', '#DD7D49', '#DDB249', '#40DD4C', '#49D2DD', '#498CDD', '#7A49DD', '#DD49B0', '#DD4978']

# sc.pl.tsne(aggr_all, color='louvain', palette=colors)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/scanpy_figures/tsne_louvain_res06.pdf', bbox_inches='tight')
sc.pl.umap(aggr_all_im, color='louvain', palette=colors)
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/scanpy_figures/umap_louvain_all_pops.pdf', bbox_inches='tight')

# sc.tl.louvain(aggr_all, resolution=0.6)
# aggr_all.write(result_file)

# sc.tl.louvain(aggr_all, resolution=0.2)

# aggr_subset = aggr_all[aggr_all.obs['louvain'].isin(['0', '1'])]

# sc.tl.tsne(aggr_subset, n_pcs=18)
# sc.pp.neighbors(aggr_subset, n_neighbors=10, n_pcs=18)
# sc.tl.umap(aggr_subset)


# colors = ['tomato', 'coral', 'darkorange', 'goldenrod', 'olive', 'limegreen', 'mediumseagreen', 'lightseagreen',
#           'deepskyblue', 'dodgerblue', 'mediumpurple', 'orchid', 'hotpink', 'fuchsia', 'crimson']

# sc.pl.tsne(aggr_subset, color='louvain', palette=colors)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/scanpy_figures/tsne_louvain_clean.png', bbox_inches='tight')
# sc.pl.umap(aggr_all, color='louvain', palette=colors)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/scanpy_figures/umap_louvain_marker.png', bbox_inches='tight')


# color = ['maroon', 'firebrick', 'red', 'darkorange', 'orange',
#          'darkolivegreen', 'olivedrab', 'olive']

# sc.pl.tsne(aggr_subset, color='LibraryID', palette=color, groups=['Ctrl - 1', 'Ctrl - 2', 'Ctrl - 3', 'CCl4 2w - 1', 'CCl4 2w - 2',
#                                                                   'CCl4 2w - 3', 'CCl4 4w - 2', 'CCl4 4w - 3'])
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/scanpy_figures/tsne_library_clean.png', bbox_inches='tight')
# sc.pl.umap(aggr_subset, color='LibraryID', palette=color, groups=['Ctrl - 1', 'Ctrl - 2', 'Ctrl - 3', 'CCl4 2w - 1', 'CCl4 2w - 2',
#                                                                   'CCl4 2w - 3', 'CCl4 4w - 2', 'CCl4 4w - 3'])
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/scanpy_figures/umap_library_clean.png', bbox_inches='tight')

sc.tl.rank_genes_groups(aggr_all, groupby='louvain', method='wilcoxon')
sc.pl.rank_genes_groups(aggr_all, n_genes=50, sharey=False, fontsize=6)
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/scanpy_figures/cluster_marker_list.pdf', bbox_inches='tight')

# sc.pl.rank_genes_groups_stacked_violin(aggr_all, n_genes=50)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/scanpy_figures/cluster_marker_list.pdf', bbox_inches='tight')

plt.pyplot.rcParams['image.cmap'] = 'jet'
sc.pl.umap(aggr_all, color='n_counts', palette='jet')
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/scanpy_figures/umap_n_counts.pdf', bbox_inches='tight')

# sc.pl.umap(aggr_all, color=['Hmgb1', 'Stmn1', 'Birc5', 'Tuba1b', 'H2afz', 'Ptma', 'Tubb5', 'Ube2c', 'Rad21', 'H2afv'])
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/scanpy_figures/umap_cluster_markers_1.pdf', bbox_inches='tight')

# cluster_markers_0 = genfromtxt('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/markers_cluster_0.csv', delimiter=',', dtype=str)
# sc.pl.umap(aggr_all, color=cluster_markers_0)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/scanpy_figures/umap_cluster_markers_plot_0.png', bbox_inches='tight')

# cluster_markers_1 = genfromtxt('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/markers_cluster_1.csv', delimiter=',', dtype=str)
# sc.pl.umap(aggr_all, color=cluster_markers_1)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/scanpy_figures/umap_cluster_markers_plot_1.png', bbox_inches='tight')

# cluster_markers_2 = genfromtxt('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/markers_cluster_2.csv', delimiter=',', dtype=str)
# sc.pl.umap(aggr_all, color=cluster_markers_2)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/scanpy_figures/umap_cluster_markers_plot_2.png', bbox_inches='tight')

# cluster_markers_3 = genfromtxt('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/markers_cluster_3.csv', delimiter=',', dtype=str)
# sc.pl.umap(aggr_all, color=cluster_markers_3)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/scanpy_figures/umap_cluster_markers_plot_3.png', bbox_inches='tight')

# cluster_markers_4 = genfromtxt('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/markers_cluster_4.csv', delimiter=',', dtype=str)
# sc.pl.umap(aggr_all, color=cluster_markers_4)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/scanpy_figures/umap_cluster_markers_plot_4.png', bbox_inches='tight')

lrat_scale = {'red': ((0.0, 0.89, 0.89),
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

adgre1_scale = {'red': ((0.0, 0.89, 0.89),
                        (0.005, 0.83, 0.83),
                        (0.10, 0.76, 0.76),
                        (0.20, 0.69, 0.69),
                        (0.30, 0.55, 0.55),
                        (0.40, 0.49, 0.49),
                        (1.00, 0.09, 0.09)),

                'green': ((0.0, 0.89, 0.89),
                          (0.005, 0.88, 0.88),
                          (0.10, 0.86, 0.86),
                          (0.20, 0.85, 0.85),
                          (0.30, 0.82, 0.82),
                          (0.40, 0.80, 0.80),
                          (1.00, 0.42, 0.42)),

                'blue': ((0.0, 0.89, 0.89),
                         (0.005, 0.83, 0.83),
                         (0.10, 0.76, 0.76),
                         (0.20, 0.69, 0.69),
                         (0.30, 0.55, 0.55),
                         (0.40, 0.49, 0.49),
                         (1.00, 0.09, 0.09))
                }

plt.pyplot.register_cmap(name='Grey_Red', data=adgre1_scale)
plt.pyplot.rcParams['image.cmap'] = 'Grey_Red'

sc.pl.umap(aggr_all, color=['Trem2'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Macrophage_norep1/scanpy_figures/umap_ms_genes_2.pdf', bbox_inches='tight')

plt.pyplot.rcParams['image.cmap'] = 'RdBu'
sc.pl.heatmap(aggr_all_im, ['Clec4f', 'Vsig4', 'Ccr2', 'Cx3cr1', 'Hmgb2', 'Stmn1', 'Birc5', 'Tuba1b', 'Ube2c', 'Hmgn2', 'Cenpa', 'Ccnb2', 'Cks2', 'Tmpo', 'Cks1b', 'Mki67', 'Racgap1', 'Hist1h2ap', 'Ccna2', 'Top2a'], groupby = 'louvain', swap_axes = True)
