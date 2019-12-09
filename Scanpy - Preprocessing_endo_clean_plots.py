import scanpy.api as sc
import matplotlib as plt
import numpy as np
from numpy import genfromtxt
import pandas as pd
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

sc.settings.set_figure_params(dpi=80, dpi_save=1200, format='pdf',
                              color_map='Reds')
sc.settings.verbosity = 3
sc.settings.autoshow = False
sc.logging.print_versions()
results_file = '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/scanpy_anal_DE.h5ad'

aggr_all = sc.read(results_file)
# np.savetxt("Endothelial_louvain.csv", aggr_all.obs, delimiter=',', fmt='%s')
len(aggr_all.obs)

sc.pl.umap(aggr_all, color='louvain')
df = pd.DataFrame(aggr_all.obs)
df.to_csv('endo_louvain_new.csv')

color = ['#A37C27', '#A7414A', '#6A8A82']
aggr_all.obs['LibraryID'].replace({'Ctrl - 1': 'Ctrl', 'Ctrl - 2': 'Ctrl', 'Ctrl - 3': 'Ctrl',
                                   'CCl4 2w - 1': 'CCl4 2w', 'CCl4 2w - 2': 'CCl4 2w', 'CCl4 2w - 3': 'CCl4 2w',
                                   'CCl4 4w - 2': 'CCl4 4w', 'CCl4 4w - 3': 'CCl4 4w'}, inplace=True)

new_treatment_color = ['#C02424', '#631313', '#FF8989']
sc.pl.umap(aggr_all, color='LibraryID', palette=new_treatment_color,
           groups=['Ctrl', 'CCl4 2w', 'CCl4 4w'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/scanpy_figures/umap_library_new_color.pdf', bbox_inches='tight')

sc.pl.tsne(aggr_all, color='LibraryID', palette=new_treatment_color,
           groups=['Ctrl', 'CCl4 2w', 'CCl4 4w'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/scanpy_figures/tsne_library_new_color.pdf', bbox_inches='tight')

# sc.tl.louvain(aggr_all, resolution=0.4)

colors = ['#B22222', '#B47746', '#50A85E', '#33689E', '#50A4A8', '#674CC7', '#C64CC7', '#C74C83']

# colors = ['#DD4949', '#DD7D49', '#DDB249', '#40DD4C', '#49D2DD', '#498CDD', '#7A49DD', '#DD49B0', '#DD4978']

# sc.pl.tsne(aggr_all, color='louvain', palette=colors)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/scanpy_figures/tsne_louvain.pdf', bbox_inches='tight')
sc.pl.umap(aggr_all, color='louvain', palette=colors)
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/scanpy_figures/umap_louvain_new_color.pdf', bbox_inches='tight')

# sc.tl.rank_genes_groups(aggr_all, groupby='louvain', method='wilcoxon')
# sc.pl.rank_genes_groups(aggr_all, n_genes=50, sharey=False)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/scanpy_figures/cluster_marker_list.pdf', bbox_inches='tight')

# color = ['maroon', 'firebrick', 'red', 'darkorange', 'orange',
#          'darkolivegreen', 'olivedrab', 'olive']

# sc.pl.tsne(aggr_all, color='LibraryID', palette=color, groups=['Ctrl - 1', 'Ctrl - 2', 'Ctrl - 3', 'CCl4 2w - 1', 'CCl4 2w - 2',
#                                                                'CCl4 2w - 3', 'CCl4 4w - 2', 'CCl4 4w - 3'])
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/scanpy_figures/tsne_library_clear.png', bbox_inches='tight')
# sc.pl.umap(aggr_all, color='LibraryID', palette=color, groups=['Ctrl - 1', 'Ctrl - 2', 'Ctrl - 3', 'CCl4 2w - 1', 'CCl4 2w - 2',
#                                                                'CCl4 2w - 3', 'CCl4 4w - 2', 'CCl4 4w - 3'])
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/scanpy_figures/umap_library_clear.png', bbox_inches='tight')

plt.pyplot.rcParams['image.cmap'] = 'jet'
sc.pl.umap(aggr_all, color='n_counts', palette='jet')
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/scanpy_figures/umap_n_counts.pdf', bbox_inches='tight')

# cluster_markers_0 = genfromtxt('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/markers_cluster_0.csv', delimiter=',', dtype=str)
# sc.pl.umap(aggr_all, color=cluster_markers_0)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/scanpy_figures/umap_cluster_markers_plot_0.png', bbox_inches='tight')

# cluster_markers_1 = genfromtxt('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/markers_cluster_1.csv', delimiter=',', dtype=str)
# sc.pl.umap(aggr_all, color=cluster_markers_1)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/scanpy_figures/umap_cluster_markers_plot_1.png', bbox_inches='tight')

# cluster_markers_2 = genfromtxt('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/markers_cluster_2.csv', delimiter=',', dtype=str)
# sc.pl.umap(aggr_all, color=cluster_markers_2)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/scanpy_figures/umap_cluster_markers_plot_2.png', bbox_inches='tight')

# cluster_markers_3 = genfromtxt('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/markers_cluster_3.csv', delimiter=',', dtype=str)
# sc.pl.umap(aggr_all, color=cluster_markers_3)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/scanpy_figures/umap_cluster_markers_plot_3.png', bbox_inches='tight')

# cluster_markers_4 = genfromtxt('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/markers_cluster_4.csv', delimiter=',', dtype=str)
# sc.pl.umap(aggr_all, color=cluster_markers_4)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/scanpy_figures/umap_cluster_markers_plot_4.png', bbox_inches='tight')

# cluster_markers_5 = genfromtxt('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/markers_cluster_5.csv', delimiter=',', dtype=str)
# sc.pl.umap(aggr_all, color=cluster_markers_5)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/scanpy_figures/umap_cluster_markers_plot_5.png', bbox_inches='tight')

# cluster_markers_6 = genfromtxt('/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/markers_cluster_6.csv', delimiter=',', dtype=str)
# sc.pl.umap(aggr_all, color=cluster_markers_6)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/scanpy_figures/umap_cluster_markers_plot_6.png', bbox_inches='tight')

ptprb_scale = {'red': ((0.0, 0.89, 0.89),
                       (0.005, 0.87, 0.87),
                       (0.10, 0.85, 0.85),
                       (0.30, 0.82, 0.82),
                       (0.50, 0.78, 0.78),
                       (0.70, 0.75, 0.75),
                       (1.00, 0.55, 0.55)),

               'green': ((0.0, 0.89, 0.89),
                         (0.005, 0.77, 0.77),
                         (0.10, 0.64, 0.64),
                         (0.30, 0.52, 0.52),
                         (0.50, 0.27, 0.27),
                         (0.70, 0.14, 0.14),
                         (1.00, 0.13, 0.13)),

               'blue': ((0.0, 0.89, 0.89),
                        (0.005, 0.77, 0.77),
                        (0.10, 0.64, 0.64),
                        (0.30, 0.52, 0.52),
                        (0.50, 0.27, 0.27),
                        (0.70, 0.14, 0.14),
                        (1.00, 0.13, 0.13))
               }

ptprb_scale2 = {'red': ((0.0, 0.89, 0.89),
                       (0.005, 0.87, 0.87),
                       (0.10, 0.85, 0.85),
                       (0.20, 0.82, 0.82),
                       (0.30, 0.78, 0.78),
                       (0.40, 0.75, 0.75),
                       (1.00, 0.41, 0.41)),

               'green': ((0.0, 0.89, 0.89),
                         (0.005, 0.77, 0.77),
                         (0.10, 0.64, 0.64),
                         (0.20, 0.52, 0.52),
                         (0.30, 0.27, 0.27),
                         (0.40, 0.14, 0.14),
                         (1.00, 0.09, 0.09)),

               'blue': ((0.0, 0.89, 0.89),
                        (0.005, 0.77, 0.77),
                        (0.10, 0.64, 0.64),
                        (0.20, 0.52, 0.52),
                        (0.30, 0.27, 0.27),
                        (0.40, 0.14, 0.14),
                        (1.00, 0.09, 0.09))
               }

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

plt.pyplot.register_cmap(name='Grey_Red', data=ptprb_scale2)
plt.pyplot.rcParams['image.cmap'] = 'Grey_Red'
sc.pl.umap(aggr_all, color=['Plvap', 'Col1a1'], vmin=0, vmax=3.85)
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/scanpy_figures/umap_plvap_blue.pdf', bbox_inches='tight')

sc.pl.umap(aggr_all, color=['Stab2', 'Rspo3', 'Adgrg6', 'Efnb1', 'Cd34', 'Ednrb', 'Csrp2', 'Ltbp4'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/scanpy_figures/umap_ms_genes.pdf', bbox_inches='tight')

sc.pl.tsne(aggr_all, color=['Plvap', 'Lrat'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/Endothelial_norep1/scanpy_figures/tsne_plvap.pdf', bbox_inches='tight')
