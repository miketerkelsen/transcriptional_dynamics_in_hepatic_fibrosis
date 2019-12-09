import scanpy.api as sc
import matplotlib as plt
import pandas as pd
import numpy as np
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

sc.settings.set_figure_params(dpi=80, dpi_save=1200, format='pdf', color_map='RdBu')
sc.settings.verbosity = 3
sc.settings.autoshow = False
sc.logging.print_versions()
results_file = '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/aggr_all_norep1/scanpy_anal_15n.h5ad'

aggr_all = sc.read(results_file)
sc.pl.umap(aggr_all, color='Clec4f')
len(aggr_all.obs)
len(aggr_all[aggr_all.obs['CellType'].isin(['Macrophages'])])
aggr_all.write_csvs('./cellphone/', skip_data=False)

print(aggr_all.var)
df = pd.DataFrame(aggr_all.obs)
df.to_csv('aggr_all_af_filt_obs.csv')

np.savetxt('aggr_all_af_filt_matrix.csv', aggr_all.X, delimiter=',')

aggr_all.obs['LibraryID'].replace({'Ctrl - 1': 'Ctrl', 'Ctrl - 2': 'Ctrl', 'Ctrl - 3': 'Ctrl',
                                   'CCl4 2w - 1': 'CCl4 2w', 'CCl4 2w - 2': 'CCl4 2w', 'CCl4 2w - 3': 'CCl4 2w',
                                   'CCl4 4w - 2': 'CCl4 4w', 'CCl4 4w - 3': 'CCl4 4w'}, inplace=True)

# sc.tl.louvain(aggr_all, resolution=0.05)
new_cluster_names = ['Endothelial cells', 'Hepatic stellate cells', 'Macrophages']
aggr_all.rename_categories('louvain', new_cluster_names)

# Count treatment distribution
print(len(aggr_all[aggr_all.obs['LibraryID'].isin(['Ctrl']) & aggr_all.obs['louvain'].isin(['Endothelial cells'])]))

# Treatment overall
color_5 = ['#4F94CD', '#274763', '#A6D6FF']
sc.pl.umap(aggr_all, color='LibraryID', palette=color_5, groups=['Ctrl', 'CCl4 2w', 'CCl4 4w'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/aggr_all_norep1/pub_figures/umap_treatment.pdf', bbox_inches='tight')


color = ['#C02424', '#4F94CD', '#7CCD7C']

sc.pl.tsne(aggr_all, color='louvain', palette=color, groups=['Hepatic stellate cells', 'Endothelial cells', 'Macrophages'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/aggr_all_norep1/pub_figures/tsne_celltype.pdf', bbox_inches='tight')
sc.pl.umap(aggr_all, color='louvain', palette=color, groups=['Hepatic stellate cells', 'Endothelial cells', 'Macrophages'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/aggr_all_norep1/pub_figures/umap_celltype.pdf', bbox_inches='tight')

# Treatment
color1 = ['#2D2D2D', '#D0D0D0', '#D0D0D0']
sc.pl.umap(aggr_all, color='LibraryID', palette=color1, groups=['Ctrl', 'CCl4 2w', 'CCl4 4w'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/aggr_all_norep1/pub_figures/umap_ccl4_2w.pdf', bbox_inches='tight')
color2 = ['#D0D0D0', '#2D2D2D', '#D0D0D0']
sc.pl.umap(aggr_all, color='LibraryID', palette=color2, groups=['Ctrl', 'CCl4 2w', 'CCl4 4w'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/aggr_all_norep1/pub_figures/umap_ccl4_4w.pdf', bbox_inches='tight')
color3 = ['#D0D0D0', '#D0D0D0', '#2D2D2D']
sc.pl.umap(aggr_all, color='LibraryID', palette=color3, groups=['Ctrl', 'CCl4 2w', 'CCl4 4w'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/aggr_all_norep1/pub_figures/umap_ctrl.pdf', bbox_inches='tight')

ptprb_scale = {'red': ((0.0, 0.89, 0.89),
                       (0.005, 0.87, 0.87),
                       (0.25, 0.85, 0.85),
                       (0.5, 0.82, 0.82),
                       (0.75, 0.78, 0.78),
                       (1.0, 0.75, 0.75)),

               'green': ((0.0, 0.89, 0.89),
                         (0.005, 0.77, 0.77),
                         (0.25, 0.64, 0.64),
                         (0.5, 0.52, 0.52),
                         (0.75, 0.27, 0.27),
                         (1.0, 0.14, 0.14)),

               'blue': ((0.0, 0.89, 0.89),
                        (0.005, 0.77, 0.77),
                        (0.25, 0.64, 0.64),
                        (0.5, 0.52, 0.52),
                        (0.75, 0.27, 0.27),
                        (1.0, 0.14, 0.14))
               }

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

adgre1_scale = {'red': ((0.0, 0.89, 0.89),
                        (0.005, 0.83, 0.83),
                        (0.25, 0.76, 0.76),
                        (0.5, 0.69, 0.69),
                        (0.75, 0.55, 0.55),
                        (1.0, 0.49, 0.49)),

                'green': ((0.0, 0.89, 0.89),
                          (0.005, 0.88, 0.88),
                          (0.25, 0.86, 0.86),
                          (0.5, 0.85, 0.85),
                          (0.75, 0.82, 0.82),
                          (1.0, 0.80, 0.80)),

                'blue': ((0.0, 0.89, 0.89),
                         (0.005, 0.83, 0.83),
                         (0.25, 0.76, 0.76),
                         (0.5, 0.69, 0.69),
                         (0.75, 0.55, 0.55),
                         (1.0, 0.49, 0.49))
                }

# Marker genes
plt.pyplot.register_cmap(name='Grey_Red', data=ptprb_scale)
plt.pyplot.rcParams['image.cmap'] = 'Grey_Red'
sc.pl.umap(aggr_all, color=['Ptprb'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/aggr_all_norep1/pub_figures/umap_marker_ptprb.pdf', bbox_inches='tight')

plt.pyplot.register_cmap(name='Grey_Blue', data=lrat_scale)
plt.pyplot.rcParams['image.cmap'] = 'Grey_Blue'
sc.pl.umap(aggr_all, color=['Lrat'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/aggr_all_norep1/pub_figures/umap_marker_lrat.pdf', bbox_inches='tight')

plt.pyplot.register_cmap(name='Grey_Green', data=adgre1_scale)
plt.pyplot.rcParams['image.cmap'] = 'Grey_Green'
sc.pl.umap(aggr_all, color=['Adgre1'])
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/aggr_all_norep1/pub_figures/umap_marker_adgre1.pdf', bbox_inches='tight')

# sc.pl.heatmap(aggr_all, var_names=['Clec4g', 'Aqp1', 'Gpihbp1', 'Ptprb', 'Fcgr2b', 'Kdr', 'Gpr182', 'Egfl7', 'Cldn5', 'Sdpr', 'Lyve1', 'Stab2', 'Nrp1', 'Pde2a', 'Stab1', 'Dcn', 'Cxcl12', 'Colec11', 'Rgs5', 'Rbp1', 'Hand2', 'Lum', 'Mmp2', 'Col1a1', 'Col1a2', 'Col3a1', 'Lrat', 'Reln', 'Sod3', 'Ifitm1', 'Cd74', 'C1qc', 'Wfdc17', 'Cd5l', 'Lyz2', 'Vsig4', 'C1qa', 'Ctss', 'C1qb', 'Tyrobp', 'Clec4f', 'Cd44', 'Aif1', 'Adgre1', 'Spi1'],
#               groupby='louvain', swap_axes=True)
# plt.pyplot.savefig(
#     '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/aggr_all_norep1/pub_figures/heatmap_markers.pdf', bbox_inches='tight')

plt.pyplot.rcParams['image.cmap'] = 'jet'
sc.pl.umap(aggr_all, color='n_counts', palette='jet')
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/aggr_all_norep1/pub_figures/umap_n_counts.pdf', bbox_inches='tight')

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

plt.pyplot.register_cmap(name='Grey_Blue', data=lrat_scale)
plt.pyplot.rcParams['image.cmap'] = 'Grey_Blue'

plt.pyplot.rcParams['image.cmap'] = 'viridis'
sc.pl.umap(aggr_all, color='Tnfrsf11b', palette='viridis')
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/aggr_all_norep1/pub_figures/umap_Tnfrsf11b_virid.pdf', bbox_inches='tight')

sc.pl.umap(aggr_all, color='Tnfrsf11b')
plt.pyplot.savefig(
    '/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/aggr_all_norep1/pub_figures/umap_Tnfrsf11b.pdf', bbox_inches='tight')
