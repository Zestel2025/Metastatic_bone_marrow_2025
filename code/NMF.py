#!/usr/bin/env python
# coding: utf-8

# In[33]:


def limit_threads(workers: int = 4):
    import os, torch
    os.environ["OMP_NUM_THREADS"] = str(workers)
    os.environ["OPENBLAS_NUM_THREADS"] = str(workers)
    os.environ["MKL_NUM_THREADS"] = str(workers)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(workers)
    os.environ["NUMEXPR_NUM_THREADS"] = str(workers)
    try:
        torch.set_num_threads(workers)
        torch.set_num_interop_threads(max(1, workers // 2))
    except Exception:
        pass
    print(f"⚙️  Using at most {workers} threads")
limit_threads(12)


# In[3]:


import scanpy as sc
import omicverse as ov
ov.plot_set()


# In[45]:


adata = sc.read("../../p_data/origin/Epi.h5ad")


# In[46]:


adata = adata.raw.to_adata()


# In[47]:


adata=ov.pp.preprocess(adata,mode='shiftlog|pearson',n_HVGs=3000,)
adata


# In[52]:


import numpy as np
# adata.X = np.nan_to_num(adata.X, nan=0.0)
ov.pp.scale(adata)
# np.isnan(adata.X).sum()
ov.pp.pca(adata)


# In[54]:


import matplotlib.pyplot as plt
from matplotlib import patheffects
# fig, ax = plt.subplots(figsize=(4,4))
ov.pl.embedding(
    adata,
    basis="X_umap",
    color=['re_annotation_2',"harmony_cluster_2","cancer_site","cancer_type"],
    frameon='small',
    title="Celltypes",
    #legend_loc='on data',
    legend_fontsize=14,
    legend_fontoutline=2,
    #size=10,
#     ax=ax,
    #legend_loc=True, 
    add_outline=False, 
    #add_outline=True,
    outline_color='black',
    outline_width=1,
    show=False,
)


# In[61]:


import numpy as np
## Initialize the cnmf object that will be used to run analyses
cnmf_obj = ov.single.cNMF(adata,components=np.arange(5,11), n_iter=20, seed=14, num_highvar_genes=3000,
                          output_dir='results', name='dg_cNMF')


# In[63]:


## Specify that the jobs are being distributed over a single worker (total_workers=1) and then launch that worker
cnmf_obj.factorize(worker_i=0, total_workers=2)


# In[ ]:


import pickle
with open("./results/cnmf_obj_1.pkl", "wb") as f:
    pickle.dump(cnmf_obj, f)


# In[ ]:


# # 读取
# import pickle
# with open("./results/cnmf_obj_1.pkl", "rb") as f:
#     cnmf_obj = pickle.load(f)


# In[ ]:


cnmf_obj.combine(skip_missing_files=True)


# In[ ]:


cnmf_obj.k_selection_plot(close_fig=False)
plt.savefig(f"./results/stability.pdf",bbox_inches='tight')


# In[ ]:


selected_K = 7
density_threshold = 2.00


# In[ ]:


cnmf_obj.consensus(k=selected_K, 
                   density_threshold=density_threshold, 
                   show_clustering=True, 
                   close_clustergram_fig=False)
plt.savefig(f"./results/Heatmap_pre_1.pdf",bbox_inches='tight')


# In[ ]:


density_threshold = 0.10


# In[ ]:


cnmf_obj.consensus(k=selected_K, 
                   density_threshold=density_threshold, 
                   show_clustering=True, 
                   close_clustergram_fig=False)
plt.savefig(f"./results/Heatmap_pre_2.pdf",bbox_inches='tight')


# In[ ]:


import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import patheffects

from matplotlib import gridspec
import matplotlib.pyplot as plt

width_ratios = [0.2, 4, 0.5, 10, 1]
height_ratios = [0.2, 4]
fig = plt.figure(figsize=(sum(width_ratios), sum(height_ratios)))
gs = gridspec.GridSpec(len(height_ratios), len(width_ratios), fig,
                        0.01, 0.01, 0.98, 0.98,
                       height_ratios=height_ratios,
                       width_ratios=width_ratios,
                       wspace=0, hspace=0)
            
D = cnmf_obj.topic_dist[cnmf_obj.spectra_order, :][:, cnmf_obj.spectra_order]
dist_ax = fig.add_subplot(gs[1,1], xscale='linear', yscale='linear',
                                      xticks=[], yticks=[],xlabel='', ylabel='',
                                      frameon=True)
dist_im = dist_ax.imshow(D, interpolation='none', cmap='viridis',
                         aspect='auto', rasterized=True)

left_ax = fig.add_subplot(gs[1,0], xscale='linear', yscale='linear', xticks=[], yticks=[],
                xlabel='', ylabel='', frameon=True)
left_ax.imshow(cnmf_obj.kmeans_cluster_labels.values[cnmf_obj.spectra_order].reshape(-1, 1),
                            interpolation='none', cmap='Spectral', aspect='auto',
                            rasterized=True)

top_ax = fig.add_subplot(gs[0,1], xscale='linear', yscale='linear', xticks=[], yticks=[],
                xlabel='', ylabel='', frameon=True)
top_ax.imshow(cnmf_obj.kmeans_cluster_labels.values[cnmf_obj.spectra_order].reshape(1, -1),
                  interpolation='none', cmap='Spectral', aspect='auto',
                    rasterized=True)

cbar_gs = gridspec.GridSpecFromSubplotSpec(3, 3, subplot_spec=gs[1, 2],
                                   wspace=0, hspace=0)
cbar_ax = fig.add_subplot(cbar_gs[1,2], xscale='linear', yscale='linear',
    xlabel='', ylabel='', frameon=True, title='Euclidean\nDistance')
cbar_ax.set_title('Euclidean\nDistance',fontsize=12)
vmin = D.min().min()
vmax = D.max().max()
fig.colorbar(dist_im, cax=cbar_ax,
        ticks=np.linspace(vmin, vmax, 3),
        )
cbar_ax.set_yticklabels(cbar_ax.get_yticklabels(),fontsize=12)
plt.savefig(f"./results/Heatmap_1.pdf",bbox_inches='tight')
