#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import pertpy as pt
import omicverse as ov
import matplotlib.pyplot as plt
import scanpy as sc
import scvelo as scv
ov.plot_set()

# import data
adata = ov.read(f"../../p_data/Total/Total_250919.h5ad")
adata.X = adata.layers["counts"].copy()
sc.pp.normalize_total(adata, target_sum=50 * 1e4)
sc.pp.log1p(adata)
adata.layers["lognorm"] = adata.X.copy()

# In[175]:


ov.utils.embedding(adata,
                   basis='X_umap',frameon='small',wspace=1,ncols=2,hspace=1,
                   color=['orig.ident','re_annotation_2',"cancer_site"],
                   show=True)


# In[176]:


print(adata.obs['re_annotation_2'].value_counts())


# In[177]:


condition = 'cancer_type'
ctrl_group='KIRC'
test_group='LIHC'
# target_cell = ["c13_Early_Erythroid_GATA1"]
target_cell_col = 're_annotation_2'


# In[10]:


adata.obs['re_annotation_2'].unique()


# In[181]:


for i in adata.obs['re_annotation_2'].unique():
    for cs in cancer_site:
        path = f"./results/lineage_All_LIHC_vs_KIRC_{cs}"
        os.makedirs(path, exist_ok=True)  
        test_adata = adata[(adata.obs[target_cell_col].isin([i])) & (adata.obs["cancer_site"].isin([cs]))].copy()
        dds=ov.bulk.pyDEG(test_adata.to_df(layer='lognorm').T)
        dds.drop_duplicates_index()
        print('... drop_duplicates_index success')
        treatment_groups=test_adata.obs[test_adata.obs[condition]==test_group].index.tolist()
        control_groups=test_adata.obs[test_adata.obs[condition]==ctrl_group].index.tolist()
        result=dds.deg_analysis(treatment_groups,control_groups,method='ttest')
        # -1 means automatically calculates
        dds.foldchange_set(fc_threshold=-1,
                       pval_threshold=0.05,
                       logp_max=10)
        dds.plot_volcano(title=f'DEG Analysis of {i}',figsize=(8,8),
                     plot_genes_num=50,plot_genes_fontsize=12,)
        i = i.replace("/", "_")
        result.sort_values('qvalue').to_csv(f"{path}/ttest_{i}.csv", index=True)
       
