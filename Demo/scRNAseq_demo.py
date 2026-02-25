import matplotlib
matplotlib.use('Agg') 

import scanpy as sc
import anndata as ad
import pandas as pd
import scanpy.external as sce
import matplotlib.pyplot as plt
import seaborn as sns
import time
import warnings
import os

warnings.filterwarnings('ignore')
# ==============================================================================
# 1. Define File Information
# ==============================================================================
files_info = [
    {"path": "./Demodata/subset_1pct_CAR_NK.h5ad", "group": "CAR_NK"},
    {"path": "./Demodata/subset_1pct_START.h5ad", "group": "START"},
    {"path": "./Demodata/subset_1pct_ECART.h5ad", "group": "ECART"}
]

# ==============================================================================
# 2. Loop through reading, processing, and filtering data
# ==============================================================================
adatas = []

for info in files_info:
    print(f"Reading: {info['group']} ...")
    adata = sc.read_h5ad(info['path'])
    adatas.append(adata)

# ==============================================================================
# 3. Concatenate (Merge)
# ==============================================================================
print("\nMerging objects...")
adata_B = ad.concat(adatas, join='outer', label='batch', keys=[f['group'] for f in files_info], index_unique='-', fill_value=0)
print(f"Merge complete: {adata_B.shape}")


import time
import scanpy.external as sce

adata_B.raw = adata_B.copy()
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),' - Start highly_variable_genes')
sc.pp.highly_variable_genes(adata_B, n_top_genes=3000, batch_key='sample', subset=False)
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),' - Start regress_out')
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),' - Start scale')
sc.pp.scale(adata_B, max_value=10)

print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),' - Start pca')
sc.tl.pca(adata_B, svd_solver="arpack", use_highly_variable=True)

print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),' - Start harmony')
sce.pp.harmony_integrate(adata_B, ['Sample_Group','sample_id'],sigma=0.6,**{'max_iter_harmony' : 30})
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),' - Start neighbors')
sc.pp.neighbors(adata_B,  use_rep='X_pca_harmony', n_neighbors=15)

print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),' - Start leiden')
sc.tl.leiden(adata_B,key_added='leiden_B')
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),' - Start louvain')
sc.tl.louvain(adata_B,key_added='louvain_B')
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),' - Start umap')
sc.tl.umap(adata_B)
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),' - Start end')

# ==============================================================================
# 4.  UMAP PDF (本地 Windows 优化版)
# ==============================================================================
print(time.strftime("%H:%M:%S"), ' - Plotting and Saving PDF')

# 解决 PDF 字体可编辑性问题
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['figure.figsize'] = (10, 8)
sns.set(style="ticks", font_scale=1.2)

# 确保保存目录存在
output_path = "./Demo/UMAP_Anonymized_Report.pdf"
os.makedirs(os.path.dirname(output_path), exist_ok=True)

# 绘图逻辑
# 注意：show=False 非常关键，配合 Agg 后端使用
fig = sc.pl.umap(
    adata_B, 
    color=['Sample_Group', 'sample'], 
    s=20, 
    ncols=2, 
    wspace=0.4,
    frameon=False,
    show=False, 
    return_fig=True
)

# 保存 PDF
fig.savefig(output_path, bbox_inches='tight', dpi=300)
plt.close(fig) # 释放内存

print(f"Success! PDF 已保存至: {os.path.abspath(output_path)}")

# 保存最终 h5ad
adata_B.write("./Demo/adata_B_final.h5ad")
print(time.strftime("%H:%M:%S"), ' - All Done.')