import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list

# === Data files ===
expression_file = "adult_genes_and_expression_value.xlsx"
orphan_gpcr_file = "Orphan_GPCR_Class_A.xlsx"

# Load gene expression data (samples as rows, genes as columns)
expr_df = pd.read_excel(expression_file, index_col=0)

# Load orphan GPCR gene list
orphan_gpcr_genes = pd.read_excel(orphan_gpcr_file, usecols=[0]).iloc[:, 0].dropna().unique()

# Filter expression data to include only orphan GPCR genes
orphanGPCR_expr_df = expr_df.loc[:, expr_df.columns.intersection(orphan_gpcr_genes)]

# === Organize samples ===
normal_samples = [
    "TCGA.06.0673", "TCGA.06.0675", "TCGA.06.0676",
    "TCGA.06.0678", "TCGA.06.0680", "TCGA.06.0681",
    "TCGA.08.0623", "TCGA.08.0625", "TCGA.08.0626", "TCGA.08.0627"
]

# Ensure sample IDs are strings and index is consistent
orphanGPCR_expr_df.index = orphanGPCR_expr_df.index.astype(str)

# Split samples into normal and tumor groups
normal_df = orphanGPCR_expr_df.loc[orphanGPCR_expr_df.index.isin(normal_samples)]
tumor_df = orphanGPCR_expr_df.loc[~orphanGPCR_expr_df.index.isin(normal_samples)]

# Transpose to genes x samples
normal_data = normal_df.T
tumor_data = tumor_df.T

# Cluster columns (samples) within each sample group
def cluster_columns(data):
    linkage_result = linkage(data.T, method='average')
    col_order = leaves_list(linkage_result)
    return data.iloc[:, col_order]

normal_data_clustered = cluster_columns(normal_data)
tumor_data_clustered = cluster_columns(tumor_data)

# Concatenate clustered columns
heatmap_data = pd.concat([normal_data_clustered, tumor_data_clustered], axis=1)

# === Plot the heatmap using clustermap ===
cg = sns.clustermap(
    heatmap_data,
    row_cluster=True,
    col_cluster=False,
    cmap="viridis",
    xticklabels=False,
    yticklabels=True,
    figsize=(14, 12),
    cbar_kws={"label": "mRNA Expression (log 2)"}
)

# Add a line separating normal and tumor samples
n_normal = normal_data_clustered.shape[1]
cg.ax_heatmap.axvline(n_normal, color='red', linestyle='--', linewidth=2)

# Finalizing the heatmap
cg.ax_heatmap.set_ylabel("Genes (Orphan GPCR Class A)")
cg.ax_heatmap.set_xlabel("Normal samples (left) vs GBM samples (right)")
plt.tight_layout()
cg.savefig("GBM_samples_orphan_gpcr_heatmap.png", dpi=300)
plt.show()
