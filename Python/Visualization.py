class Visualization(WorkflowStep):
    def run(self):
        import scanpy as sc
        import os

        adata = sc.read_h5ad(os.path.join(self.input_path, "adata_pseudotime.h5ad"))
        sc.pl.umap(adata, color=['leiden'], save="_clusters.png", show=False)
        sc.pl.umap(adata, color=['dpt_pseudotime'], save="_pseudotime.png", show=False)

        top_genes = adata.uns['rank_genes_groups']['names'][0][:5]
        sc.pl.violin(adata, top_genes, groupby='leiden', save="_markers.png", show=False)
        sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, groupby='leiden', save="_heatmap.png", show=False)

        adata.obs.to_csv(os.path.join(self.output_path, "cell_metadata.csv"))
