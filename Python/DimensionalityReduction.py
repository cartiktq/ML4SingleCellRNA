class DimensionalityReduction(WorkflowStep):
    def run(self):
        import scanpy as sc

        adata = sc.read_h5ad(os.path.join(self.input_path, "adata_qc.h5ad"))
        sc.tl.pca(adata, svd_solver='arpack')
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        sc.tl.umap(adata)
        adata.write_h5ad(os.path.join(self.output_path, "adata_reduced.h5ad"))
