class ClusteringStratification(WorkflowStep):
    def run(self):
        import scanpy as sc

        adata = sc.read_h5ad(os.path.join(self.input_path, "adata_reduced.h5ad"))
        sc.tl.leiden(adata, resolution=0.5)
        adata.write_h5ad(os.path.join(self.output_path, "adata_clustered.h5ad"))
