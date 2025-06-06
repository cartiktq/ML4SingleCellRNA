class BiomarkerIdentification(WorkflowStep):
    def run(self):
        import scanpy as sc

        adata = sc.read_h5ad(os.path.join(self.input_path, "adata_clustered.h5ad"))
        sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test')
        adata.write_h5ad(os.path.join(self.output_path, "adata_markers.h5ad"))
