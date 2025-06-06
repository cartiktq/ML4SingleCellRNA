class TrajectoryAnalysis(WorkflowStep):
    def run(self):
        import scanpy as sc

        adata = sc.read_h5ad(os.path.join(self.input_path, "adata_markers.h5ad"))
        sc.tl.diffmap(adata)
        sc.tl.dpt(adata, n_dcs=10)
        adata.write_h5ad(os.path.join(self.output_path, "adata_pseudotime.h5ad"))
