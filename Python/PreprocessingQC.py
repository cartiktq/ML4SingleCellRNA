class PreprocessingQC(WorkflowStep):
    def run(self):
        import scanpy as sc
        import os

        adata = sc.read_h5ad(os.path.join(self.input_path, "adata_raw.h5ad"))
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
        adata = adata[adata.obs.n_genes_by_counts < 2500, :]
        adata = adata[adata.obs.pct_counts_mt < 5, :]
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata = adata[:, adata.var.highly_variable]
        sc.pp.scale(adata, max_value=10)

        os.makedirs(self.output_path, exist_ok=True)
        adata.write_h5ad(os.path.join(self.output_path, "adata_qc.h5ad"))
