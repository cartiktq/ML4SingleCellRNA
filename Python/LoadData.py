import scanpy as sc
import os

class LoadData(WorkflowStep):
    def __init__(self, output_path, dataset_name="pbmc3k", file_format="h5ad"):
        super().__init__(input_path=None, output_path=output_path)
        self.dataset_name = dataset_name
        self.file_format = file_format

    def run(self):
        print(f"Loading dataset: {self.dataset_name}")
        adata = sc.datasets.pbmc3k() if self.dataset_name == "pbmc3k" else sc.read(self.dataset_name)
        os.makedirs(self.output_path, exist_ok=True)
        adata.write_h5ad(os.path.join(self.output_path, f"adata_raw.{self.file_format}"))
