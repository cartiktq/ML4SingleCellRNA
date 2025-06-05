import scanpy as sc
import pandas as pd
import os

# --- CONFIGURATION PARAMETERS (OBFUSCATED) ---
DATASET_ID = "GSE123456"  # Replace with real accession
OUTPUT_DIR = "./output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- STEP A: EXTRACT DATA FROM PUBLIC REPOSITORY ---
def download_scrna_data(dataset_id: str) -> sc.AnnData:
    """
    Download scRNA-seq dataset using GEO ID via scanpy's built-in functions.
    """
    print(f"Downloading dataset {dataset_id}...")
    adata = sc.datasets.pbmc3k()  # Substitute with actual download logic if needed
    return adata

# --- STEP B: QUALITY CONTROL AND PREPROCESSING ---
def preprocess_data(adata: sc.AnnData) -> sc.AnnData:
    """
    Performs standard QC and preprocessing:
    - Filters cells and genes
    - Normalizes and logs the data
    """
    print("Starting preprocessing...")

    # Basic QC metrics
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # Filter cells and genes
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs.pct_counts_mt < 5, :]  # Remove cells with >5% mitochondrial counts

    # Normalize and log-transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    print("QC and normalization complete.")
    return adata

# --- STEP C: OUTPUT TO CSV ---
def export_to_csv(adata: sc.AnnData, output_dir: str) -> None:
    """
    Exports gene expression matrix and metadata to CSV.
    """
    print(f"Exporting data to {output_dir}...")

    # Expression matrix
    expr_df = pd.DataFrame(adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
                           index=adata.obs_names, columns=adata.var_names)
    expr_df.to_csv(os.path.join(output_dir, "expression_matrix.csv"))

    # Metadata
    adata.obs.to_csv(os.path.join(output_dir, "metadata.csv"))

    print("Export complete.")

# --- MAIN WORKFLOW ---
if __name__ == "__main__":
    adata = download_scrna_data(DATASET_ID)
    adata = preprocess_data(adata)
    export_to_csv(adata, OUTPUT_DIR)
