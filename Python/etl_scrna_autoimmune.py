import os
import json
import scanpy as sc
import pandas as pd
import numpy as np
import GEOparse

# === CONFIGURATION === #
GEO_ID = "GSE123456"  # Replace with actual autoimmune-related scRNA dataset
WORKDIR = "./data"
OUTPUT_JSON = "preprocessed_scrna.json"

# === EXTRACT === #
def extract_geo_data(geo_id: str, workdir: str) -> str:
    os.makedirs(workdir, exist_ok=True)
    print(f"Fetching GEO dataset: {geo_id}")
    gse = GEOparse.get_GEO(geo=geo_id, destdir=workdir)
    # Assume supplementary files contain raw matrix or processed h5ad files
    # For example, download from supplementary URL or use scanpy to read directly
    data_path = os.path.join(workdir, "matrix.h5ad")  # placeholder
    return data_path

# === TRANSFORM === #
def preprocess_scrna(filepath: str) -> sc.AnnData:
    print(f"Reading scRNA-seq data from: {filepath}")
    adata = sc.read(filepath)

    # Basic quality control
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # Filter cells
    adata = adata[adata.obs["pct_counts_mt"] < 10, :]

    # Normalize & log transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Optional: high variable gene selection
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]

    return adata

# === LOAD === #
def save_to_json(adata: sc.AnnData, output_file: str):
    print(f"Saving preprocessed data to: {output_file}")

    # Downsample to avoid huge JSON (e.g., first 100 cells)
    n_cells = min(100, adata.n_obs)
    obs_subset = adata.obs.iloc[:n_cells].to_dict(orient="records")
    X_subset = adata.X[:n_cells].todense().tolist() if hasattr(adata.X, "todense") else adata.X[:n_cells].tolist()

    json_data = {
        "obs": obs_subset,
        "var_names": adata.var_names.tolist(),
        "X": X_subset
    }

    with open(output_file, "w") as f:
        json.dump(json_data, f, indent=2)

# === MAIN === #
if __name__ == "__main__":
    print("Starting ETL pipeline for scRNA-seq autoimmune data")

    # Extract
    data_path = extract_geo_data(GEO_ID, WORKDIR)

    # Transform
    adata = preprocess_scrna(data_path)

    # Load
    save_to_json(adata, OUTPUT_JSON)

    print("ETL step completed.")
