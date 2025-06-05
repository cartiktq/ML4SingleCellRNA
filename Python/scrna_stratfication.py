import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import os

# --- CONFIGURATION ---
INPUT_DIR = "./output"
EXPR_FILE = os.path.join(INPUT_DIR, "expression_matrix.csv")
META_FILE = os.path.join(INPUT_DIR, "metadata.csv")
OUTPUT_FILE = os.path.join(INPUT_DIR, "metadata_with_clusters.csv")

# --- LOAD DATA ---
print("Loading data...")
expr_df = pd.read_csv(EXPR_FILE, index_col=0)
meta_df = pd.read_csv(META_FILE, index_col=0)

# --- STANDARDIZE AND REDUCE DIMENSIONS ---
print("Standardizing and applying PCA...")
scaler = StandardScaler()
X_scaled = scaler.fit_transform(expr_df)

# Use PCA to reduce dimensions to retain 95% variance or fixed components
pca = PCA(n_components=50)
X_pca = pca.fit_transform(X_scaled)

# --- CLUSTERING ---
# Try different values of k or use silhouette analysis for optimal k
k = 5
print(f"Clustering cells into {k} groups using KMeans...")
kmeans = KMeans(n_clusters=k, random_state=42)
clusters = kmeans.fit_predict(X_pca)

# --- ATTACH CLUSTER LABELS TO METADATA ---
meta_df["cluster"] = clusters

# --- SAVE UPDATED METADATA ---
print(f"Saving stratified metadata to {OUTPUT_FILE}...")
meta_df.to_csv(OUTPUT_FILE)

print("Stratification complete.")
