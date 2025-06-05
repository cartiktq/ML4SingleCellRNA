import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# --- CONFIGURATION ---
INPUT_DIR = "./output"
EXPR_FILE = os.path.join(INPUT_DIR, "expression_matrix.csv")
META_FILE = os.path.join(INPUT_DIR, "metadata.csv")
PLOT_DIR = "./plots"
os.makedirs(PLOT_DIR, exist_ok=True)

# --- LOAD DATA ---
print("Loading expression matrix and metadata...")
expr_df = pd.read_csv(EXPR_FILE, index_col=0)
meta_df = pd.read_csv(META_FILE, index_col=0)

# --- SUMMARY STATISTICS ---
print("Generating summary statistics...")

summary_stats = {
    "Number of cells": expr_df.shape[0],
    "Number of genes": expr_df.shape[1],
    "Mean total counts per cell": expr_df.sum(axis=1).mean(),
    "Median total counts per cell": expr_df.sum(axis=1).median(),
    "Max counts in a single gene across all cells": expr_df.max().max()
}

summary_df = pd.DataFrame.from_dict(summary_stats, orient='index', columns=["Value"])
print(summary_df)

# --- VISUALIZATION FUNCTIONS ---
def plot_total_counts_per_cell(df, save_path):
    total_counts = df.sum(axis=1)
    plt.figure(figsize=(8, 5))
    sns.histplot(total_counts, bins=50, kde=True)
    plt.title("Total Counts per Cell")
    plt.xlabel("Total Counts")
    plt.ylabel("Number of Cells")
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

def plot_top_variable_genes(df, top_n, save_path):
    gene_variance = df.var(axis=0).sort_values(ascending=False).head(top_n)
    plt.figure(figsize=(10, 6))
    sns.barplot(x=gene_variance.index, y=gene_variance.values, palette="viridis")
    plt.xticks(rotation=90)
    plt.title(f"Top {top_n} Most Variable Genes")
    plt.ylabel("Variance")
    plt.xlabel("Gene")
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

def plot_mitochondrial_content(meta_df, save_path):
    if "pct_counts_mt" in meta_df.columns:
        plt.figure(figsize=(8, 5))
        sns.histplot(meta_df["pct_counts_mt"], bins=40, color="tomato", kde=True)
        plt.title("Mitochondrial Gene Expression Percentage")
        plt.xlabel("% Mitochondrial Counts")
        plt.ylabel("Number of Cells")
        plt.tight_layout()
        plt.savefig(save_path)
        plt.close()

# --- GENERATE PLOTS ---
print("Generating visualizations...")

plot_total_counts_per_cell(expr_df, os.path.join(PLOT_DIR, "total_counts_per_cell.png"))
plot_top_variable_genes(expr_df, top_n=20, save_path=os.path.join(PLOT_DIR, "top_variable_genes.png"))
plot_mitochondrial_content(meta_df, save_path=os.path.join(PLOT_DIR, "mitochondrial_content.png"))

print(f"Summary and visualizations saved in '{PLOT_DIR}' directory.")
