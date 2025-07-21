import argparse
import scanpy as sc
import pandas as pd
import numpy as np

def compute_top_markers(adata, cluster_col, top_n=10):
    sc.tl.rank_genes_groups(adata, groupby=cluster_col, method='wilcoxon')
    result = adata.uns['rank_genes_groups']
    
    clusters = result['names'].dtype.names
    markers = {}
    
    for cluster in clusters:
        df = pd.DataFrame({
            'gene': result['names'][cluster],
            'pval': result['pvals'][cluster],
            'logfc': result['logfoldchanges'][cluster]
        })
        df = df.dropna()
        df = df[df['logfc'] > 0]  # optional: positive logFC only
        df = df.sort_values(by=['pval', 'logfc'], ascending=[True, False])
        top_genes = df['gene'].head(top_n).tolist()
        markers[cluster] = top_genes
    
    return markers

def format_prompt(markers, tissue_name):
    header = f"Identify cell types of {tissue_name} cells using the following markers separately for each row. Only provide the cell type name. Do not show numbers before the name. Some can be a mixture of multiple cell types.\n"
    body = "\n".join([f"{cluster}:{','.join(genes)}" for cluster, genes in markers.items()])
    return header + body

def main():
    parser = argparse.ArgumentParser(description="Generate GPT prompt for cell type annotation from h5ad")
    parser.add_argument("h5ad_file", help="Path to .h5ad file")
    parser.add_argument("--cluster_col", required=True, help="Column in .obs for cluster labels")
    parser.add_argument("--top_n", type=int, default=10, help="Number of top genes per cluster")
    parser.add_argument("--tissue", default="the tissue", help="Tissue name for prompt")
    
    args = parser.parse_args()

    # Load data
    adata = sc.read(args.h5ad_file)

    # Compute markers
    markers = compute_top_markers(adata, cluster_col=args.cluster_col, top_n=args.top_n)

    # Create prompt
    prompt = format_prompt(markers, tissue_name=args.tissue)

    # Print prompt
    print("\n=== COPY THIS PROMPT INTO CHATGPT ===\n")
    print(prompt)
    print("\n=== END OF PROMPT ===\n")

if __name__ == "__main__":
    main()
