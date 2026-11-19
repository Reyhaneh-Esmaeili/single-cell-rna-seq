# Single-Cell RNA-seq Analysis: Clustering Comparison and Differential Expression

A Scanpy workflow for single-cell RNA-seq data covering quality control, normalization, dimensionality reduction, comparison of three unsupervised clustering methods, and differential expression analysis of selected clusters.

## Dataset

UMI count matrix of 7,551 cells with cell-level metadata, including reference cell-type annotations (`celltype`, `redefined_cluster`). The data are human hematopoietic and immune cells (T cells, B cells, plasma cells and progenitors such as HSC, MPP, CMP, GMP and MEP).


The notebook reads the data from `data/single_cell_dataset.zip`, which must contain:

- `ABC_umi_matrix_7551_cells.csv`: cells x genes UMI counts
- `ABC_Meta.txt`: tab-separated cell metadata

The archive is extracted automatically to `data/extracted/` the first time the notebook runs.

## Workflow

1. **Preprocessing**: merge matrix and metadata into an AnnData object; QC on detected genes, total UMIs and mitochondrial fraction; cell filtering (1,000-7,000 genes, 2,000-200,000 UMIs); 7,070 cells passed QC
2. **Normalization**: library-size normalization to 10,000 counts per cell, `log1p`, 2,000 highly variable genes (Seurat method), scaling (clipped at ±10)
3. **Clustering**: PCA (30 PCs), then K-means, Ward hierarchical and spectral clustering, each with 15 clusters
4. **Evaluation**: silhouette score (in PCA space), ARI and NMI against the reference annotations
5. **Visualization**: UMAP and t-SNE colored by spectral cluster and reference cell type
6. **Differential expression**: Wilcoxon rank-sum test (cluster vs. rest) for three selected clusters, with interpretation of top genes

## Results

| Method | Silhouette | ARI | NMI |
|---|---|---|---|
| K-means | 0.330 | 0.298 | 0.650 |
| Hierarchical (Ward) | 0.305 | 0.262 | 0.636 |
| Spectral | 0.312 | 0.309 | 0.655 |

Spectral clustering had the highest ARI and NMI against the reference annotations, although the margin over K-means is small, and K-means had the highest silhouette score. Spectral clusters were used for downstream analysis.

Top differentially expressed genes of the selected clusters:

| Cluster | Dominant cell types | Top genes |
|---|---|---|
| 2 | CD4+ / CD8+ T cells | EEF1A1, IL32, RPL13, RPLP2, CD3D |
| 5 | Naive / regulatory B cells | MS4A1, CD74, BANK1, CD79A, CD52 |
| 10 | Plasma cells | SSR4, MZB1, DERL3, FKBP11, SEC11C |

The marker genes match the expected identities of each cell population. Hematopoietic progenitors form a continuous region in UMAP and t-SNE rather than discrete clusters.

## Repository structure

```
.
├── scrna_clustering_de_analysis.ipynb   # full analysis with outputs and interpretation
├── data/
│   └── single_cell_dataset.zip          # input data (see Dataset)
├── requirements.txt
├── LICENSE
└── README.md
```

## Running the notebook

```bash
git clone https://github.com/<your-username>/<repo-name>.git
cd <repo-name>
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
jupyter notebook scrna_clustering_de_analysis.ipynb
```

Developed with Python 3.10.

## Notes

- Mitochondrial genes (`MT-` prefix) were not present in the matrix, so mitochondrial filtering was not informative and no cells were removed on that basis.
- Random seeds are fixed (`random_state=42`) for K-means, spectral clustering and t-SNE.

## License

The code and notebook are released under the [MIT License](LICENSE). The dataset is not covered by this license; see its original source for terms of use.
