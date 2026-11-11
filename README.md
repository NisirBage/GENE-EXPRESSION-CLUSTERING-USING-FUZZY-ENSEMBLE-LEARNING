# Gene-Expression Clustering with a Fuzzy Ensemble

Unsupervised clustering of leukaemia gene-expression profiles (ALL vs AML)
using a consensus ensemble: five k-means runs, agglomerative clustering and
spectral clustering are combined into a co-association matrix, and fuzzy
c-means on that matrix yields the final clusters.

## Data

- **Golub et al. (1999) leukaemia training set**: 38 samples × 3,051 genes
  (Affymetrix), loaded from the R object `leukemia_data_Golub99_3051.rda`
  (objects `golub_train_3051` and `golub_train_response`). The `.rda` file is
  not in this repository; `converted_gene_expression.csv` holds the expression
  matrix only, without labels.
- `pipeline_unified_3.py` additionally pulls five GEO series (GSE297413,
  GSE253086, GSE201492, GSE216738, GSE181157) from NCBI and harmonises them.

## Pipeline (`fullpipelinefinal.py`)

1. Select the top 500 genes by ANOVA F-score (`SelectKBest(f_classif)`).
2. Shift to non-negative, `log1p`, KNN imputation, z-scoring, and removal of
   near-constant features.
3. PCA to at most 50 components.
4. Choose k in {2, 3, 4} by silhouette score.
5. Fuzzy ensemble: build the co-association matrix from the base clusterings,
   min–max normalise it, run fuzzy c-means, and take the argmax membership.
6. Evaluate with silhouette, Davies–Bouldin and adjusted Rand index (ARI)
   against the ALL/AML labels. Run an ablation against single methods.

## Results

Silhouette score by method, from `ablation_study.png`:

| Method | Silhouette |
|---|---:|
| k-means | 0.272 |
| Fuzzy c-means | 0.248 |
| Ensemble (no fuzzy step) | 0.272 |
| Fuzzy ensemble | 0.272 |

On this dataset the fuzzy ensemble **does not improve** on plain k-means by
silhouette. Its ARI values are printed by the script; they are not recorded
in this repository.

## Limitations

- **Label leakage into an unsupervised evaluation.** Step 1 selects genes
  using the class labels, and step 6 scores clusters against those same
  labels. ARI is therefore optimistic. A fair evaluation would select genes
  without labels (e.g. by variance) or inside a held-out split.
- **n = 38.** Differences in the third decimal place are not meaningful at
  this sample size; no confidence intervals or repeated-seed statistics are
  reported.
- The data path is hard-coded in the `__main__` block of both scripts.

## Running

```bash
pip install -r requirements.txt
# edit DATA_PATH at the bottom of fullpipelinefinal.py, then:
python fullpipelinefinal.py
```
