
# GENE EXPRESSION CLUSTERING USING FUZZY ENSEMBLE LEARNING


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -------------------- SKLEARN --------------------
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans, AgglomerativeClustering, SpectralClustering
from sklearn.metrics import (
    silhouette_score,
    silhouette_samples,
    davies_bouldin_score,
    adjusted_rand_score,
)
from sklearn.impute import KNNImputer
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_classif
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier

# -------------------- OTHER --------------------
from skfuzzy.cluster import cmeans
import pyreadr


# 1. DATA LOADING
#    Reads an .rda file and extracts the Golub leukemia dataset.
#    Transposes if necessary so rows = samples, cols = genes.


def load_data(path):
    print("\n[INFO] Loading data from:", path)

    result = pyreadr.read_r(path)
    print("[INFO] Objects found:", list(result.keys()))

    # Extract expression matrix and class labels
    X = result['golub_train_3051']
    y = result['golub_train_response'].values.ravel()

    # Keep numeric columns only; drop columns with all-NaN
    X = X.select_dtypes(include=[np.number]).dropna(axis=1)

    # Ensure orientation: rows = samples, cols = genes
    if X.shape[0] > X.shape[1]:
        X = X.T

    print(f"[INFO] Final data shape: {X.shape}")
    return X.values, y



# 2. GENE SELECTION
#    Uses ANOVA F-score (f_classif) to select the top-k most
#    informative genes w.r.t. the class labels.
#    Reduces noise and dimensionality before preprocessing.


def select_genes(X, y, k=500):
    print(f"[INFO] Selecting top {k} genes via f_classif...")
    selector = SelectKBest(score_func=f_classif, k=k)
    X_new = selector.fit_transform(X, y)
    print(f"[INFO] Shape after selection: {X_new.shape}")
    return X_new



# 3. PREPROCESSING
#    Pipeline:
#      1) Shift negative values to >= 0  [FIX from final version]
#      2) log1p transform (variance stabilisation for expression data)
#      3) KNN imputation for any remaining NaNs
#      4) Z-score standardisation
#      5) Remove near-zero variance features


def preprocess(X):
    print("[INFO] Preprocessing...")

    # Step 1 — Shift: log1p is only valid for x >= 0
    min_val = np.min(X)
    if min_val < 0:
        X = X - min_val  # shift entire matrix up

    # Step 2 — Log transform: compresses dynamic range
    X = np.log1p(X)

    # Step 3 — Imputation: fills any residual NaN values
    X = KNNImputer(n_neighbors=5).fit_transform(X)

    # Step 4 — Standardise: zero mean, unit variance
    X = StandardScaler().fit_transform(X)

    # Step 5 — Remove near-constant features (threshold = 0.01)
    X = VarianceThreshold(0.01).fit_transform(X)

    print(f"[INFO] Shape after preprocessing: {X.shape}")
    return X



# 4. DIMENSIONALITY REDUCTION (PCA)
#    Retains up to 50 principal components (or min(shape)-1 if
#    the dataset is smaller) to feed compact features to
#    downstream clustering.


def reduce_dim(X):
    n_comp = min(50, min(X.shape) - 1)
    print(f"[INFO] PCA: retaining {n_comp} components")
    return PCA(n_components=n_comp).fit_transform(X)


# 5. INDIVIDUAL CLUSTERING METHODS
#    Used both independently (ablation) and as base learners
#    for the ensemble.


def run_kmeans(X, k):
    """Standard KMeans clustering."""
    return KMeans(n_clusters=k, random_state=42, n_init=10).fit_predict(X)


def run_hierarchical(X, k):
    """Agglomerative (Ward) hierarchical clustering."""
    return AgglomerativeClustering(n_clusters=k).fit_predict(X)


def run_spectral(X, k):
    """
    Spectral clustering with nearest-neighbour affinity.
    n_neighbors is clamped to (n_samples - 1) to avoid errors
    on small datasets.  [FIX from final version]
    """
    n_neighbors = min(10, X.shape[0] - 1)
    return SpectralClustering(
        n_clusters=k,
        affinity='nearest_neighbors',
        n_neighbors=n_neighbors,
        random_state=42
    ).fit_predict(X)


# 6. FUZZY C-MEANS
#    Returns the full soft-membership matrix U (shape: k x n).
#    m=2 is the standard fuzziness exponent.


def run_fuzzy_cmeans(X, k):
    """Apply Fuzzy C-Means; returns membership matrix U."""
    _cntr, u, _u0, _d, _jm, _p, _fpc = cmeans(
        X.T,          # skfuzzy expects (features x samples)
        c=k,
        m=2,
        error=0.005,
        maxiter=1000
    )
    return u  # shape: (k, n_samples)


# 7. CO-ASSOCIATION MATRIX
#    Counts how often each pair of samples is assigned to the
#    same cluster across all base clustering runs, then
#    normalises to [0, 1].  Higher value = more co-clustered.


def coassociation(clusterings):
    """
    Build a co-association (consensus) matrix from a list of
    hard label arrays.  Entry [i, j] = fraction of clusterings
    in which samples i and j are in the same cluster.
    """
    n = len(clusterings[0])
    mat = np.zeros((n, n))

    for labels in clusterings:
        for i in range(n):
            for j in range(n):
                if labels[i] == labels[j]:
                    mat[i][j] += 1

    return mat / len(clusterings)  # normalise to [0, 1]


# 8. FUZZY ENSEMBLE
#    Strategy:
#      - Run KMeans 5× (different random inits) for stability
#        [improvement from final version]
#      - Add Hierarchical + Spectral as diverse base learners
#      - Build co-association matrix from all 7 runs
#      - Normalise via MinMaxScaler  [FIX from final version]
#      - Apply Fuzzy C-Means on the consensus matrix
#      - Return hard labels via argmax over membership matrix


def fuzzy_ensemble(X, k):
    """
    Ensemble clustering combining KMeans (×5), Hierarchical,
    and Spectral via a fuzzy co-association approach.
    """
    print("[INFO] Building fuzzy ensemble...")

    # Base learners: 5 KMeans runs + hierarchical + spectral
    base_clusterings = [run_kmeans(X, k) for _ in range(5)]
    base_clusterings.append(run_hierarchical(X, k))
    base_clusterings.append(run_spectral(X, k))

    # Co-association matrix (soft similarity between samples)
    co_mat = coassociation(base_clusterings)

    # Normalise similarity values to [0, 1] for stable FCM
    co_mat = MinMaxScaler().fit_transform(co_mat)

    # Apply Fuzzy C-Means on the consensus similarity space
    u = run_fuzzy_cmeans(co_mat, k)

    # Defuzzify: assign each sample to its highest-membership cluster
    return np.argmax(u, axis=0)


# 9. EVALUATION METRICS
#    Silhouette Score  : higher is better  (range: -1 to 1)
#    Davies-Bouldin    : lower is better   (range: >= 0)
#    Adjusted Rand Index: higher is better (range: -1 to 1)


def evaluate(X, labels, y):
    """Return a dict of three external/internal clustering metrics."""
    return {
        "Silhouette":      silhouette_score(X, labels),
        "Davies-Bouldin":  davies_bouldin_score(X, labels),
        "ARI":             adjusted_rand_score(y, labels),
    }



# 10. BEST K SELECTION
#     Sweeps k in [2, 4] using KMeans + Silhouette score to
#     pick the number of clusters before running the full
#     ensemble (faster proxy search).


def find_best_k(X, k_range=(2, 4)):
    print("\n[INFO] Searching for best K...")
    best_k, best_score = 2, -1

    for k in range(k_range[0], k_range[1] + 1):
        labels = run_kmeans(X, k)
        score = silhouette_score(X, labels)
        print(f"  K={k}  Silhouette={score:.4f}")

        if score > best_score:
            best_k, best_score = k, score

    print(f"[INFO] Best K = {best_k} (Silhouette = {best_score:.4f})")
    return best_k


# 11. CLASSIFIER VALIDATION
#     Uses cluster labels as pseudo-class labels and trains
#     three classifiers to assess cluster separability.
#     High accuracy indicates well-separated, consistent clusters.


def classifier_test(X, labels):
    """Train RF, SVM, KNN on cluster labels; report test accuracy."""
    X_tr, X_te, y_tr, y_te = train_test_split(
        X, labels, test_size=0.3, random_state=42
    )

    models = {
        "RandomForest": RandomForestClassifier(random_state=42),
        "SVM":          SVC(),
        "KNN":          KNeighborsClassifier(),
    }

    print("\n[INFO] Classifier Validation (cluster label prediction):")
    for name, model in models.items():
        model.fit(X_tr, y_tr)
        acc = model.score(X_te, y_te)
        print(f"  {name:>14}: {acc:.4f}")



# 12. ABLATION STUDY
#     Measures the contribution of each component by comparing:
#       - KMeans alone
#       - Fuzzy C-Means alone
#       - Ensemble without Fuzzy (KMeans on co-assoc matrix)
#       - Full Fuzzy Ensemble (proposed method)
#     [Retained from original pipeline — key for paper/report]


def ablation(X, k, y):
    """
    Compare four configurations to assess each component's value.
    Returns a dict mapping method name -> evaluation dict.
    """
    results = {}

    # Configuration 1: Plain KMeans
    results["KMeans"] = evaluate(X, run_kmeans(X, k), y)

    # Configuration 2: Fuzzy C-Means only (no ensemble)
    u_fcm = run_fuzzy_cmeans(X, k)
    results["Fuzzy C-Means"] = evaluate(X, np.argmax(u_fcm, axis=0), y)

    # Configuration 3: Ensemble without Fuzzy
    #   — co-association built from KMeans + Hierarchical + Spectral,
    #     then hard KMeans applied on the consensus matrix
    c1 = run_kmeans(X, k)
    c2 = run_hierarchical(X, k)
    c3 = run_spectral(X, k)
    co_mat = coassociation([c1, c2, c3])
    co_mat_norm = MinMaxScaler().fit_transform(co_mat)
    results["Ensemble (no Fuzzy)"] = evaluate(X, run_kmeans(co_mat_norm, k), y)

    # Configuration 4: Full Fuzzy Ensemble (proposed)
    results["Fuzzy Ensemble"] = evaluate(X, fuzzy_ensemble(X, k), y)

    return results


def plot_ablation(results):
    """Bar chart comparing Silhouette scores across ablation configurations."""
    methods = list(results.keys())
    scores  = [results[m]["Silhouette"] for m in methods]

    plt.figure(figsize=(8, 5))
    bars = plt.bar(methods, scores, color=["#4C72B0", "#55A868", "#C44E52", "#8172B2"])
    plt.xticks(rotation=20, ha='right')
    plt.ylabel("Silhouette Score")
    plt.title("Ablation Study — Silhouette Score by Method")
    plt.tight_layout()

    # Annotate bar values
    for bar, score in zip(bars, scores):
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.005,
            f"{score:.3f}",
            ha='center', va='bottom', fontsize=9
        )

    plt.savefig("ablation_study.png", dpi=150, bbox_inches='tight')
    plt.show()
    print("[INFO] Ablation plot saved to ablation_study.png")


# 13. VISUALISATION
#     Three diagnostic plots:
#       a) Heatmap — expression patterns per cluster
#       b) t-SNE   — 2D cluster topology
#       c) Silhouette distribution — per-sample cohesion


def plot_heatmap(X, labels):
    """Sorted gene expression heatmap (samples grouped by cluster)."""
    idx = np.argsort(labels)
    plt.figure(figsize=(10, 6))
    plt.imshow(X[idx], aspect='auto', interpolation='nearest')
    plt.colorbar(label="Scaled Expression")
    plt.xlabel("Genes (PCs)")
    plt.ylabel("Samples (sorted by cluster)")
    plt.title("Gene Expression Heatmap — Sorted by Cluster")
    plt.tight_layout()
    plt.savefig("heatmap.png", dpi=150, bbox_inches='tight')
    plt.show()
    print("[INFO] Heatmap saved to heatmap.png")


def plot_tsne(X, labels):
    """t-SNE embedding coloured by cluster assignment."""
    print("[INFO] Running t-SNE (may take a moment)...")
    emb = TSNE(n_components=2, random_state=42).fit_transform(X)

    plt.figure(figsize=(7, 6))
    scatter = plt.scatter(emb[:, 0], emb[:, 1], c=labels, cmap='tab10', s=40)
    plt.colorbar(scatter, label="Cluster")
    plt.title("t-SNE Cluster Visualisation")
    plt.xlabel("t-SNE 1")
    plt.ylabel("t-SNE 2")
    plt.tight_layout()
    plt.savefig("tsne.png", dpi=150, bbox_inches='tight')
    plt.show()
    print("[INFO] t-SNE plot saved to tsne.png")


def plot_silhouette(X, labels):
    """Histogram of per-sample silhouette values."""
    vals = silhouette_samples(X, labels)
    plt.figure(figsize=(7, 4))
    plt.hist(vals, bins=20, color="#4C72B0", edgecolor="white")
    plt.axvline(vals.mean(), color='red', linestyle='--',
                label=f"Mean = {vals.mean():.3f}")
    plt.xlabel("Silhouette Value")
    plt.ylabel("Count")
    plt.title("Per-Sample Silhouette Distribution")
    plt.legend()
    plt.tight_layout()
    plt.savefig("silhouette_dist.png", dpi=150, bbox_inches='tight')
    plt.show()
    print("[INFO] Silhouette distribution saved to silhouette_dist.png")

# 14. MAIN PIPELINE
#     Orchestrates all steps in order and prints a summary.


def run(path):
    # --- Data ---
    X, y = load_data(path)

    # --- Feature engineering ---
    X = select_genes(X, y, k=500)
    X = preprocess(X)
    X = reduce_dim(X)

    # --- Optimal K ---
    k = find_best_k(X, k_range=(2, 4))

    # --- Final model: Fuzzy Ensemble ---
    labels = fuzzy_ensemble(X, k)

    # --- Evaluation ---
    print("\n[FINAL RESULTS — Fuzzy Ensemble]")
    metrics = evaluate(X, labels, y)
    for metric, val in metrics.items():
        print(f"  {metric:>16}: {val:.4f}")

    # --- Classifier validation ---
    classifier_test(X, labels)

    # --- Ablation study ---
    print("\n[ABLATION STUDY]")
    ab_results = ablation(X, k, y)
    for method, res in ab_results.items():
        sil = res["Silhouette"]
        db  = res["Davies-Bouldin"]
        ari = res["ARI"]
        print(f"  {method:<22} | Silhouette={sil:.4f}  DB={db:.4f}  ARI={ari:.4f}")

    plot_ablation(ab_results)

    # --- Visualisation ---
    plot_heatmap(X, labels)
    plot_tsne(X, labels)
    plot_silhouette(X, labels)


# ENTRY POINT


if __name__ == "__main__":
    # Update this path to point to your local copy of the dataset
    DATA_PATH = r"C:\Users\BIT\Downloads\archive (6)\leukemia_data_Golub99_3051.rda"
    run(DATA_PATH)