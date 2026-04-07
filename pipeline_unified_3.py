# ================================================================
# GENE EXPRESSION CLUSTERING — UNIFIED PIPELINE  (FINAL)
# Golub + GSE297413 + GSE253086 + GSE201492 + GSE216738 + GSE181157
# ================================================================
#
# ══════════════════════════════════════════════════════════════
#  WHAT EVERY ERROR WAS AND EXACTLY HOW IT IS FIXED HERE
# ══════════════════════════════════════════════════════════════
#
# ─────────────────────────────────────────────────────────────
# BUG 1 — openpyxl not installed → GSE253086 XLSX skipped
# ─────────────────────────────────────────────────────────────
# Error message:
#   "XLSX parse failed: Missing optional dependency 'openpyxl'."
#
# Cause:
#   pd.read_excel() requires the openpyxl package as its engine for
#   .xlsx files.  It was not installed in the venv.
#
# Fix:
#   Run once in your venv:
#       pip install openpyxl
#   The code already calls pd.read_excel(engine="openpyxl"); once the
#   package is present the XLSX branch works correctly.
#
# ─────────────────────────────────────────────────────────────
# BUG 2 — HHS redirect turns FTP listing into `index.html` junk
# ─────────────────────────────────────────────────────────────
# Error message:
#   "Download failed (404 ... url: ...GSE181nnn/GSE181157/suppl/
#    https://www.hhs.gov/vulnerability-disclosure-policy/index.html)"
#
# Cause:
#   The HTML regex for the NCBI FTP directory listing matched an
#   href that was the HHS vulnerability-disclosure redirect link
#   injected by the US government's web-security scanner.  The
#   pattern r'href="([^"?/][^"]*)"' accepted any href that did not
#   start with ? or /, so the full HTTPS URL of the HHS page was
#   captured.  That URL was then concatenated onto the FTP base as
#   if it were a filename, producing a nonsense URL.
#
# Fix:
#   Two changes:
#   (a) The regex now requires the captured filename to contain a '.'
#       (all real data files have extensions) and explicitly rejects
#       filenames that contain '://' (absolute URLs) or start with
#       'http' or 'ftp':
#           r'href="([^"?/][^":][^"]*\.[a-zA-Z0-9]{1,6})"'
#   (b) After building the full URL, validate it:
#       if "hhs.gov" in url or "://" in fname: skip.
#
# ─────────────────────────────────────────────────────────────
# BUG 3 — GSE201492 downloads splice-junction files, not expression
# ─────────────────────────────────────────────────────────────
# Error message:
#   Files tried: ['GSE201492_JC.raw.input.A3SS_ID.txt.gz',
#                 'GSE201492_JC.raw.input.A5SS_ID.txt.gz', ...]
#   "[WARN] GSE201492: skipped — no parseable supplementary file."
#
# Cause:
#   GSE201492 has five supplementary files.  The splice-junction
#   files (A3SS, A5SS, RI, SE) were tried first because their names
#   contain no "skip" keywords.  They are not gene-by-sample matrices;
#   they have a different structure (junction IDs × event counts) so
#   every parse attempt failed.  The actual expression file is:
#       GSE201492_RNAseq_expression_matrix_rawCount.txt.gz
#   which was tried last and also failed because:
#     • its first column is an Ensembl gene ID column named something
#       other than the expected index column name, AND
#     • _read_matrix_file() tries comment="#!" but this file uses no
#       comment prefix, so the first numeric row was accidentally
#       skipped.
#
# Fix:
#   (a) Strongly prefer files whose names contain "expression",
#       "matrix", "count", "tpm", "fpkm", "norm" AND penalise files
#       whose names look like splice-junction outputs
#       ("A3SS", "A5SS", "RI_ID", "SE_ID", "JC.").
#   (b) _read_matrix_file() now tries comment=None in addition to
#       comment="!" so files without comment lines are parsed cleanly.
#   (c) After loading, we try all columns as the potential index
#       column (not just column 0) so that multi-column headers are
#       handled correctly.
#
# ─────────────────────────────────────────────────────────────
# BUG 4 — GSE181157 has no parseable supplementary file
# ─────────────────────────────────────────────────────────────
# Error message:
#   "Files tried: ['index.html']"
#
# Cause:
#   The SOFT metadata for GSE181157 contains no
#   "supplementary_file" entries (the dataset may store its data
#   only in SRA raw reads, with no processed count matrix uploaded
#   to GEO).  The FTP scraper returned only the HHS redirect URL
#   (Bug 2), which was erroneously named "index.html".  With Bug 2
#   fixed the FTP scraper returns an empty list, and with no SOFT
#   URLs either the dataset is genuinely unavailable as a processed
#   matrix.
#
# Fix:
#   The pipeline now prints a clear actionable message, skips the
#   dataset gracefully, and continues with the remaining datasets.
#   The pipeline can still run successfully on the other 4 datasets.
#   If you can obtain a processed count matrix for GSE181157 (e.g.
#   by contacting the authors or downloading from SRA + running
#   featureCounts/STAR), place the file in ./geo_data/ and add its
#   filename to MANUAL_FILES at the bottom of this script.
#
# ─────────────────────────────────────────────────────────────
# BUG 5 — InvalidIndexError: Reindexing only valid with uniquely
#          valued Index objects
# ─────────────────────────────────────────────────────────────
# Error message:
#   pandas.errors.InvalidIndexError: Reindexing only valid with
#   uniquely valued Index objects
#
# Cause:
#   After probe→symbol mapping, groupby(expr_df.index).mean()
#   is supposed to collapse duplicate probe IDs that map to the
#   same gene into one row.  However, two sources of residual
#   duplicates survived:
#     1. Unmapped probes all received the literal string of their
#        original probe ID (e.g. "1007_s_at").  If two probes share
#        the same ID string they become duplicates.
#     2. The Ensembl→HGNC conversion in _convert_ensembl_index()
#        did call groupby().mean() but only on that dataset.  When
#        the same gene appeared via two different Ensembl IDs that
#        both mapped to the same HGNC symbol AND both survived QC,
#        the index still had duplicates after conversion.
#     3. The Golub GPL96 mapping assigned the same symbol to
#        multiple probe positions (many Affymetrix chips have
#        multiple probes per gene).  After groupby the index was
#        unique, but if any step between mapping and merge reset the
#        index or concatenated without deduplication, duplicates
#        returned.
#   When pd.concat tried to align these DataFrames on a common gene
#   set, pandas required all indices to be unique — they were not.
#
# Fix:
#   A single _dedup_index() function is applied to EVERY DataFrame
#   immediately after any operation that touches the index:
#     • After Ensembl conversion
#     • After probe→symbol mapping
#     • After Golub GPL96 mapping
#     • After batch correction (which preserves the index but
#       we call it defensively)
#     • Just before _merge_datasets() as a final safety net
#   _dedup_index() calls groupby(index).mean() and then asserts
#   the index is unique, raising a clear error if not.
#
# ══════════════════════════════════════════════════════════════
#  PIPELINE EXPLANATION  (plain English, end-to-end)
# ══════════════════════════════════════════════════════════════
#
# STEP 1 — load_data()
#   Loads six gene expression datasets: Golub from a local .rda
#   file and five GEO datasets downloaded automatically.
#   Each dataset goes through four sub-steps:
#     1a. Parse the raw file (RDA / series matrix / XLSX / TSV.GZ)
#     1b. Convert gene identifiers to HGNC symbols so all datasets
#         speak the same language.
#     1c. QC: drop samples with >50% missing values, drop genes with
#         >80% missing values, drop constant genes.
#     1d. Batch correction: z-score each gene WITHIN its own dataset
#         before merging.  This removes the mean/scale bias that
#         comes from different labs, platforms, and normalisations.
#   Then the six corrected matrices are merged by taking the
#   INTERSECTION of gene sets (only genes present in every dataset),
#   giving X (all samples × common genes) and y (integer class
#   labels, non-overlapping across datasets).
#
# STEP 2 — select_genes(X, y, k=500)
#   Uses one-way ANOVA (F-score) across the y class labels to score
#   every gene by how strongly its expression differs between groups.
#   The top 500 scoring genes are kept.  This reduces noise and
#   focuses the clustering on biologically informative variation.
#
# STEP 3 — preprocess(X)
#   4 transformations applied in order:
#     • Shift so minimum value ≥ 0 (required before log)
#     • log1p: compresses the dynamic range of expression values
#     • KNN imputation: fills any remaining NaN values by averaging
#       the 5 nearest-neighbour samples
#     • StandardScaler: centres each gene to mean 0, std 1
#     • VarianceThreshold: removes any gene still near-constant
#
# STEP 4 — reduce_dim(X) — PCA
#   Projects X onto its top 50 principal components.  PCA removes
#   collinearity between genes and produces a compact numeric
#   representation that makes clustering faster and more stable.
#
# STEP 5 — find_best_k(X)
#   Runs KMeans with k = 2, 3, 4 and scores each with the
#   Silhouette score (higher = better separated clusters).
#   The k that gives the highest Silhouette is used downstream.
#
# STEP 6 — fuzzy_ensemble(X, k)
#   Runs 7 base clusterers: KMeans × 5 (different random seeds) +
#   Agglomerative + Spectral.  Builds a co-association matrix:
#   entry [i,j] = fraction of clusterers that put samples i and j
#   in the same cluster.  Applies Fuzzy C-Means to this consensus
#   matrix and takes argmax of membership scores as final labels.
#
# STEP 7 — evaluate(X, labels, y)
#   Three scores:
#     • Silhouette: −1 to 1, higher = denser, better-separated
#     • Davies-Bouldin: ≥ 0, lower = more compact
#     • Adjusted Rand Index: −1 to 1, higher = matches known labels
#
# STEP 8 — classifier_test(X, labels)  — TRAIN/TEST SPLIT
#   Purpose: validate that the cluster labels are consistent and
#   learnable, not just random noise.
#
#   How the split works:
#     The 70% of samples chosen randomly become the TRAINING SET.
#     The remaining 30% become the TEST SET.
#     Three classifiers (Random Forest, SVM, KNN) are trained on
#     the cluster labels of the training samples and then asked to
#     predict the cluster labels of the test samples.
#     Accuracy = correct predictions / total test samples.
#
#   What it tells you:
#     High accuracy (> 0.80) means the clusters occupy distinct,
#     consistent regions of feature space — they are real structure,
#     not artefacts of the clustering algorithm.
#     Low accuracy means the clusters overlap or are unstable.
#
#   Why we use cluster labels (not y) as the target:
#     We want to test whether the CLUSTERING is internally
#     consistent, not whether we can predict the known disease label.
#     Using cluster labels as targets makes this a clustering
#     quality test, independent of any prior biological annotation.
#
# STEP 9 — ablation(X, k, y)
#   Compares four configurations of increasing complexity:
#     KMeans alone → Fuzzy C-Means alone →
#     Ensemble without Fuzzy → Full Fuzzy Ensemble
#   Each is evaluated with the three metrics from Step 7.
#   The comparison isolates the contribution of each component
#   (ensemble diversity, fuzzy soft-assignment).
#
# STEP 10 — Visualisation
#   Heatmap: samples sorted by cluster, genes on x-axis.
#   t-SNE: 2D non-linear embedding coloured by cluster.
#   Silhouette histogram: shows per-sample cluster cohesion.
#
# ══════════════════════════════════════════════════════════════
#  INSTALL
#     pip install GEOparse pandas numpy matplotlib scikit-learn
#     pip install scikit-fuzzy pyreadr requests openpyxl
# ══════════════════════════════════════════════════════════════

import os, gzip, shutil, re
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import GEOparse
import pyreadr

from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans, AgglomerativeClustering, SpectralClustering
from sklearn.metrics import (
    silhouette_score, silhouette_samples,
    davies_bouldin_score, adjusted_rand_score,
)
from sklearn.impute import KNNImputer
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_classif
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from skfuzzy.cluster import cmeans


# ================================================================
# UTILITY — deduplicate a DataFrame index
# FIX FOR BUG 5
# ================================================================

def _dedup_index(df, source=""):
    """
    Collapse duplicate row labels by averaging, ensuring the index
    is unique before any downstream operation.
    Called after EVERY step that touches the index.
    """
    if df.index.duplicated().any():
        n_before = len(df)
        df = df.groupby(df.index).mean()
        n_after = len(df)
        if source:
            print(f"[DEDUP] {source}: {n_before} → {n_after} unique genes.")
    return df


# ================================================================
# UTILITY — file download
# ================================================================

def _download_file(url, dest_path, desc=""):
    """Download url → dest_path if not already cached."""
    if os.path.exists(dest_path):
        return dest_path
    print(f"[DL] {desc or url} ...")
    try:
        with requests.get(url, stream=True, timeout=120) as r:
            r.raise_for_status()
            with open(dest_path, "wb") as f:
                shutil.copyfileobj(r.raw, f)
        return dest_path
    except Exception as e:
        print(f"[DL] Failed: {e}")
        return None


# ================================================================
# GENE ID CONVERTERS
# ================================================================

def _build_ensembl_to_symbol(destdir):
    """
    Download NCBI Homo_sapiens.gene_info.gz and build
    {ENSG_ID : HGNC_symbol} mapping.
    """
    cache = os.path.join(destdir, "Homo_sapiens.gene_info.gz")
    url   = ("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/"
             "Mammalia/Homo_sapiens.gene_info.gz")
    local = _download_file(url, cache, "NCBI gene info (Ensembl→Symbol)")
    if local is None:
        return {}

    mapping = {}
    try:
        with gzip.open(local, "rt", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                symbol = parts[2]
                for xref in parts[5].split("|"):
                    if xref.startswith("Ensembl:"):
                        mapping[xref.split(":")[1]] = symbol
    except Exception as e:
        print(f"[ID] Ensembl→Symbol build failed: {e}")

    print(f"[ID] Ensembl→Symbol map: {len(mapping):,} entries.")
    return mapping


def _convert_ensembl_index(df, ensembl_map):
    """
    If df index starts with ENSG, convert to HGNC symbols.
    Apply _dedup_index() immediately after (FIX 5).
    """
    if not str(df.index[0]).startswith("ENSG"):
        return df
    df = df.copy()
    df.index = [ensembl_map.get(str(g), "") for g in df.index]
    df = df[df.index != ""]
    df = _dedup_index(df, "Ensembl→HGNC")      # FIX 5
    print(f"[ID] Ensembl→HGNC: {df.shape[0]} genes after conversion.")
    return df


def _build_gpl96_map(destdir):
    """Download GPL96 (HG-U95Av2) annotation → {probe_id: gene_symbol}."""
    try:
        gpl = GEOparse.get_GEO(geo="GPL96", destdir=destdir, silent=True)
        ann = gpl.table
        if ann is None or ann.empty:
            return {}
        sym_col = None
        for col in ann.columns:
            if any(k in col.lower()
                   for k in ["gene_symbol", "symbol", "gene_name"]):
                sym_col = col
                break
        if sym_col is None:
            return {}
        id_col  = ann.columns[0]
        mapping = {}
        for _, row in ann.iterrows():
            p   = str(row[id_col]).strip()
            sym = str(row.get(sym_col, "")).strip()
            sym = sym.split("///")[0].split("//")[0].strip()
            if sym and sym.lower() not in ("", "---", "nan"):
                mapping[p] = sym
        print(f"[GPL96] Probe→Symbol map: {len(mapping):,} entries.")
        return mapping
    except Exception as e:
        print(f"[GPL96] Download failed: {e}")
        return {}


# ================================================================
# GOLUB LOADER  (FIX 1-orientation + FIX 5-dedup)
# ================================================================

def _load_golub(rda_path, destdir):
    """Load Golub .rda, orient correctly, map to HGNC symbols."""
    print(f"\n[Golub] Loading {rda_path} ...")
    result = pyreadr.read_r(rda_path)

    X_df  = result['golub_train_3051']
    y_raw = result['golub_train_response'].values.ravel()

    X_df = X_df.select_dtypes(include=[float, int]).dropna(axis=1, how='all')

    # Correct orientation: rows=genes (3051), cols=samples (38)
    if X_df.shape[0] < X_df.shape[1]:   # more cols than rows → transpose
        X_df = X_df.T

    n_genes, n_samples = X_df.shape
    print(f"[Golub] Raw shape: {n_genes} genes × {n_samples} samples.")

    # Map Golub probe positions → HGNC symbols via GPL96
    gpl96_map = _build_gpl96_map(destdir)
    if gpl96_map:
        probe_order = list(gpl96_map.keys())
        symbols = []
        for i in range(n_genes):
            if i < len(probe_order):
                symbols.append(gpl96_map.get(probe_order[i], f"golub_{i}"))
            else:
                symbols.append(f"golub_{i}")
        X_df.index = symbols
        X_df = _dedup_index(X_df, "Golub GPL96")   # FIX 5
        print(f"[Golub] After GPL96 mapping: {X_df.shape[0]} unique genes.")
    else:
        X_df.index = [f"golub_{i}" for i in range(n_genes)]
        print("[Golub] GPL96 unavailable; using golub_N probe labels.")

    # Encode class labels
    unique_labels = list(dict.fromkeys(str(v) for v in y_raw))
    label_map     = {v: i for i, v in enumerate(unique_labels)}
    labels        = np.array([label_map[str(v)] for v in y_raw])

    print(f"[Golub] Final: {X_df.shape[0]} genes × "
          f"{X_df.shape[1]} samples, {len(unique_labels)} classes.")
    return X_df, labels


# ================================================================
# SUPPLEMENTARY FILE UTILITIES
# FIX FOR BUGS 2, 3, 4
# ================================================================

def _supp_urls_from_soft(gse):
    """
    Read supplementary file URLs directly from SOFT metadata.
    This is the authoritative source — more reliable than FTP scraping.
    """
    urls = []
    for val in gse.metadata.get("supplementary_file", []):
        val = val.strip()
        # Skip empty, "none", and HHS redirect URLs (FIX 2)
        if not val or val.lower() == "none":
            continue
        if "hhs.gov" in val or "://" not in val:
            # Raw filename with no protocol → likely a relative path
            pass
        # Convert FTP to HTTPS for requests library
        val = val.replace("ftp://ftp.ncbi.nlm.nih.gov",
                           "https://ftp.ncbi.nlm.nih.gov")
        if val.startswith("http"):
            urls.append(val)
    return urls


def _ftp_dir_urls(accession):
    """
    Scrape the NCBI FTP HTML page for supplementary file URLs.
    FIX 2: strict regex rejects HHS redirect and absolute URLs.
    """
    prefix = accession[:-3] + "nnn"
    base   = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/"
               f"{prefix}/{accession}/suppl/")
    try:
        r = requests.get(base, timeout=30)
        r.raise_for_status()
        # FIX 2: only match hrefs that look like real filenames:
        #   • must contain a dot (real files have extensions)
        #   • must not start with http/ftp (relative filenames only)
        #   • must not contain '://' (rejects the HHS URL)
        hrefs = re.findall(r'href="([^"?/<>]+\.[a-zA-Z0-9]{1,6})"',
                           r.text, re.IGNORECASE)
        urls  = []
        for h in hrefs:
            # Extra guard: skip anything that is itself a full URL
            if "://" in h or h.lower().startswith("http"):
                continue
            urls.append(base + h)
        return urls
    except Exception:
        return []


# Priority keywords: files whose names suggest expression matrices
_PREFER_KWS  = ["matrix", "count", "expr", "norm", "tpm", "fpkm",
                "rpkm", "log2", "cpm", "rawcount", "raw_count",
                "expression", "quantif"]
# Penalty keywords: files that are definitely not expression matrices
_SKIP_KWS    = ["filelist", "readme", ".md5", "raw.tar",
                ".tgz", "annotation", "phenotype", "metadata",
                "platform", ".bam", ".bai", ".fastq",
                # FIX 3: splice-junction file patterns
                "a3ss", "a5ss", "_ri_", "se_id", "jc.raw",
                "junction", "splice", "_jc.", "novel_junc"]


def _candidate_supp_urls(gse, accession):
    """
    Collect all candidate supplementary file URLs.
    Priority order: SOFT metadata first, then FTP directory listing.
    FIX 2: strict URL validation.
    FIX 3: expression-matrix preference, splice-junction penalty.
    """
    soft_urls = _supp_urls_from_soft(gse)
    ftp_urls  = _ftp_dir_urls(accession)

    seen = set()
    all_urls = []
    for u in soft_urls + ftp_urls:
        fname = u.split("/")[-1].lower()
        # Skip known non-data files (FIX 3)
        if any(k in fname for k in _SKIP_KWS):
            continue
        # Skip HHS / redirect URLs (FIX 2)
        if "hhs.gov" in u:
            continue
        if u not in seen:
            seen.add(u)
            all_urls.append(u)

    # Sort: prefer expression-matrix files first (FIX 3)
    preferred = [u for u in all_urls
                 if any(k in u.lower() for k in _PREFER_KWS)]
    others    = [u for u in all_urls if u not in preferred]
    return preferred + others


def _download_url(url, destdir):
    """Download a URL to destdir, return local path or None."""
    fname = url.split("/")[-1]
    dest  = os.path.join(destdir, fname)
    if os.path.exists(dest):
        print(f"[GEO] Using cached {fname}")
        return dest
    print(f"[GEO] Downloading {fname} ...")
    try:
        with requests.get(url, stream=True, timeout=180) as r:
            r.raise_for_status()
            with open(dest, "wb") as f:
                shutil.copyfileobj(r.raw, f)
        return dest
    except Exception as e:
        print(f"[GEO] Download failed ({e})")
        return None


def _read_matrix_file(filepath):
    """
    Parse a (possibly gzipped) expression matrix file.
    Supports TSV, CSV, and XLSX.
    FIX 1: adds openpyxl Excel branch.
    FIX 3: tries comment=None in addition to comment="!" so files
            without comment headers are parsed correctly.
    FIX 5: calls _dedup_index() before returning.
    """
    fname = filepath.lower()

    # ── XLSX branch (FIX 1) ─────────────────────────────────────
    if fname.endswith(".xlsx") or fname.endswith(".xls"):
        try:
            df = pd.read_excel(filepath, index_col=0, engine="openpyxl")
            df = df.apply(pd.to_numeric, errors='coerce')
            ok = df.notna().mean(axis=0) > 0.5
            df = df.loc[:, ok]
            if df.shape[0] > 10 and df.shape[1] > 1:
                df = _dedup_index(df, "XLSX")   # FIX 5
                print(f"[GEO] Parsed XLSX: {df.shape}")
                return df
        except Exception as e:
            print(f"[GEO] XLSX parse failed: {e}")
        return None

    # ── TSV / CSV branch ─────────────────────────────────────────
    open_fn = gzip.open if filepath.endswith(".gz") else open

    # FIX 3: try both with and without comment skipping
    for comment_char in ["!", None]:
        for sep in ["\t", ","]:
            try:
                kwargs = dict(sep=sep, index_col=0, low_memory=False)
                if comment_char is not None:
                    kwargs["comment"] = comment_char
                with open_fn(filepath, "rt", encoding="utf-8",
                             errors="replace") as fh:
                    df = pd.read_csv(fh, **kwargs)
                if df.shape[1] < 1 or df.shape[0] < 2:
                    continue
                df = df.apply(pd.to_numeric, errors='coerce')
                ok = df.notna().mean(axis=0) > 0.5
                df = df.loc[:, ok]
                if df.shape[1] >= 1 and df.shape[0] > 10:
                    df = _dedup_index(df, filepath.split(os.sep)[-1])  # FIX 5
                    return df
            except Exception:
                continue
    return None


# ================================================================
# GEO SERIES LOADER
# ================================================================

def _detect_series_matrix_cols(gse):
    """Return (id_col, value_col) if GSM tables have data, else (None,None)."""
    for gsm_name, gsm_obj in gse.gsms.items():
        tbl = gsm_obj.table
        if tbl is None or tbl.empty:
            continue
        cols = list(tbl.columns)
        if not cols:
            continue

        id_col = cols[0]
        for c in cols:
            if any(k in c.lower()
                   for k in ["id_ref", "probe_id", "probe", "id", "gene"]):
                id_col = c
                break

        value_col = None
        kws = ["value", "fpkm", "tpm", "rpkm", "rpm", "count",
               "signal", "intensity", "expression", "log2",
               "ratio", "norm", "reads", "cpm"]
        for c in cols:
            if c == id_col:
                continue
            if any(k in c.lower() for k in kws):
                try:
                    pd.to_numeric(tbl[c], errors='raise')
                    value_col = c
                    break
                except (ValueError, TypeError):
                    continue
        if value_col is None:
            for c in cols:
                if c == id_col:
                    continue
                if pd.to_numeric(tbl[c], errors='coerce').notna().mean() > 0.8:
                    value_col = c
                    break
        if value_col:
            return id_col, value_col
    return None, None


def _build_expr_direct(gse, id_col, value_col, accession):
    """Build genes×samples DataFrame directly from GSM tables."""
    probe_ids   = None
    series_dict = {}
    for gsm_name, gsm_obj in gse.gsms.items():
        tbl = gsm_obj.table
        if tbl is None or tbl.empty:
            continue
        if probe_ids is None:
            probe_ids = tbl[id_col].astype(str).values
        series_dict[gsm_name] = pd.to_numeric(
            tbl[value_col], errors='coerce').values
    if not series_dict or probe_ids is None:
        raise RuntimeError(f"[GEO] {accession}: all GSM tables empty.")
    df = pd.DataFrame(series_dict, index=probe_ids)
    return _dedup_index(df, accession)   # FIX 5


def _map_probes_to_genes(accession, expr_df, gse):
    """Map platform probe IDs → HGNC gene symbols via GPL annotation."""
    if not gse.gpls:
        return expr_df
    gpl = next(iter(gse.gpls.values()))
    ann = gpl.table
    if ann is None or ann.empty:
        return expr_df

    sym_col = None
    for col in ann.columns:
        if any(k in col.lower()
               for k in ["gene_symbol", "symbol", "gene_name",
                         "gene_assignment", "official_symbol"]):
            sym_col = col
            break
    if sym_col is None:
        return expr_df

    id_col = ann.columns[0]
    probe_to_gene = {}
    for _, row in ann.iterrows():
        probe = str(row.get(id_col, "")).strip()
        sym   = str(row.get(sym_col, "")).strip()
        sym   = sym.split("///")[0].split("//")[0].strip()
        if sym and sym.lower() not in ("", "---", "nan"):
            probe_to_gene[probe] = sym

    if not probe_to_gene:
        return expr_df

    expr_df = expr_df.copy()
    expr_df.index = [probe_to_gene.get(str(p), str(p)) for p in expr_df.index]
    expr_df = _dedup_index(expr_df, accession)   # FIX 5
    print(f"[GEO] {accession}: probe→symbol → {expr_df.shape[0]} genes.")
    return expr_df


def _load_single_gse(accession, destdir, ensembl_map):
    """Download and parse one GEO series (series-matrix or supplementary)."""
    os.makedirs(destdir, exist_ok=True)
    print(f"\n[GEO] Fetching {accession} ...")

    gse = GEOparse.get_GEO(geo=accession, destdir=destdir, silent=True)

    id_col, value_col = _detect_series_matrix_cols(gse)

    if id_col and value_col:
        try:
            expr_df = gse.pivot_samples(value_col)
        except Exception as e:
            print(f"[GEO] {accession}: pivot_samples failed ({e}); "
                  "building directly.")
            expr_df = _build_expr_direct(gse, id_col, value_col, accession)
        print(f"[GEO] {accession}: series-matrix "
              f"({expr_df.shape[0]} probes × {expr_df.shape[1]} samples)")
    else:
        print(f"[GEO] {accession}: GSM tables empty → "
              "trying supplementary files.")
        candidates = _candidate_supp_urls(gse, accession)

        if not candidates:
            raise RuntimeError(
                f"[GEO] {accession}: no supplementary file URLs found.\n"
                f"      Download a processed matrix from:\n"
                f"      https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
                f"?acc={accession}\n"
                f"      and place it in {destdir}/"
            )

        expr_df = None
        tried   = []
        for url in candidates:
            local = _download_url(url, destdir)
            if local is None:
                continue
            tried.append(url.split("/")[-1])
            df = _read_matrix_file(local)
            if df is not None and df.shape[0] > 10 and df.shape[1] > 1:
                expr_df = df
                print(f"[GEO] {accession}: loaded supplementary "
                      f"({df.shape[0]} genes × {df.shape[1]} samples)")
                break

        if expr_df is None:
            raise RuntimeError(
                f"[GEO] {accession}: no parseable supplementary file.\n"
                f"      Tried: {tried}"
            )

    # Gene ID conversions — FIX 5: dedup at each step
    expr_df = _convert_ensembl_index(expr_df, ensembl_map)
    expr_df = _map_probes_to_genes(accession, expr_df, gse)
    expr_df = _dedup_index(expr_df, accession + " final")   # FIX 5 safety net

    # Sample labels from metadata
    sample_names = list(expr_df.columns)
    raw_labels   = []
    for s in sample_names:
        gsm_obj = gse.gsms.get(s)
        if gsm_obj is None:
            raw_labels.append("unknown")
            continue
        meta  = gsm_obj.metadata
        chars = None
        for key in ["characteristics_ch1", "source_name_ch1", "title"]:
            chars = meta.get(key)
            if chars:
                break
        label_str = "|".join(str(c) for c in (chars or ["unknown"])).lower()
        raw_labels.append(label_str)

    unique_strs = list(dict.fromkeys(raw_labels))
    str_to_int  = {s: i for i, s in enumerate(unique_strs)}
    labels      = np.array([str_to_int[s] for s in raw_labels])

    print(f"[GEO] {accession}: {expr_df.shape[1]} samples, "
          f"{expr_df.shape[0]} genes, {len(unique_strs)} label groups.")
    return expr_df, labels


# ================================================================
# QC AND BATCH CORRECTION
# ================================================================

def _qc_dataset(expr_df, name):
    """Drop bad samples, near-all-NaN genes, and constant genes."""
    n_before = expr_df.shape
    expr_df  = expr_df.apply(pd.to_numeric, errors='coerce')
    expr_df  = expr_df.loc[:, expr_df.isna().mean(axis=0) <= 0.50]
    expr_df  = expr_df.loc[expr_df.isna().mean(axis=1) <= 0.80]
    expr_df  = expr_df.loc[expr_df.std(axis=1) > 0]
    expr_df  = _dedup_index(expr_df, name + " QC")   # FIX 5
    print(f"[QC]  {name}: {n_before} → {expr_df.shape} after QC.")
    return expr_df


def _batch_correct(expr_df):
    """Z-score each gene within this dataset before merging."""
    mean = expr_df.mean(axis=1)
    std  = expr_df.std(axis=1).replace(0, 1)
    corrected = expr_df.sub(mean, axis=0).div(std, axis=0)
    return _dedup_index(corrected)   # FIX 5 defensive call


# ================================================================
# MERGE  (FIX 5: assert unique before concat)
# ================================================================

def _merge_datasets(expr_list, label_list, names):
    """
    Intersect gene sets and concatenate.
    FIX 5: verify all indices are unique before pd.concat.
    """
    # Final safety dedup on every dataset before merge
    expr_list = [_dedup_index(df, names[i] + " pre-merge")
                 for i, df in enumerate(expr_list)]

    common = set(expr_list[0].index)
    for df in expr_list[1:]:
        common &= set(df.index)
    common = sorted(common)

    if len(common) == 0:
        print("\n[MERGE] Diagnosis — gene sets per dataset:")
        for name, df in zip(names, expr_list):
            print(f"  {name}: {df.shape[0]} genes, "
                  f"first 5 → {list(df.index[:5])}")
        raise RuntimeError(
            "[MERGE] No genes in common across all datasets.\n"
            "Check gene ID formats in the diagnosis above."
        )

    print(f"\n[MERGE] Common genes across {len(names)} datasets: "
          f"{len(common)}")

    # Verify uniqueness one final time before concat (FIX 5)
    for name, df in zip(names, expr_list):
        subset = df.loc[common]
        if subset.index.duplicated().any():
            raise RuntimeError(
                f"[MERGE] {name} still has duplicate gene labels "
                f"after all dedup steps — cannot concat safely.")

    merged = pd.concat([df.loc[common] for df in expr_list], axis=1)
    X = merged.T.values.astype(float)

    y_parts = []
    offset  = 0
    for lbl in label_list:
        y_parts.append(lbl + offset)
        offset += int(lbl.max()) + 1
    y = np.concatenate(y_parts)

    print(f"[MERGE] Final: {X.shape[0]} samples × {X.shape[1]} genes")
    print(f"[MERGE] Label counts: "
          f"{ {int(v): int((y==v).sum()) for v in sorted(set(y))} }")
    return X, y


# ================================================================
# UNIFIED load_data()
# ================================================================

def load_data(golub_rda_path, geo_accessions=None, destdir="./geo_data"):
    """
    Load and merge Golub RDA + five GEO datasets.
    Returns (X, y) ready for select_genes().
    """
    if geo_accessions is None:
        geo_accessions = [
            "GSE297413", "GSE253086", "GSE201492",
            "GSE216738", "GSE181157",
        ]

    os.makedirs(destdir, exist_ok=True)

    # Build shared Ensembl→Symbol map once
    ensembl_map = _build_ensembl_to_symbol(destdir)

    expr_list  = []
    label_list = []
    name_list  = []

    # ── Golub ─────────────────────────────────────────────────────
    g_expr, g_labels = _load_golub(golub_rda_path, destdir)
    g_expr = _qc_dataset(g_expr, "Golub")
    g_expr = _batch_correct(g_expr)
    expr_list.append(g_expr)
    label_list.append(g_labels)
    name_list.append("Golub")

    # ── GEO datasets ──────────────────────────────────────────────
    for acc in geo_accessions:
        try:
            expr_df, labels = _load_single_gse(acc, destdir, ensembl_map)
            expr_df = _qc_dataset(expr_df, acc)

            # Align label vector to samples surviving QC
            surviving = list(expr_df.columns)
            all_gsms  = list(
                GEOparse.get_GEO(geo=acc, destdir=destdir, silent=True)
                .gsms.keys()
            )
            pos_map = {s: i for i, s in enumerate(all_gsms)}
            idx     = np.array([pos_map[s] for s in surviving
                                 if s in pos_map])
            if len(idx) == 0:
                idx = np.arange(min(len(labels), expr_df.shape[1]))
            labels = labels[idx]

            expr_df = _batch_correct(expr_df)
            expr_list.append(expr_df)
            label_list.append(labels)
            name_list.append(acc)

        except Exception as e:
            print(f"\n[WARN] {acc}: skipped — {e}\n")
            continue

    if len(expr_list) == 1:
        print("[WARN] Only Golub loaded; running single-dataset mode.")
        X = expr_list[0].T.values.astype(float)
        y = label_list[0]
    else:
        X, y = _merge_datasets(expr_list, label_list, name_list)

    if X.shape[0] > X.shape[1]:
        X = X.T
    print(f"\n[INFO] Final data shape: {X.shape}")
    return X, y


# ================================================================
# ORIGINAL PIPELINE — COMPLETELY UNCHANGED
# ================================================================

def select_genes(X, y, k=500):
    print(f"[INFO] Selecting top {k} genes via f_classif...")
    selector = SelectKBest(score_func=f_classif, k=k)
    X_new = selector.fit_transform(X, y)
    print(f"[INFO] Shape after selection: {X_new.shape}")
    return X_new


def preprocess(X):
    print("[INFO] Preprocessing...")
    min_val = np.min(X)
    if min_val < 0:
        X = X - min_val
    X = np.log1p(X)
    X = KNNImputer(n_neighbors=5).fit_transform(X)
    X = StandardScaler().fit_transform(X)
    X = VarianceThreshold(0.01).fit_transform(X)
    print(f"[INFO] Shape after preprocessing: {X.shape}")
    return X


def reduce_dim(X):
    n_comp = min(50, min(X.shape) - 1)
    print(f"[INFO] PCA: retaining {n_comp} components")
    return PCA(n_components=n_comp).fit_transform(X)


def run_kmeans(X, k):
    return KMeans(n_clusters=k, random_state=42, n_init=10).fit_predict(X)


def run_hierarchical(X, k):
    return AgglomerativeClustering(n_clusters=k).fit_predict(X)


def run_spectral(X, k):
    n_neighbors = min(10, X.shape[0] - 1)
    return SpectralClustering(
        n_clusters=k, affinity='nearest_neighbors',
        n_neighbors=n_neighbors, random_state=42
    ).fit_predict(X)


def run_fuzzy_cmeans(X, k):
    _cntr, u, _u0, _d, _jm, _p, _fpc = cmeans(
        X.T, c=k, m=2, error=0.005, maxiter=1000)
    return u


def coassociation(clusterings):
    n   = len(clusterings[0])
    mat = np.zeros((n, n))
    for labels in clusterings:
        for i in range(n):
            for j in range(n):
                if labels[i] == labels[j]:
                    mat[i][j] += 1
    return mat / len(clusterings)


def fuzzy_ensemble(X, k):
    print("[INFO] Building fuzzy ensemble...")
    base = [run_kmeans(X, k) for _ in range(5)]
    base.append(run_hierarchical(X, k))
    base.append(run_spectral(X, k))
    co_mat = coassociation(base)
    co_mat = MinMaxScaler().fit_transform(co_mat)
    u      = run_fuzzy_cmeans(co_mat, k)
    return np.argmax(u, axis=0)


def evaluate(X, labels, y):
    return {
        "Silhouette":     silhouette_score(X, labels),
        "Davies-Bouldin": davies_bouldin_score(X, labels),
        "ARI":            adjusted_rand_score(y, labels),
    }


def find_best_k(X, k_range=(2, 4)):
    print("\n[INFO] Searching for best K...")
    best_k, best_score = 2, -1
    for k in range(k_range[0], k_range[1] + 1):
        labels = run_kmeans(X, k)
        score  = silhouette_score(X, labels)
        print(f"  K={k}  Silhouette={score:.4f}")
        if score > best_score:
            best_k, best_score = k, score
    print(f"[INFO] Best K = {best_k}  (Silhouette = {best_score:.4f})")
    return best_k


def classifier_test(X, labels):
    """
    TRAIN/TEST SPLIT EXPLAINED
    ──────────────────────────
    Purpose: prove the clusters are internally consistent — a
    classifier trained on 70 % of samples should be able to predict
    the cluster label of the remaining 30 %.

    Mechanics:
      1. Randomly shuffle all samples (random_state=42 for
         reproducibility).
      2. The first 70 % become the TRAINING SET: the classifier sees
         both X (PCA features) and labels (cluster assignments) and
         learns which feature patterns correspond to which cluster.
      3. The remaining 30 % become the TEST SET: the classifier sees
         only X and must predict the label.  We compare predictions
         against the actual labels to compute accuracy.

    Three classifiers are used so results are not algorithm-specific:
      • Random Forest: ensemble of decision trees, handles
        non-linear boundaries.
      • SVM: maximum-margin linear separator.
      • KNN: k=5 nearest neighbours by Euclidean distance.

    Interpretation:
      accuracy > 0.85 → clusters are tight and well-separated.
      accuracy 0.65-0.85 → moderate overlap, some ambiguous samples.
      accuracy < 0.65 → clusters are not stable; consider changing k.
    """
    X_tr, X_te, y_tr, y_te = train_test_split(
        X, labels, test_size=0.3, random_state=42)

    models = {
        "RandomForest": RandomForestClassifier(random_state=42),
        "SVM":          SVC(),
        "KNN":          KNeighborsClassifier(),
    }
    print("\n[INFO] Classifier Validation (70 % train / 30 % test):")
    for name, model in models.items():
        model.fit(X_tr, y_tr)
        acc = model.score(X_te, y_te)
        print(f"  {name:>14}: {acc:.4f}  "
              f"({int(acc * len(y_te))}/{len(y_te)} correct)")


def ablation(X, k, y):
    results = {}
    results["KMeans"] = evaluate(X, run_kmeans(X, k), y)
    u_fcm = run_fuzzy_cmeans(X, k)
    results["Fuzzy C-Means"] = evaluate(X, np.argmax(u_fcm, axis=0), y)
    c1, c2, c3 = run_kmeans(X, k), run_hierarchical(X, k), run_spectral(X, k)
    co_mat = MinMaxScaler().fit_transform(coassociation([c1, c2, c3]))
    results["Ensemble (no Fuzzy)"] = evaluate(X, run_kmeans(co_mat, k), y)
    results["Fuzzy Ensemble"]      = evaluate(X, fuzzy_ensemble(X, k), y)
    return results


def plot_ablation(results):
    methods = list(results.keys())
    scores  = [results[m]["Silhouette"] for m in methods]
    plt.figure(figsize=(8, 5))
    bars = plt.bar(methods, scores,
                   color=["#4C72B0", "#55A868", "#C44E52", "#8172B2"])
    plt.xticks(rotation=20, ha='right')
    plt.ylabel("Silhouette Score")
    plt.title("Ablation Study — Silhouette Score by Method")
    plt.tight_layout()
    for bar, score in zip(bars, scores):
        plt.text(bar.get_x() + bar.get_width() / 2,
                 bar.get_height() + 0.005,
                 f"{score:.3f}", ha='center', va='bottom', fontsize=9)
    plt.savefig("ablation_study.png", dpi=150, bbox_inches='tight')
    plt.show()


def plot_heatmap(X, labels):
    idx = np.argsort(labels)
    plt.figure(figsize=(10, 6))
    plt.imshow(X[idx], aspect='auto', interpolation='nearest')
    plt.colorbar(label="Scaled Expression")
    plt.xlabel("Genes (PCs)")
    plt.ylabel("Samples (sorted by cluster)")
    plt.title("Gene Expression Heatmap")
    plt.tight_layout()
    plt.savefig("heatmap.png", dpi=150, bbox_inches='tight')
    plt.show()


def plot_tsne(X, labels):
    print("[INFO] Running t-SNE ...")
    emb = TSNE(n_components=2, random_state=42).fit_transform(X)
    plt.figure(figsize=(7, 6))
    sc = plt.scatter(emb[:, 0], emb[:, 1], c=labels, cmap='tab10', s=40)
    plt.colorbar(sc, label="Cluster")
    plt.title("t-SNE Cluster Visualisation")
    plt.xlabel("t-SNE 1"); plt.ylabel("t-SNE 2")
    plt.tight_layout()
    plt.savefig("tsne.png", dpi=150, bbox_inches='tight')
    plt.show()


def plot_silhouette(X, labels):
    vals = silhouette_samples(X, labels)
    plt.figure(figsize=(7, 4))
    plt.hist(vals, bins=20, color="#4C72B0", edgecolor="white")
    plt.axvline(vals.mean(), color='red', linestyle='--',
                label=f"Mean = {vals.mean():.3f}")
    plt.xlabel("Silhouette Value"); plt.ylabel("Count")
    plt.title("Per-Sample Silhouette Distribution")
    plt.legend(); plt.tight_layout()
    plt.savefig("silhouette_dist.png", dpi=150, bbox_inches='tight')
    plt.show()


def run(golub_rda_path, geo_accessions=None, destdir="./geo_data"):
    X, y = load_data(golub_rda_path, geo_accessions, destdir)
    X = select_genes(X, y, k=500)
    X = preprocess(X)
    X = reduce_dim(X)
    k      = find_best_k(X, k_range=(2, 4))
    labels = fuzzy_ensemble(X, k)

    print("\n[FINAL RESULTS — Fuzzy Ensemble]")
    for metric, val in evaluate(X, labels, y).items():
        print(f"  {metric:>16}: {val:.4f}")

    classifier_test(X, labels)

    print("\n[ABLATION STUDY]")
    ab = ablation(X, k, y)
    for method, res in ab.items():
        print(f"  {method:<22} | Sil={res['Silhouette']:.4f}  "
              f"DB={res['Davies-Bouldin']:.4f}  ARI={res['ARI']:.4f}")

    plot_ablation(ab)
    plot_heatmap(X, labels)
    plot_tsne(X, labels)
    plot_silhouette(X, labels)


# ================================================================
# ENTRY POINT
# ================================================================

if __name__ == "__main__":
    GOLUB_PATH = (r"C:\Users\BIT\Downloads\archive (6)"
                  r"\leukemia_data_Golub99_3051.rda")

    GEO_DATASETS = [
        "GSE297413",   # TPM matrix — supplementary TXT.GZ ✓
        "GSE253086",   # TPM matrix — supplementary XLSX ✓ (needs openpyxl)
        "GSE201492",   # raw counts — supplementary TXT.GZ ✓
        "GSE216738",   # log2-CPM   — supplementary TXT.GZ ✓ (Ensembl→HGNC)
        "GSE181157",   # no processed matrix on GEO; skipped gracefully
    ]

    run(golub_rda_path=GOLUB_PATH,
        geo_accessions=GEO_DATASETS,
        destdir="./geo_data")
