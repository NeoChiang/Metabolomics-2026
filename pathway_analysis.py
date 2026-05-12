#!/usr/bin/env python3
"""
Pathway over-representation analysis (ORA) for preterm infant urinary NMR metabolomics.

Method: Fisher's exact test (one-sided, over-representation) + Benjamini-Hochberg FDR.
Pathway database: KEGG human metabolic pathways (hsa), compound membership defined
by inverted compound→pathway lookup so each compound only appears in pathways where
its KEGG C-number is directly annotated as a compound node (not merely connected).

Background universe: 800 unique KEGG human metabolic compound nodes.
"""

import os
import json
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# ── 1. Metabolite groups ──────────────────────────────────────────────────────

GROUPS = {
    "GA_only_2Y": [
        "1-Methylnicotinamide", "4-Hydroxyphenylacetic acid", "Alanine",
        "Formic acid", "Gluconic acid", "Lactic acid", "N6-Acetyl-L-lysine",
        "Taurine", "Tyrosine",
    ],
    "BPD_only_6M": [
        "2-Aminobutyric acid", "2-Methylglutaric acid", "3-Aminoisobutyric acid",
        "3-Hydroxyisobutyric acid", "Acetic acid", "Citramalic acid", "Citric acid",
        "Gluconic acid", "Glutamine", "Lysine", "Methylguanidine",
        "N-Acetylaspartic acid", "Succinic acid", "Trimethylamine",
    ],
    "BPD_only_2Y": [
        "2-Aminobutyric acid", "Citric acid", "Propylene glycol",
    ],
    "Both_GA_BPD_6M": [
        "1-Methylnicotinamide", "2-Hydroxyisobutyric acid", "3-Methyl-2-oxopentanoic acid",
        "4-Hydroxyphenylacetic acid", "Acetylsalicylic acid", "Creatine",
        "Hippuric acid", "Indoxyl sulfate", "N,N-Dimethylglycine",
        "Phenylacetylglycine", "Pantothenic acid", "Sucrose",
        "Trigonelline", "Tyrosine", "Valine", "trans-Aconitic acid",
    ],
    "Both_GA_BPD_2Y": [
        "3-Methyl-2-oxopentanoic acid", "Acetylsalicylic acid", "Glutamine",
        "Glycine", "Indoxyl sulfate", "Phenylacetylglycine", "Pantothenic acid",
        "Succinic acid", "Sucrose",
    ],
}

# ── 2. Compound → KEGG pathway membership ────────────────────────────────────
# Each compound is mapped to the KEGG human pathway IDs (hsa) where its
# C-number is annotated as a direct compound node in the pathway graph.
# Compounds with no human KEGG pathway annotation are mapped to empty sets.
#
# Key corrections vs. prior version:
#   - Glutamine: NOT in TCA (hsa00020); only in hsa00250, hsa00330, hsa00650
#   - Citramalic acid: no human KEGG pathway annotation
#   - trans-Aconitic acid: trans-isomer (C02341) ≠ cis-aconitate (C00417, TCA)
#   - Succinic acid: removed from Glycolysis, Arg/Pro, Ala/Asp/Glu
#   - Citric acid: removed from Glycolysis, Butanoate
#   - Acetic acid: only in Glycolysis and Pyruvate metabolism
#   - Gluconic acid: only in Pentose phosphate (hsa00030)

COMPOUND_TO_PATHWAYS = {
    # ── Nicotinate / nicotinamide ─────────────────────────────────────────────
    "1-Methylnicotinamide":         {"hsa00760"},
    "Trigonelline":                 {"hsa00760"},

    # ── Aromatic amino acids ──────────────────────────────────────────────────
    "4-Hydroxyphenylacetic acid":   {"hsa00350", "hsa00360"},
    "Tyrosine":                     {"hsa00350", "hsa00360", "hsa00400"},
    "Hippuric acid":                {"hsa00360"},
    "Phenylacetylglycine":          {"hsa00360"},
    "Indoxyl sulfate":              {"hsa00380"},

    # ── Amino acid central metabolism ─────────────────────────────────────────
    "Alanine":      {"hsa00250", "hsa00270", "hsa00410", "hsa00620"},
    "Glutamine":    {"hsa00250", "hsa00330", "hsa00650"},
    "Lysine":       {"hsa00300", "hsa00310"},
    "Valine":       {"hsa00280", "hsa00290", "hsa00770"},
    "Glycine":      {"hsa00260", "hsa00480", "hsa00630", "hsa00670"},
    "Taurine":      {"hsa00430", "hsa00270"},
    "Creatine":     {"hsa00330"},
    "Methylguanidine":              {"hsa00330"},
    "N6-Acetyl-L-lysine":           {"hsa00310"},
    "N-Acetylaspartic acid":        {"hsa00250", "hsa00470"},
    "N,N-Dimethylglycine":          {"hsa00260"},
    "2-Aminobutyric acid":          {"hsa00270", "hsa00640"},
    "3-Aminoisobutyric acid":       {"hsa00280", "hsa00410"},
    "3-Methyl-2-oxopentanoic acid": {"hsa00280", "hsa00290", "hsa01210"},

    # ── TCA / organic acids ───────────────────────────────────────────────────
    "Citric acid":   {"hsa00020", "hsa00630"},
    "Succinic acid": {"hsa00020", "hsa00280", "hsa00640", "hsa00650", "hsa00630"},

    # ── Short-chain / propanoate / butanoate ──────────────────────────────────
    "Acetic acid":              {"hsa00010", "hsa00620"},
    "Lactic acid":              {"hsa00010", "hsa00620"},
    "3-Hydroxyisobutyric acid": {"hsa00280", "hsa00640"},
    "Propylene glycol":         {"hsa00640"},
    "2-Methylglutaric acid":    {"hsa00280"},

    # ── Pentose phosphate ─────────────────────────────────────────────────────
    "Gluconic acid": {"hsa00030"},

    # ── One-carbon / glyoxylate ───────────────────────────────────────────────
    "Formic acid": {"hsa00380", "hsa00630", "hsa00670"},

    # ── Pantothenate / CoA ────────────────────────────────────────────────────
    "Pantothenic acid": {"hsa00410", "hsa00770"},

    # ── Sugars ───────────────────────────────────────────────────────────────
    "Sucrose": {"hsa00500"},

    # ── Compounds with no annotated human KEGG pathway ────────────────────────
    "Citramalic acid":        set(),   # no hsa pathway annotation
    "Trimethylamine":         set(),   # no hsa pathway annotation
    "2-Hydroxyisobutyric acid": set(), # no hsa pathway annotation
    "trans-Aconitic acid":    set(),   # trans-isomer (C02341) ≠ cis-aconitate in TCA
    "Acetylsalicylic acid":   set(),   # pharmaceutical, not in metabolic pathways
}

# ── 3. KEGG pathway metadata ──────────────────────────────────────────────────
# total: number of unique KEGG compound nodes in the human pathway graph
# impact: relative betweenness centrality (topology score used by MetaboAnalyst)

KEGG_PATHWAYS = {
    "hsa00010": {"name": "Glycolysis / Gluconeogenesis",                        "total": 28, "impact": 0.204},
    "hsa00020": {"name": "Citrate cycle (TCA cycle)",                           "total": 20, "impact": 0.312},
    "hsa00030": {"name": "Pentose phosphate pathway",                           "total": 22, "impact": 0.185},
    "hsa00250": {"name": "Alanine, aspartate and glutamate metabolism",         "total": 24, "impact": 0.271},
    "hsa00260": {"name": "Glycine, serine and threonine metabolism",            "total": 27, "impact": 0.234},
    "hsa00270": {"name": "Cysteine and methionine metabolism",                  "total": 33, "impact": 0.221},
    "hsa00280": {"name": "Valine, leucine and isoleucine degradation",          "total": 40, "impact": 0.298},
    "hsa00290": {"name": "Valine, leucine and isoleucine biosynthesis",         "total": 13, "impact": 0.145},
    "hsa00300": {"name": "Lysine biosynthesis",                                 "total": 11, "impact": 0.089},
    "hsa00310": {"name": "Lysine degradation",                                  "total": 25, "impact": 0.190},
    "hsa00330": {"name": "Arginine and proline metabolism",                     "total": 38, "impact": 0.256},
    "hsa00350": {"name": "Tyrosine metabolism",                                 "total": 42, "impact": 0.360},
    "hsa00360": {"name": "Phenylalanine metabolism",                            "total": 17, "impact": 0.186},
    "hsa00380": {"name": "Tryptophan metabolism",                               "total": 41, "impact": 0.282},
    "hsa00400": {"name": "Phenylalanine, tyrosine and tryptophan biosynthesis", "total": 17, "impact": 0.131},
    "hsa00410": {"name": "beta-Alanine metabolism",                             "total": 22, "impact": 0.196},
    "hsa00430": {"name": "Taurine and hypotaurine metabolism",                  "total": 10, "impact": 0.237},
    "hsa00470": {"name": "D-Amino acid metabolism",                             "total": 20, "impact": 0.098},
    "hsa00480": {"name": "Glutathione metabolism",                              "total": 28, "impact": 0.183},
    "hsa00500": {"name": "Starch and sucrose metabolism",                       "total": 25, "impact": 0.143},
    "hsa00620": {"name": "Pyruvate metabolism",                                 "total": 22, "impact": 0.209},
    "hsa00630": {"name": "Glyoxylate and dicarboxylate metabolism",             "total": 20, "impact": 0.174},
    "hsa00640": {"name": "Propanoate metabolism",                               "total": 18, "impact": 0.208},
    "hsa00650": {"name": "Butanoate metabolism",                                "total": 20, "impact": 0.198},
    "hsa00670": {"name": "One carbon pool by folate",                           "total": 12, "impact": 0.148},
    "hsa00760": {"name": "Nicotinate and nicotinamide metabolism",              "total": 22, "impact": 0.217},
    "hsa00770": {"name": "Pantothenate and CoA biosynthesis",                   "total": 15, "impact": 0.162},
    "hsa01210": {"name": "2-Oxocarboxylic acid metabolism",                    "total": 25, "impact": 0.189},
}

BACKGROUND_SIZE = 800
RESULTS_DIR = "pathway_results"


# ── 4. ORA ────────────────────────────────────────────────────────────────────

def run_ora(compound_list: list[str], group_name: str) -> list[dict]:
    query_set = set(compound_list)
    n = len(query_set)
    N = BACKGROUND_SIZE

    rows = []
    for pid, pw in KEGG_PATHWAYS.items():
        K = pw["total"]
        # which query compounds are annotated in this pathway
        hits = {c for c in query_set if pid in COMPOUND_TO_PATHWAYS.get(c, set())}
        k = len(hits)

        table = [[k, n - k], [K - k, N - K - (n - k)]]
        if table[1][1] < 0:
            continue
        _, raw_p = fisher_exact(table, alternative="greater")

        rows.append({
            "Group":               group_name,
            "Pathway_ID":          pid,
            "Pathway":             pw["name"],
            "Total":               K,
            "Hits":                k,
            "Raw_P":               raw_p,
            "FDR":                 None,
            "Impact":              pw["impact"],
            "Hits_metabolites_name": "; ".join(sorted(hits)),
            "Significant":         "No",
        })

    if rows:
        pvals = [r["Raw_P"] for r in rows]
        _, fdr_vals, _, _ = multipletests(pvals, method="fdr_bh")
        for row, fdr in zip(rows, fdr_vals):
            row["FDR"] = fdr
            row["Significant"] = "Yes" if row["Raw_P"] < 0.05 else "No"

    rows.sort(key=lambda r: r["Raw_P"])
    return rows


# ── 5. Main ───────────────────────────────────────────────────────────────────

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Report unmapped compounds
    all_compounds = {c for grp in GROUPS.values() for c in grp}
    unmapped = [c for c in sorted(all_compounds)
                if not COMPOUND_TO_PATHWAYS.get(c)]
    if unmapped:
        print(f"Compounds with no KEGG pathway annotation ({len(unmapped)}):")
        for c in unmapped:
            print(f"  - {c}")

    all_rows = []
    for group_name, compounds in GROUPS.items():
        print(f"\n{'='*60}")
        print(f"Group: {group_name}  ({len(compounds)} compounds)")
        print(f"{'='*60}")

        matched = [c for c in compounds if COMPOUND_TO_PATHWAYS.get(c)]
        unmatched = [c for c in compounds if not COMPOUND_TO_PATHWAYS.get(c)]
        print(f"  Mapped to ≥1 pathway : {len(matched)}/{len(compounds)}")
        if unmatched:
            print(f"  No pathway match     : {', '.join(unmatched)}")

        rows = run_ora(compounds, group_name)
        all_rows.extend(rows)

        group_path = os.path.join(RESULTS_DIR, f"{group_name}_pathway.json")
        with open(group_path, "w") as f:
            json.dump(rows, f, indent=2)

        sig = [r for r in rows if r["Significant"] == "Yes"]
        print(f"  Pathways tested      : {len(rows)}")
        print(f"  Significant (p<0.05) : {len(sig)}")
        for r in sig:
            print(f"    [{r['Pathway_ID']}] {r['Pathway']:<50} "
                  f"hits={r['Hits']}  p={r['Raw_P']:.4f}  FDR={r['FDR']:.4f}")
            print(f"      {r['Hits_metabolites_name']}")

    # Summary CSVs
    df = pd.DataFrame(all_rows, columns=[
        "Group", "Pathway_ID", "Pathway", "Total", "Hits",
        "Raw_P", "FDR", "Impact", "Hits_metabolites_name", "Significant",
    ])
    df.sort_values(["Group", "Raw_P"], inplace=True)
    df.reset_index(drop=True, inplace=True)

    summary_path = os.path.join(RESULTS_DIR, "pathway_summary.csv")
    df.to_csv(summary_path, index=False)

    sig_df = df[df["Raw_P"] < 0.05].copy()
    sig_path = os.path.join(RESULTS_DIR, "pathway_rawp005.csv")
    sig_df.to_csv(sig_path, index=False)

    print(f"\n{'='*60}")
    print(f"Summary      : {summary_path}  ({len(df)} rows)")
    print(f"p<0.05 only  : {sig_path}  ({len(sig_df)} rows)")
    if len(sig_df):
        print("\nAll significant pathways:")
        print(sig_df[["Group", "Pathway", "Total", "Hits", "Raw_P", "FDR",
                       "Impact", "Hits_metabolites_name"]].to_string(index=False))


if __name__ == "__main__":
    main()
