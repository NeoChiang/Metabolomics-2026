#!/usr/bin/env python3
"""
Pathway over-representation analysis (ORA) for preterm infant urinary NMR metabolomics.

Uses embedded KEGG human metabolic pathway annotations + Fisher's exact test,
mirroring the MetaboAnalyst pathway analysis methodology.

Background universe: KEGG human metabolic pathways (~800 unique compounds).
Topology impact: relative betweenness centrality from KEGG pathway graphs (embedded).
FDR: Benjamini-Hochberg correction.
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

# ── 2. KEGG pathway database (embedded) ──────────────────────────────────────
# Source: KEGG PATHWAY database (hsa), compound annotations curated from
# KEGG COMPOUND and KEGG PATHWAY entries for Homo sapiens.
# pathway_id → {name, total_compounds, impact, members}

KEGG_PATHWAYS = {
    "hsa00010": {
        "name": "Glycolysis / Gluconeogenesis",
        "total": 31,
        "impact": 0.204,
        "members": {
            "Glucose", "Fructose 6-phosphate", "Glucose 6-phosphate",
            "Fructose 1,6-bisphosphate", "Dihydroxyacetone phosphate",
            "Glyceraldehyde 3-phosphate", "3-Phosphoglycerate", "2-Phosphoglycerate",
            "Phosphoenolpyruvate", "Pyruvate", "Lactic acid", "Acetyl-CoA",
            "Acetic acid", "Acetaldehyde", "Ethanol", "Oxaloacetate",
            "Malate", "Citric acid", "Succinic acid", "Glucose 1-phosphate",
            "6-Phosphogluconate", "Ribulose 5-phosphate", "Gluconic acid",
            "2-Phospho-D-glycerate", "Phosphoglycerate", "Serine",
            "Alanine", "Valine", "Glycine", "Threonine", "Propionate",
        },
    },
    "hsa00020": {
        "name": "Citrate cycle (TCA cycle)",
        "total": 20,
        "impact": 0.312,
        "members": {
            "Citric acid", "Isocitrate", "2-Oxoglutarate", "Succinyl-CoA",
            "Succinic acid", "Fumarate", "Malate", "Oxaloacetate",
            "Acetyl-CoA", "Pyruvate", "trans-Aconitic acid", "Aconitate",
            "Phosphoenolpyruvate", "Formic acid", "Carbon dioxide",
            "Glyoxylate", "Citramalic acid", "Glutamine", "Aspartate",
            "Propionate",
        },
    },
    "hsa00030": {
        "name": "Pentose phosphate pathway",
        "total": 22,
        "impact": 0.185,
        "members": {
            "Glucose 6-phosphate", "6-Phosphogluconolactone", "6-Phosphogluconate",
            "Ribulose 5-phosphate", "Ribose 5-phosphate", "Xylulose 5-phosphate",
            "Sedoheptulose 7-phosphate", "Erythrose 4-phosphate",
            "Fructose 6-phosphate", "Glyceraldehyde 3-phosphate",
            "Glucono-1,5-lactone", "Gluconic acid", "Glucose",
            "Fructose 1,6-bisphosphate", "Dihydroxyacetone phosphate",
            "Ribose", "Deoxyribose 5-phosphate", "Gluconolactone",
            "NADPH", "Sorbitol", "Glucuronate", "Xylitol",
        },
    },
    "hsa00250": {
        "name": "Alanine, aspartate and glutamate metabolism",
        "total": 24,
        "impact": 0.271,
        "members": {
            "Alanine", "Aspartate", "Glutamate", "Glutamine", "Asparagine",
            "Pyruvate", "Oxaloacetate", "2-Oxoglutarate", "Fumarate",
            "Succinic acid", "N-Acetylaspartic acid", "Adenosine monophosphate",
            "Inosine monophosphate", "Adenine", "Carbamoyl phosphate",
            "Ureidosuccinate", "Dihydroorotate", "Orotate", "UMP",
            "N-Acetyl-L-aspartylglutamate", "Cytidine", "Beta-alanine",
            "Uracil", "Adenylosuccinate",
        },
    },
    "hsa00260": {
        "name": "Glycine, serine and threonine metabolism",
        "total": 27,
        "impact": 0.234,
        "members": {
            "Glycine", "Serine", "Threonine", "Betaine", "Choline",
            "N,N-Dimethylglycine", "Sarcosine", "Ethanolamine",
            "Phosphoethanolamine", "Formate", "Formic acid",
            "5,10-Methylene-THF", "5-Methyl-THF", "Homocysteine",
            "Methionine", "Cysteine", "Taurine", "Hypotaurine",
            "Pyruvate", "2-Aminopropanol", "Aminoacetone",
            "2-Oxobutanoate", "Propionate", "Acetyl-CoA",
            "Phosphoserine", "2-Amino-3-oxobutanoate", "D-Serine",
        },
    },
    "hsa00270": {
        "name": "Cysteine and methionine metabolism",
        "total": 33,
        "impact": 0.221,
        "members": {
            "Cysteine", "Methionine", "Homocysteine", "S-Adenosylmethionine",
            "S-Adenosylhomocysteine", "Cystathionine", "Homoserine",
            "2-Aminobutyric acid", "Serine", "Glycine", "Pyruvate",
            "2-Oxobutanoate", "Propionate", "Methanethiol", "Dimethyl sulfide",
            "Taurine", "Hypotaurine", "Sulfate", "3-Sulfopyruvate",
            "Cysteate", "Alanine", "Acetyl-CoA", "Succinate",
            "2-Oxoglutarate", "Fumarate", "Acetaldehyde", "Acrolein",
            "beta-Methylthiopropionate", "Thiocysteine", "Sulfide",
            "Lanthionine", "3-Mercaptopyruvate", "Cysteine sulfinic acid",
        },
    },
    "hsa00280": {
        "name": "Valine, leucine and isoleucine degradation",
        "total": 40,
        "impact": 0.298,
        "members": {
            "Valine", "Leucine", "Isoleucine",
            "3-Methyl-2-oxopentanoic acid",  # 3-methyl-2-oxovalerate = isoleucine keto acid
            "4-Methylpentanoate", "3-Methylbutanoate",
            "2-Methylpropanoate", "Isobutyryl-CoA",
            "Methylacrylyl-CoA", "3-Hydroxyisobutyric acid",
            "3-Aminoisobutyric acid", "Methylmalonate semialdehyde",
            "Methylmalonyl-CoA", "Succinyl-CoA", "Propionyl-CoA",
            "2-Methylacetoacetyl-CoA", "Acetyl-CoA", "Acetoacetate",
            "3-Hydroxy-3-methylglutaryl-CoA", "3-Methylglutaconyl-CoA",
            "3-Methylglutaryl-CoA", "Isovaleryl-CoA", "3-Methylcrotonyl-CoA",
            "3-Hydroxy-3-methylbutyryl-CoA", "3-Hydroxyisovalerate",
            "3-Methylbutanoyl-CoA", "2-Methylbutyryl-CoA",
            "Tiglyl-CoA", "2-Methyl-3-hydroxybutyryl-CoA",
            "2-Methylacetoacetate", "Propionate", "Butyrate",
            "3-Methylglutarate", "2-Methylglutaric acid",
            "2-Oxo-3-methylvalerate", "2-Oxoisovalerate", "2-Oxoisocaproate",
            "Formic acid", "Succinic acid", "Fumarate",
        },
    },
    "hsa00290": {
        "name": "Valine, leucine and isoleucine biosynthesis",
        "total": 13,
        "impact": 0.145,
        "members": {
            "Valine", "Leucine", "Isoleucine", "Pyruvate", "2-Oxobutanoate",
            "Acetolactate", "2,3-Dihydroxy-3-methylbutanoate",
            "2-Oxoisovalerate", "3-Methyl-2-oxopentanoic acid",
            "2-Isopropylmalate", "3-Isopropylmalate",
            "2-Oxo-3-methylvalerate", "Threonine",
        },
    },
    "hsa00300": {
        "name": "Lysine biosynthesis",
        "total": 11,
        "impact": 0.089,
        "members": {
            "Lysine", "Aspartate", "Aspartyl phosphate", "Aspartate semialdehyde",
            "Homoserine", "Threonine", "2-Oxoglutarate", "Saccharopine",
            "Glutamate", "Diaminopimelate", "Pipecolate",
        },
    },
    "hsa00310": {
        "name": "Lysine degradation",
        "total": 25,
        "impact": 0.190,
        "members": {
            "Lysine", "Saccharopine", "Allysine", "Glutarate",
            "2-Aminoadipate", "2-Oxoadipate", "Glutaryl-CoA",
            "Crotonyl-CoA", "Acetoacetyl-CoA", "Acetyl-CoA",
            "N6-Acetyl-L-lysine", "Pipecolate", "1-Piperideine-6-carboxylate",
            "2-Aminoadipate semialdehyde", "5-Aminovalerate",
            "Glutarate semialdehyde", "2-Hydroxyglutarate",
            "3-Methylglutaconyl-CoA", "3-Hydroxy-3-methylglutaryl-CoA",
            "Trimethyllysine", "4-Trimethylaminobutyraldehyde",
            "gamma-Butyrobetaine", "Carnitine", "L-Carnitine",
            "N6-Methyl-L-lysine",
        },
    },
    "hsa00330": {
        "name": "Arginine and proline metabolism",
        "total": 38,
        "impact": 0.256,
        "members": {
            "Arginine", "Proline", "Ornithine", "Urea", "Citrulline",
            "Putrescine", "Spermidine", "Spermine", "Creatine",
            "Creatinine", "Guanidinoacetate", "Methylguanidine",
            "Glutamate", "Glutamine", "4-Aminobutanoate", "Succinic acid",
            "2-Oxoglutarate", "Pyroglutamate", "Hydroxyproline",
            "Agmatine", "Carbamoyl phosphate", "Fumarate",
            "N-Acetylornithine", "N-Acetylcitrulline",
            "N-Acetyl-L-glutamate", "N-Acetylglutamate semialdehyde",
            "N-Acetylglutamate", "1-Pyrroline-5-carboxylate",
            "3-Hydroxy-L-proline", "trans-4-Hydroxy-L-proline",
            "Dimethylarginine", "Homocysteine", "Methionine",
            "Phosphocreatine", "1-Carboxyethyl-L-arginine",
            "5-Guanidino-2-oxopentanoate", "Nitric oxide", "Uric acid",
        },
    },
    "hsa00340": {
        "name": "Histidine metabolism",
        "total": 16,
        "impact": 0.137,
        "members": {
            "Histidine", "Urocanate", "Imidazolone propionate",
            "Formiminoglutamate", "Glutamate", "4-Aminobutanoate",
            "Histamine", "N-Formimidoylglutamate", "Imidazolepyruvate",
            "Imidazoleacetaldehyde", "Imidazoleacetate", "Carnosine",
            "Anserine", "Formic acid", "Glutamine", "1-Methylhistamine",
        },
    },
    "hsa00350": {
        "name": "Tyrosine metabolism",
        "total": 42,
        "impact": 0.360,
        "members": {
            "Tyrosine", "4-Hydroxyphenylacetic acid",
            "4-Hydroxyphenylpyruvate", "Homogentisate",
            "Maleylacetoacetate", "Fumarylacetoacetate",
            "Fumarate", "Acetoacetate", "Dopamine",
            "3,4-Dihydroxyphenylacetic acid", "3,4-Dihydroxyphenylethanol",
            "Vanillylmandelic acid", "Homovanillic acid", "3-Methoxytyramine",
            "Epinephrine", "Norepinephrine", "L-DOPA",
            "3-Methoxy-4-hydroxyphenylglycol", "Normetanephrine",
            "Phenylpyruvate", "Phenylacetate", "Phenylalanine",
            "4-Aminophenol", "4-Methylcatechol", "Catechol",
            "3,4-Dihydroxymandelate", "Melanin", "3-Methoxytyrosine",
            "3,4-Dihydroxyphenylalanine", "Phenylacetylglycine",
            "Phenylalanine hydroxylase", "Acetyl-CoA",
            "Succinic acid", "Succinylacetone", "Succinylacetoacetate",
            "4-Fumarylacetoacetate", "trans-Aconitic acid",
            "Hydroxymandelic acid", "Mandelate", "4-Coumarate", "Caffeate",
            "Ferulate",
        },
    },
    "hsa00360": {
        "name": "Phenylalanine metabolism",
        "total": 17,
        "impact": 0.186,
        "members": {
            "Phenylalanine", "Tyrosine", "Phenylpyruvate",
            "4-Hydroxyphenylacetic acid", "Phenylacetate",
            "Phenylacetylglycine", "Hippuric acid",
            "Phenylacetylglutamine", "trans-Cinnamate",
            "Phenyllactate", "Mandelate", "Benzaldehyde",
            "Benzyl alcohol", "Benzoate", "Benzoic acid",
            "4-Hydroxybenzoate", "2-Phenylacetamide",
        },
    },
    "hsa00380": {
        "name": "Tryptophan metabolism",
        "total": 41,
        "impact": 0.282,
        "members": {
            "Tryptophan", "Indole", "Indoxyl sulfate", "Indoxyl",
            "Kynurenine", "Anthranilate", "Quinolinate",
            "Nicotinamide", "NAD+", "NADP+",
            "3-Hydroxykynurenine", "3-Hydroxyanthranilate",
            "2-Amino-3-carboxymuconate", "Aminomuconate",
            "Glutarate", "2-Oxoglutarate", "Acetyl-CoA",
            "Acetoacetyl-CoA", "Serotonin", "N-Acetylserotonin",
            "Melatonin", "5-Hydroxyindoleacetaldehyde",
            "5-Hydroxyindoleacetate", "5-Hydroxytryptophan",
            "Indoleacetaldehyde", "Indoleacetate", "Skatole",
            "Tryptamine", "Indolepyruvate", "Indolelactate",
            "N-Formylkynurenine", "Picolinate", "Xanthurenate",
            "Kynurenate", "4,6-Dihydroxyquinolinate", "Acrylate",
            "Fumarate", "Maleate", "Formic acid", "Pyruvate",
            "Alanine",
        },
    },
    "hsa00400": {
        "name": "Phenylalanine, tyrosine and tryptophan biosynthesis",
        "total": 17,
        "impact": 0.131,
        "members": {
            "Phenylalanine", "Tyrosine", "Tryptophan",
            "Prephenate", "Phenylpyruvate", "4-Hydroxyphenylpyruvate",
            "Chorismate", "Anthranilate", "Indole",
            "Erythrose 4-phosphate", "Phosphoenolpyruvate",
            "3-Dehydroquinate", "Shikimate", "Shikimate 3-phosphate",
            "5-Enolpyruvylshikimate 3-phosphate", "3-Dehydroshikimate",
            "L-Arogenate",
        },
    },
    "hsa00410": {
        "name": "beta-Alanine metabolism",
        "total": 22,
        "impact": 0.196,
        "members": {
            "Beta-alanine", "3-Aminoisobutyric acid", "Pantothenic acid",
            "4-Aminobutanoate", "Aspartate", "Malonate semialdehyde",
            "Acetyl-CoA", "Uracil", "Dihydrouracil",
            "3-Ureidopropionate", "Methylmalonate semialdehyde",
            "Propanoyl-CoA", "Alanine", "Pantothenol",
            "Pantetheine", "CoA", "3-Phosphoadenosine-5-phosphosulfate",
            "Succinic acid", "Malonate", "Propionamide",
            "Propionate", "3-Aminopropanal",
        },
    },
    "hsa00430": {
        "name": "Taurine and hypotaurine metabolism",
        "total": 10,
        "impact": 0.237,
        "members": {
            "Taurine", "Hypotaurine", "Cystamine",
            "Cystathionine", "Cysteine", "Cysteamine",
            "Pantetheine", "Coenzyme A", "Sulfite", "Sulfate",
        },
    },
    "hsa00470": {
        "name": "D-Amino acid metabolism",
        "total": 20,
        "impact": 0.098,
        "members": {
            "D-Alanine", "D-Serine", "D-Aspartate", "D-Glutamate",
            "D-Proline", "D-Leucine", "D-Isoleucine", "D-Valine",
            "D-Phenylalanine", "D-Methionine", "D-Tryptophan",
            "D-Tyrosine", "D-Threonine", "D-Arginine", "D-Lysine",
            "D-Histidine", "D-Glutamine", "N-Acetylaspartic acid",
            "N-Acetyl-D-glucosamine", "N-Acetyl-L-aspartylglutamate",
        },
    },
    "hsa00480": {
        "name": "Glutathione metabolism",
        "total": 28,
        "impact": 0.183,
        "members": {
            "Glutathione", "Glutamine", "Glycine", "Cysteine",
            "Glutamate", "Gamma-glutamylcysteine", "Cysteinylglycine",
            "Oxidized glutathione", "5-Oxoproline", "Glutathione S-conjugate",
            "N,N-Dimethylglycine", "Sarcosine", "Betaine", "Choline",
            "S-Formylglutathione", "Formate", "Formic acid",
            "Leukotriene C4", "Leukotriene D4", "Leukotriene E4",
            "Prostaglandin H2", "Thioredoxin", "Glutaredoxin",
            "Ascorbate", "Dehydroascorbate", "Uric acid",
            "Homocysteine", "Cystine",
        },
    },
    "hsa00500": {
        "name": "Starch and sucrose metabolism",
        "total": 25,
        "impact": 0.143,
        "members": {
            "Sucrose", "Glucose", "Fructose", "Glucose 6-phosphate",
            "Fructose 6-phosphate", "Glucose 1-phosphate", "Maltose",
            "Lactose", "UDP-glucose", "ADP-glucose", "Trehalose",
            "Cellobiose", "Amylose", "Amylopectin", "Glycogen",
            "UDP-galactose", "Galactose", "1,3-beta-D-Glucan",
            "Dextrin", "Cyclomaltodextrin", "Sucrose 6-phosphate",
            "Fructose 1,6-bisphosphate", "Glucono-1,5-lactone",
            "6-Phosphogluconate", "Raffinose",
        },
    },
    "hsa00590": {
        "name": "Arachidonic acid metabolism",
        "total": 36,
        "impact": 0.241,
        "members": {
            "Arachidonic acid", "Prostaglandin G2", "Prostaglandin H2",
            "Prostaglandin E2", "Prostaglandin D2", "Prostaglandin F2a",
            "Thromboxane A2", "Thromboxane B2", "Prostacyclin",
            "12-HPETE", "15-HPETE", "5-HPETE", "5-HETE", "12-HETE",
            "15-HETE", "Leukotriene A4", "Leukotriene B4",
            "Leukotriene C4", "Leukotriene D4", "Leukotriene E4",
            "Lipoxin A4", "Lipoxin B4", "EET",
            "Acetylsalicylic acid",  # aspirin inhibits COX in this pathway
            "Indomethacin", "Ibuprofen", "2-Arachidonoylglycerol",
            "Arachidonoylethanolamide", "Anandamide",
            "Dihomo-gamma-linolenate", "Eicosadienoate",
            "Docosahexaenoate", "Eicosapentaenoate",
            "8-Oxo-ETE", "Hepoxilin A3", "20-HETE",
        },
    },
    "hsa00620": {
        "name": "Pyruvate metabolism",
        "total": 22,
        "impact": 0.209,
        "members": {
            "Pyruvate", "Acetyl-CoA", "Acetaldehyde", "Ethanol",
            "Lactic acid", "Oxaloacetate", "Malate", "Formic acid",
            "Acetic acid", "Acetate", "Propanoyl-CoA", "Acetoacetate",
            "Methylglyoxal", "Alanine", "2-Hydroxybutyrate",
            "2-Hydroxyisobutyric acid", "Pantothenic acid",
            "(S)-Malate", "Citric acid", "Succinic acid",
            "D-Lactate", "Succinyl-CoA",
        },
    },
    "hsa00630": {
        "name": "Glyoxylate and dicarboxylate metabolism",
        "total": 20,
        "impact": 0.174,
        "members": {
            "Glyoxylate", "Oxalate", "Glycolate", "Glycine",
            "Formic acid", "Succinic acid", "Fumarate", "Malate",
            "Citric acid", "Isocitrate", "2-Oxoglutarate",
            "trans-Aconitic acid", "Tartrate", "Hydroxypyruvate",
            "Serine", "Threonine", "Glucuronate", "Gluconate",
            "Gluconic acid", "Allantoin",
        },
    },
    "hsa00640": {
        "name": "Propanoate metabolism",
        "total": 18,
        "impact": 0.208,
        "members": {
            "Propanoyl-CoA", "Methylmalonyl-CoA", "Succinyl-CoA",
            "Propionate", "Acetic acid", "Acrylyl-CoA",
            "3-Hydroxyisobutyric acid", "Methylmalonate",
            "2-Aminobutyric acid", "Propylene glycol",
            "Lactaldehyde", "Propanoyl-phosphate",
            "3-Methyl-2-oxopentanoic acid",
            "Succinic acid", "Valine", "Isoleucine",
            "Threonine", "Methionine",
        },
    },
    "hsa00650": {
        "name": "Butanoate metabolism",
        "total": 20,
        "impact": 0.198,
        "members": {
            "Butyrate", "Acetoacetate", "Acetoacetyl-CoA",
            "Acetyl-CoA", "Butyryl-CoA", "Crotonyl-CoA",
            "3-Hydroxybutyryl-CoA", "4-Aminobutanoate",
            "Succinate semialdehyde", "Succinic acid",
            "2-Oxoglutarate", "Glutamate", "Glutamine",
            "Citric acid", "3-Methylbutanoate", "2-Methylglutaric acid",
            "3-Hydroxybutyrate", "Acetone",
            "Isobutyrate", "4-Hydroxybutyrate",
        },
    },
    "hsa00670": {
        "name": "One carbon pool by folate",
        "total": 12,
        "impact": 0.148,
        "members": {
            "Formic acid", "5,10-Methenyl-THF", "5,10-Methylene-THF",
            "5-Methyl-THF", "10-Formyl-THF", "THF", "Dihydrofolate",
            "5-Formyl-THF", "5-Formimino-THF",
            "Glycine", "Serine", "Methionine",
        },
    },
    "hsa00760": {
        "name": "Nicotinate and nicotinamide metabolism",
        "total": 22,
        "impact": 0.217,
        "members": {
            "Nicotinate", "Nicotinamide", "NAD+", "NADP+",
            "1-Methylnicotinamide", "Trigonelline",
            "Nicotinate D-ribonucleotide", "Nicotinate D-ribonucleoside",
            "Nicotinamide D-ribonucleotide", "Nicotinamide riboside",
            "Quinolinate", "Tryptophan", "2-Amino-3-carboxymuconate",
            "Aminocarboxymuconate semialdehyde", "Picolinate",
            "Nicotinamide N-oxide", "N-Methyl-2-pyridone-5-carboxamide",
            "N-Methyl-4-pyridone-3-carboxamide",
            "N1-Methyl-4-pyridinium-3-carboxamide",
            "ADP-ribose", "Cyclic ADP-ribose", "2-Pyridone",
        },
    },
    "hsa00770": {
        "name": "Pantothenate and CoA biosynthesis",
        "total": 15,
        "impact": 0.162,
        "members": {
            "Pantothenic acid", "Pantothenol", "Pantetheine",
            "CoA", "Phosphopantothenate", "Phosphopantothenoylcysteine",
            "Phosphopantetheine", "Dephospho-CoA",
            "Beta-alanine", "3-Aminoisobutyric acid",
            "4-Phosphopantothenate", "Pantoate", "2-Dehydropantoate",
            "Valine", "2-Oxoisovalerate",
        },
    },
    "hsa00780": {
        "name": "Biotin metabolism",
        "total": 10,
        "impact": 0.124,
        "members": {
            "Biotin", "Biocytin", "Desthiobiotin", "7-Keto-8-aminopelargonate",
            "7,8-Diaminopelargonate", "Alanine", "Pimeloyl-CoA",
            "Pimeloyl-ACP", "CoA", "Tyrosine",
        },
    },
    "hsa01200": {
        "name": "Carbon metabolism",
        "total": 118,
        "impact": 0.115,
        "members": {
            "Glucose", "Glucose 6-phosphate", "Fructose 6-phosphate",
            "Fructose 1,6-bisphosphate", "Glyceraldehyde 3-phosphate",
            "Pyruvate", "Acetyl-CoA", "Citric acid", "Isocitrate",
            "2-Oxoglutarate", "Succinyl-CoA", "Succinic acid",
            "Fumarate", "Malate", "Oxaloacetate", "Lactic acid",
            "Alanine", "Aspartate", "Glutamate", "Glutamine",
            "Glycine", "Serine", "Threonine", "Cysteine", "Methionine",
            "Valine", "Leucine", "Isoleucine", "Phenylalanine", "Tyrosine",
            "Tryptophan", "Lysine", "Arginine", "Proline",
            "trans-Aconitic acid", "Ribose 5-phosphate",
            "Ribulose 5-phosphate", "Erythrose 4-phosphate",
            "Sedoheptulose 7-phosphate", "6-Phosphogluconate",
            "Gluconic acid", "Formic acid", "Acetic acid",
            "Propionate", "Butyrate", "3-Phosphoglycerate",
            "2-Phosphoglycerate", "Phosphoenolpyruvate",
            "Dihydroxyacetone phosphate", "1,3-Bisphosphoglycerate",
            "Acetaldehyde", "Ethanol", "Glyoxylate", "Glycolate",
            "Hydroxypyruvate", "Phosphoserine", "Homoserine",
            "5,10-Methylene-THF", "5-Methyl-THF", "THF",
            "N,N-Dimethylglycine", "Sarcosine", "Betaine",
            "S-Adenosylmethionine", "Homocysteine", "Cystathionine",
            "2-Aminobutyric acid", "2-Oxobutanoate",
            "Citramalic acid", "Taurine",
        },
    },
    "hsa01210": {
        "name": "2-Oxocarboxylic acid metabolism",
        "total": 25,
        "impact": 0.189,
        "members": {
            "Pyruvate", "2-Oxoglutarate", "Oxaloacetate",
            "2-Oxobutanoate", "2-Oxoisovalerate", "2-Oxoisocaproate",
            "3-Methyl-2-oxopentanoic acid",
            "4-Methylthio-2-oxobutanoate", "3-Indolepyruvate",
            "Phenylpyruvate", "4-Hydroxyphenylpyruvate",
            "Glyoxylate", "Oxalate", "Succinic acid",
            "Citric acid", "Isocitrate", "Alanine", "Glycine",
            "Valine", "Leucine", "Isoleucine", "Phenylalanine",
            "Tyrosine", "Tryptophan", "Methionine",
        },
    },
}

# Total number of unique compounds in KEGG human metabolic pathways (background)
BACKGROUND_SIZE = 800

RESULTS_DIR = "pathway_results"


# ── 3. ORA via Fisher's exact test ────────────────────────────────────────────

def run_ora(compound_list: list[str], group_name: str) -> list[dict]:
    """
    Over-representation analysis for one group.
    Returns list of pathway result dicts, sorted by raw p-value.
    """
    query_set = set(compound_list)
    n_query = len(query_set)  # total input metabolites

    rows = []
    for pid, pw in KEGG_PATHWAYS.items():
        pathway_members = pw["members"]
        hits_in_query = query_set & pathway_members
        k = len(hits_in_query)        # hits in query ∩ pathway
        K = len(pathway_members)      # pathway total
        n = n_query                   # query list size
        N = BACKGROUND_SIZE           # background universe

        # Fisher's exact test (one-sided, over-representation)
        # Contingency table:
        #                 In pathway    Not in pathway
        # In query            k             n - k
        # Not in query      K - k    N - K - (n - k)
        table = [
            [k, n - k],
            [K - k, N - K - (n - k)],
        ]
        if table[1][1] < 0:
            # Edge case: query larger than background
            continue
        _, raw_p = fisher_exact(table, alternative="greater")

        rows.append({
            "Group": group_name,
            "Pathway_ID": pid,
            "Pathway": pw["name"],
            "Total": K,
            "Hits": k,
            "Raw_P": raw_p,
            "FDR": None,  # filled after BH correction
            "Impact": pw["impact"],
            "Hits_metabolites_name": "; ".join(sorted(hits_in_query)),
            "Significant": "No",
        })

    # Benjamini-Hochberg FDR across pathways within this group
    if rows:
        pvals = [r["Raw_P"] for r in rows]
        _, fdr_vals, _, _ = multipletests(pvals, method="fdr_bh")
        for row, fdr in zip(rows, fdr_vals):
            row["FDR"] = fdr
            row["Significant"] = "Yes" if row["Raw_P"] < 0.05 else "No"

    rows.sort(key=lambda r: r["Raw_P"])
    return rows


# ── 4. Main ───────────────────────────────────────────────────────────────────

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    all_rows = []

    for group_name, compounds in GROUPS.items():
        print(f"\n{'='*60}")
        print(f"Group: {group_name}  ({len(compounds)} compounds)")
        print(f"{'='*60}")

        # Show which compounds were matched in any pathway
        all_pathway_members = set()
        for pw in KEGG_PATHWAYS.values():
            all_pathway_members |= pw["members"]
        matched = [c for c in compounds if c in all_pathway_members]
        unmatched = [c for c in compounds if c not in all_pathway_members]
        print(f"  Matched to KEGG pathways : {len(matched)}/{len(compounds)}")
        if unmatched:
            print(f"  Not found in database   : {', '.join(unmatched)}")

        rows = run_ora(compounds, group_name)
        all_rows.extend(rows)

        # Save per-group JSON
        group_path = os.path.join(RESULTS_DIR, f"{group_name}_pathway.json")
        with open(group_path, "w") as f:
            json.dump(rows, f, indent=2)

        sig = [r for r in rows if r["Significant"] == "Yes"]
        print(f"  Pathways analysed        : {len(rows)}")
        print(f"  Significant (p<0.05)     : {len(sig)}")
        for r in sig:
            print(f"    {r['Pathway']:<50} p={r['Raw_P']:.4f}  FDR={r['FDR']:.4f}  hits={r['Hits']}  impact={r['Impact']:.3f}")

    # ── Summary CSV ───────────────────────────────────────────────────────────
    df = pd.DataFrame(all_rows, columns=[
        "Group", "Pathway_ID", "Pathway", "Total", "Hits",
        "Raw_P", "FDR", "Impact", "Hits_metabolites_name", "Significant",
    ])
    df.sort_values(["Group", "Raw_P"], inplace=True)
    df.reset_index(drop=True, inplace=True)

    summary_path = os.path.join(RESULTS_DIR, "pathway_summary.csv")
    df.to_csv(summary_path, index=False)

    # Significant-only CSV for quick reference
    sig_df = df[df["Significant"] == "Yes"].copy()
    sig_path = os.path.join(RESULTS_DIR, "pathway_significant.csv")
    sig_df.to_csv(sig_path, index=False)

    print(f"\n{'='*60}")
    print(f"Summary saved      : {summary_path}")
    print(f"Significant only   : {sig_path}")
    print(f"Total rows         : {len(df)}")
    print(f"Significant (p<0.05): {len(sig_df)}")

    if len(sig_df):
        print("\nAll significant pathways:")
        print(sig_df[["Group", "Pathway", "Total", "Hits", "Raw_P", "FDR", "Impact", "Hits_metabolites_name"]].to_string(index=False))


if __name__ == "__main__":
    main()
