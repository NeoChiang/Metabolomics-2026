#!/usr/bin/env python3
"""
Create Metabolomics 2026 conference poster (A0 portrait, 84.1 × 118.9 cm).
Layout:
  - Header (navy): title / authors / affiliations
  - Text zone (3 cols): Introduction+Methods | Key Findings 1-4 | Findings 5-6+Table 1
  - Figure zone (full-width split): Figure 1A GA heatmap | Figure 1B BPD heatmap
  - Pathway zone: Table 2 pathway analysis | Conclusions
  - Footer (navy): references + contact
Output: poster/Metabolomics2026_poster.pptx
"""

import os
from pptx import Presentation
from pptx.util import Cm, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from lxml import etree
from pptx.oxml.ns import qn

# ── Palette ───────────────────────────────────────────────────────────────────
NAVY     = RGBColor(0x00, 0x2B, 0x5C)
WHITE    = RGBColor(0xFF, 0xFF, 0xFF)
BLACK    = RGBColor(0x1A, 0x1A, 0x1A)
TEAL     = RGBColor(0x00, 0x7B, 0xA7)
DGRAY    = RGBColor(0x2C, 0x3E, 0x50)
LGRAY    = RGBColor(0xEC, 0xF0, 0xF1)
THEAD    = RGBColor(0x00, 0x2B, 0x5C)
TROW_A   = RGBColor(0xD6, 0xEA, 0xF8)
TROW_B   = RGBColor(0xD5, 0xF5, 0xE3)
TROW_C   = RGBColor(0xFD, 0xF2, 0xCC)
TROW_D   = RGBColor(0xFA, 0xD7, 0xA0)

# ── Geometry ──────────────────────────────────────────────────────────────────
A0_W, A0_H = Cm(84.1), Cm(118.9)
MAR  = Cm(2.0)
GUT  = Cm(1.0)
CW   = Cm(26.0)
CX1, CX2, CX3 = MAR, MAR + CW + GUT, MAR + 2 * (CW + GUT)

HDR_H   = Cm(14.0)
TEXT_Y  = HDR_H + Cm(0.5)           # 14.5
TEXT_H  = Cm(41.0)                  # text zone height
FIG_Y   = TEXT_Y + TEXT_H + Cm(0.5) # 56.0
FIG_H   = Cm(33.5)                  # figure zone height
PTH_Y   = FIG_Y + FIG_H + Cm(0.5)  # 90.0
PTH_H   = Cm(16.5)
FTR_Y   = PTH_Y + PTH_H             # 106.5
FTR_H   = A0_H - FTR_Y             # ~12.4

POSTER_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT   = os.path.join(POSTER_DIR, "Metabolomics2026_poster.pptx")
IMG_GA   = os.path.join(POSTER_DIR, "GA_heatmap_v6_cre.png")
IMG_BPD  = os.path.join(POSTER_DIR, "BPD_heatmap_v6_cre.png")

# GA: 1800×1400 → ratio w/h = 1.286  →  h = w / 1.286
# BPD: 1800×1600 → ratio w/h = 1.125  →  h = w / 1.125
GA_RATIO  = 1800 / 1400
BPD_RATIO = 1800 / 1600


# ── Low-level helpers ─────────────────────────────────────────────────────────

def rect(slide, x, y, w, h, fill=None, line=False, line_color=None):
    shp = slide.shapes.add_shape(1, int(x), int(y), int(w), int(h))
    if fill:
        shp.fill.solid()
        shp.fill.fore_color.rgb = fill
    else:
        shp.fill.background()
    if line and line_color:
        shp.line.color.rgb = line_color
        shp.line.width = Pt(0.75)
    else:
        shp.line.fill.background()
    return shp


def txbox(slide, x, y, w, h):
    tb = slide.shapes.add_textbox(int(x), int(y), int(w), int(h))
    tb.word_wrap = True
    tf = tb.text_frame
    tf.word_wrap = True
    return tb, tf


def para(tf, text, size=22, bold=False, italic=False, color=BLACK,
         align=PP_ALIGN.LEFT, space_before=0, space_after=0, first=False,
         underline=False):
    p = tf.paragraphs[0] if first else tf.add_paragraph()
    p.alignment = align
    if space_before:
        p.space_before = Pt(space_before)
    if space_after:
        p.space_after = Pt(space_after)
    if text or text == "":
        run = p.add_run()
        run.text = text
        run.font.size = Pt(size)
        run.font.bold = bold
        run.font.italic = italic
        run.font.color.rgb = color
        if underline:
            run.font.underline = True
    return p


def heading(slide, text, x, y, w, size=30, color=NAVY):
    """Section heading with teal underline. Returns next y."""
    tb, tf = txbox(slide, x, y, w, Cm(1.8))
    para(tf, text.upper(), size=size, bold=True, color=color, first=True)
    rect(slide, x, y + Cm(1.5), w, Cm(0.15), fill=TEAL)
    return y + Cm(2.0)


def finding(slide, cx, cy, num_title, body, cw=CW):
    """Add a Key Finding block. Returns next y."""
    tb, tf = txbox(slide, cx, cy, cw, Cm(1.0))
    para(tf, num_title, size=20, bold=True, color=NAVY, first=True)
    cy += Cm(0.95)
    nchars = len(body)
    est_lines = max(2, nchars // 72 + 1)
    tb, tf = txbox(slide, cx + Cm(0.4), cy, cw - Cm(0.4), Cm(est_lines * 0.85 + 0.3))
    para(tf, body, size=18, color=DGRAY, first=True)
    cy += Cm(est_lines * 0.75 + 0.5)
    return cy


def set_cell_bg(cell, color):
    tc = cell._tc
    tcPr = tc.get_or_add_tcPr()
    for old in tcPr.findall(qn("a:solidFill")):
        tcPr.remove(old)
    sf = etree.SubElement(tcPr, qn("a:solidFill"))
    clr = etree.SubElement(sf, qn("a:srgbClr"))
    clr.set("val", str(color))


def set_cell_text(cell, text, size=16, bold=False, italic=False,
                  color=BLACK, align=PP_ALIGN.LEFT):
    tf = cell.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.alignment = align
    p.clear()
    run = p.add_run()
    run.text = text
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = color


# ── Build poster ──────────────────────────────────────────────────────────────

def main():
    prs = Presentation()
    prs.slide_width  = A0_W
    prs.slide_height = A0_H

    slide = prs.slides.add_slide(prs.slide_layouts[6])  # blank layout
    for shp in list(slide.shapes):
        shp._element.getparent().remove(shp._element)

    # White background
    rect(slide, 0, 0, A0_W, A0_H, fill=WHITE)

    # ──────────────────────────────────────────────────────────────────────────
    # HEADER
    # ──────────────────────────────────────────────────────────────────────────
    rect(slide, 0, 0, A0_W, HDR_H, fill=NAVY)

    # Accent stripe at header bottom
    rect(slide, 0, HDR_H - Cm(0.4), A0_W, Cm(0.4), fill=TEAL)

    # Title
    tb, tf = txbox(slide, Cm(2), Cm(0.8), Cm(78), Cm(5.5))
    para(tf,
         "Longitudinal Urinary Metabolomics in Preterm Infants",
         size=52, bold=True, color=WHITE, align=PP_ALIGN.CENTER, first=True)
    para(tf,
         "Associated with Gestational Age and BPD Severity from 6 Months to 2 Years of Age",
         size=44, bold=True, color=WHITE, align=PP_ALIGN.CENTER)

    # Authors
    tb, tf = txbox(slide, Cm(2), Cm(7.0), Cm(78), Cm(1.8))
    para(tf,
         "Meng-Han Chiang¹  │  Chih-Yung Chiu¹²*",
         size=30, bold=True, color=WHITE, align=PP_ALIGN.CENTER, first=True)

    # Affiliations
    tb, tf = txbox(slide, Cm(2), Cm(9.0), Cm(78), Cm(3.2))
    para(tf,
         "¹Multi-omics Core Laboratory, Chang Gung Memorial Hospital at Linkou, "
         "Taoyuan, Taiwan",
         size=21, italic=True, color=WHITE, align=PP_ALIGN.CENTER, first=True)
    para(tf,
         "²Division of Pediatric Pulmonology, Department of Pediatrics, "
         "Chang Gung Memorial Hospital at Linkou, and Chang Gung University, Taiwan",
         size=21, italic=True, color=WHITE, align=PP_ALIGN.CENTER)
    para(tf,
         "*Corresponding author: pedchestIC@gmail.com  │  Metabolomics 2026",
         size=19, color=RGBColor(0xAA, 0xCC, 0xFF), align=PP_ALIGN.CENTER)

    # CGMH label (top-right)
    tb, tf = txbox(slide, Cm(72), Cm(1.0), Cm(10), Cm(4.0))
    para(tf, "CGMH\nLINKOU", size=24, bold=True, color=WHITE,
         align=PP_ALIGN.CENTER, first=True)

    # ──────────────────────────────────────────────────────────────────────────
    # TEXT ZONE — 3 columns
    # ──────────────────────────────────────────────────────────────────────────

    # ── Column 1: Introduction + Methods ─────────────────────────────────────
    cx, cy = CX1, TEXT_Y

    cy = heading(slide, "Introduction", cx, cy, CW, size=28)

    intro = (
        "Bronchopulmonary dysplasia (BPD) is a chronic lung disease predominantly "
        "affecting preterm infants requiring mechanical ventilation or supplemental oxygen. "
        "Despite clinical improvement within the first two years, the longitudinal trajectory "
        "of urinary metabolic alterations — and their associations with gestational age (GA) "
        "and BPD severity — remains insufficiently characterized.\n\n"
        "This study identified urinary metabolomic profiles at 6 months (6M) and 2 years (2Y) "
        "of corrected age across GA and BPD severity groups using ¹H-NMR spectroscopy."
    )
    tb, tf = txbox(slide, cx, cy, CW, Cm(13.5))
    para(tf, intro, size=19, color=DGRAY, first=True)
    cy += Cm(14.0)

    cy = heading(slide, "Methods", cx, cy, CW, size=28)

    method_items = [
        ("Study Population:",
         "140 preterm (<32 wk) and full-term (≥37 wk) infants prospectively enrolled "
         "(2019 NICHD BPD criteria). 254 urine samples: 6M (n=139), 2Y (n=115)."),
        ("Stratification:",
         "GA: ≥37 wk (controls), 28–32 wk, <28 wk\n"
         "BPD: Healthy controls (HC), No+Mild BPD, Moderate+Severe (M+S) BPD"),
        ("¹H-NMR Processing:",
         "Spectra processed via NMRProcFlow. Metabolites identified with Chenomx 8.1. "
         "Intensities normalized to urinary creatinine."),
        ("Statistics:",
         "Kruskal-Wallis test • Mfuzz clustering • MetaboAnalyst pathway ORA "
         "(Fisher’s exact + Benjamini-Hochberg FDR)"),
    ]

    for label, text in method_items:
        tb, tf = txbox(slide, cx, cy, CW, Cm(0.9))
        para(tf, label, size=19, bold=True, color=NAVY, first=True)
        cy += Cm(0.85)
        nchars = len(text)
        est = max(1, nchars // 65 + 1)
        tb, tf = txbox(slide, cx + Cm(0.5), cy, CW - Cm(0.5),
                       Cm(est * 0.8 + 0.3))
        para(tf, text, size=18, color=DGRAY, first=True)
        cy += Cm(est * 0.72 + 0.45)

    # ── Column 2: Key Findings 1–4 ────────────────────────────────────────────
    cx, cy = CX2, TEXT_Y

    cy = heading(slide, "Key Findings", cx, cy, CW, size=28)

    cy = finding(slide, cx, cy,
        "1 — Expanding metabolic burden",
        "16 and 18 metabolites associated with GA, and 30 and 12 with BPD severity, "
        "at 6M and 2Y (P<0.05). Marked contraction of BPD-associated metabolites "
        "from 6M to 2Y suggests a shifting metabolic landscape despite clinical improvement.")
    cy = finding(slide, cx, cy,
        "2 — Gut microbiome signature at 6M",
        "3-MOV, 4-hydroxyphenylacetate, indoxyl sulfate, and N-PAG significantly elevated "
        "in GA<28 wk vs. full-term (P<0.001); also elevated in M+S BPD vs. HC (P<0.001), "
        "pointing to early gut dysbiosis in extremely preterm infants.")
    cy = finding(slide, cx, cy,
        "3 — Persistence and amplification to 2Y",
        "N-PAG and indoxyl sulfate remained elevated at 2Y with higher fold changes "
        "(N-PAG: 1.82 vs 1.76; indoxyl sulfate: 1.73 vs 1.57), suggesting metabolic "
        "divergence intensifies over time despite apparent clinical recovery.")
    cy = finding(slide, cx, cy,
        "4 — Complete metabolic overlap at 6M",
        "All 16 GA-significant metabolites at 6M were also BPD-significant — zero "
        "GA-exclusive metabolites. Prematurity and BPD severity share an indistinguishable "
        "urinary metabolic signature at 6 months corrected age.")

    # ── Column 3: Findings 5-6 + Table 1 (GA demographics) ───────────────────
    cx, cy = CX3, TEXT_Y

    cy = heading(slide, "Key Findings (cont.)", cx, cy, CW, size=28)

    cy = finding(slide, cx, cy,
        "5 — Consistent biomarker panel",
        "N-PAG, indoxyl sulfate, acetylsalicylate, glutamine, pantothenic acid, and "
        "succinate were associated with GA and/or BPD at both 6M and 2Y — the most "
        "temporally consistent metabolic alterations identified.")
    cy = finding(slide, cx, cy,
        "6 — Persistent pathway disruption",
        "Alanine/aspartate/glutamate (FDR<0.05 at 6M) and glyoxylate/dicarboxylate "
        "metabolism remained enriched at both timepoints. Valine/leucine/isoleucine "
        "biosynthesis enriched at both 6M and 2Y, consistent with 3-MOV signal.")

    cy += Cm(0.4)
    cy = heading(slide, "Table 1 — Demographics (GA Groups)", cx, cy, CW, size=24)

    tbl1_rows = [
        ["Characteristic",      "≥37 wk\n(n=50/40)",  "28–32 wk\n(n=48/37)", "<28 wk\n(n=41/38)", "P"],
        ["GA (wk)",              "39.1±0.8 / 39.0±0.9", "30.7±0.9 / 30.8±0.9", "26.1±1.4 / 26.2±1.3", "<0.001"],
        ["Birth wt (g)",         "3136±348 / 3147±385", "1352±276 / 1374±224", "782±205 / 772±183",   "<0.001"],
        ["Body wt (kg)",         "7.9±0.8 / 12.3±1.6",  "7.8±1.4 / 11.4±1.7",  "7.0±1.5 / 10.6±1.5",  "0.003 / <0.001"],
        ["BMI (kg/m²)",     "17.2±1.1 / 16.0±1.3", "17.7±2.4 / 15.6±1.4", "16.5±1.8 / 14.9±1.4", "0.019 / 0.002"],
        ["Breastfed ≥6M (%)", "56.8 / 60.0",           "30.0 / 43.2",           "44.7 / 44.7",           "0.047 / 0.261"],
        ["Sepsis (%)",           "2.3 / 0.0",             "30.0 / 32.4",           "81.1 / 81.1",           "<0.001"],
    ]

    n_rows, n_cols = len(tbl1_rows), 5
    col_ws1 = [Cm(7.5), Cm(4.8), Cm(4.8), Cm(4.8), Cm(4.1)]
    row_h   = Cm(1.6)
    tbl_h1  = n_rows * row_h

    tbl1 = slide.shapes.add_table(
        n_rows, n_cols, int(cx), int(cy), int(CW), int(tbl_h1)
    ).table
    for ci, cw in enumerate(col_ws1):
        tbl1.columns[ci].width = int(cw)
    for ri, row_data in enumerate(tbl1_rows):
        for ci, cell_text in enumerate(row_data):
            cell = tbl1.cell(ri, ci)
            set_cell_text(cell, cell_text,
                          size=14 if ri == 0 else 13,
                          bold=(ri == 0),
                          color=WHITE if ri == 0 else BLACK,
                          align=PP_ALIGN.CENTER if ci > 0 else PP_ALIGN.LEFT)
            if ri == 0:
                set_cell_bg(cell, THEAD)
            elif ri % 2 == 1:
                set_cell_bg(cell, LGRAY)
            else:
                set_cell_bg(cell, WHITE)

    cy += tbl_h1 + Cm(0.3)
    tb, tf = txbox(slide, cx, cy, CW, Cm(1.5))
    para(tf,
         "Values: mean±SD or % (6M / 2Y). Bold P < 0.05 is significant. "
         "GA, gestational age; BMI, body mass index.",
         size=15, italic=True, color=DGRAY, first=True)

    # ──────────────────────────────────────────────────────────────────────────
    # FIGURE ZONE — GA and BPD heatmaps side by side
    # ──────────────────────────────────────────────────────────────────────────
    fig_zone_y = FIG_Y

    # Thin accent bar separating text from figures
    rect(slide, MAR, fig_zone_y, A0_W - 2 * MAR, Cm(0.12), fill=TEAL)

    cap_y   = fig_zone_y + Cm(0.3)
    img_y   = cap_y + Cm(1.2)
    img_h   = Cm(30.0)                          # fixed height for both

    ga_w    = img_h * GA_RATIO                  # 30 × 1.286 = 38.6 cm
    bpd_w   = img_h * BPD_RATIO                 # 30 × 1.125 = 33.75 cm

    ga_x    = MAR
    bpd_x   = ga_x + ga_w + Cm(2.0)

    # Captions above each figure
    tb, tf = txbox(slide, ga_x, cap_y, ga_w, Cm(1.1))
    para(tf,
         "Figure 1A. Urinary metabolites significantly associated with gestational age (GA).",
         size=18, bold=True, color=NAVY, first=True)

    tb, tf = txbox(slide, bpd_x, cap_y, bpd_w, Cm(1.1))
    para(tf,
         "Figure 1B. Urinary metabolites significantly associated with BPD severity.",
         size=18, bold=True, color=NAVY, first=True)

    # Images
    slide.shapes.add_picture(IMG_GA,  int(ga_x),  int(img_y), int(ga_w),  int(img_h))
    slide.shapes.add_picture(IMG_BPD, int(bpd_x), int(img_y), int(bpd_w), int(img_h))

    # Thin accent bar below figures
    rect(slide, MAR, PTH_Y - Cm(0.2), A0_W - 2 * MAR, Cm(0.12), fill=TEAL)

    # ──────────────────────────────────────────────────────────────────────────
    # PATHWAY + CONCLUSIONS ZONE
    # ──────────────────────────────────────────────────────────────────────────
    pth_x  = MAR
    pth_w  = Cm(54.0)
    con_x  = pth_x + pth_w + GUT
    con_w  = A0_W - con_x - MAR

    pcy = PTH_Y + Cm(0.3)
    pcy = heading(slide, "Pathway Enrichment Analysis", pth_x, pcy, pth_w, size=26)

    pathway_rows = [
        ["Group",       "Pathway",                              "Hits / Total", "Raw P",   "FDR",    "Function"],
        ["6M BPD only", "Ala, Asp & Glu metabolism",           "4 / 28",       "<0.001",  "<0.001", "Amino acid"],
        ["6M BPD only", "Glyoxylate & dicarboxylate metabolism","3 / 32",       "<0.001",  "0.023",  "Carbohydrate"],
        ["6M GA+BPD",   "Val/Leu/Ile biosynthesis",            "2 / 8",        "<0.001",  "0.063",  "Amino acid"],
        ["2Y GA only",  "Tyrosine metabolism",                  "2 / 42",       "0.017",   "0.792",  "Amino acid"],
        ["2Y GA only",  "Taurine & hypotaurine metabolism",     "1 / 10",       "0.040",   "0.792",  "Other AA"],
        ["2Y GA+BPD",   "Ala, Asp & Glu metabolism",           "2 / 28",       "0.004",   "0.223",  "Amino acid"],
        ["2Y GA+BPD",   "Glyoxylate & dicarboxylate metabolism","2 / 32",       "0.006",   "0.223",  "Carbohydrate"],
    ]

    n_rows2 = len(pathway_rows)
    col_ws2 = [Cm(9.0), Cm(19.5), Cm(7.5), Cm(6.0), Cm(6.0), Cm(6.0)]
    row_h2  = Cm(1.65)
    tbl_h2  = n_rows2 * row_h2

    row_colors = [None, TROW_A, TROW_A, TROW_B, TROW_C, TROW_C, TROW_D, TROW_D]

    tbl2 = slide.shapes.add_table(
        n_rows2, 6, int(pth_x), int(pcy), int(pth_w), int(tbl_h2)
    ).table
    for ci, cw in enumerate(col_ws2):
        tbl2.columns[ci].width = int(cw)
    for ri, row_data in enumerate(pathway_rows):
        for ci, cell_text in enumerate(row_data):
            cell = tbl2.cell(ri, ci)
            set_cell_text(cell, cell_text,
                          size=15 if ri == 0 else 14,
                          bold=(ri == 0),
                          color=WHITE if ri == 0 else BLACK,
                          align=PP_ALIGN.CENTER if ci >= 2 else PP_ALIGN.LEFT)
            if ri == 0:
                set_cell_bg(cell, THEAD)
            elif row_colors[ri]:
                set_cell_bg(cell, row_colors[ri])

    pcy += tbl_h2 + Cm(0.3)
    tb, tf = txbox(slide, pth_x, pcy, pth_w, Cm(1.3))
    para(tf,
         "ORA: Fisher’s exact test, Benjamini-Hochberg FDR. Background: 800 KEGG human metabolic compounds. "
         "6M GA only: all GA metabolites overlap BPD group. 2Y BPD only: too few metabolites for analysis.",
         size=14, italic=True, color=DGRAY, first=True)

    # Conclusions (right of pathway table)
    ccy = PTH_Y + Cm(0.3)
    ccy = heading(slide, "Conclusions", con_x, ccy, con_w, size=26)

    conclusions = [
        ("•  Metabolic signatures of prematurity and BPD persist from 6M to 2Y "
         "and may intensify over time."),
        ("•  Gut microbiome-derived metabolites (N-PAG, indoxyl sulfate) "
         "consistently elevated in extremely preterm and M+S BPD infants."),
        ("•  Complete GA–BPD metabolic overlap at 6M followed by divergence "
         "at 2Y underscores evolving, BPD-specific dysregulation."),
        ("•  Persistently enriched amino acid pathways implicate disrupted nitrogen "
         "metabolism beyond the acute BPD phase, with implications for biomarker "
         "identification and therapeutic intervention."),
    ]
    for line in conclusions:
        nchars = len(line)
        est = max(1, nchars // 52 + 1)
        tb, tf = txbox(slide, con_x, ccy, con_w, Cm(est * 0.85 + 0.2))
        para(tf, line, size=18, color=DGRAY, first=True)
        ccy += Cm(est * 0.78 + 0.35)

    # ──────────────────────────────────────────────────────────────────────────
    # FOOTER
    # ──────────────────────────────────────────────────────────────────────────
    rect(slide, 0, FTR_Y, A0_W, FTR_H, fill=NAVY)

    # References (left ~55cm)
    ref_x, ref_y = MAR, FTR_Y + Cm(0.8)
    tb, tf = txbox(slide, ref_x, ref_y, Cm(55), Cm(1.0))
    para(tf, "REFERENCES", size=22, bold=True, color=WHITE, first=True)
    rect(slide, ref_x, ref_y + Cm(0.9), Cm(55), Cm(0.1), fill=TEAL)

    refs = (
        "1. Jobe AH, Bancalari E. Bronchopulmonary Dysplasia. AJRCCM 2001;163(7):1723-9.\n"
        "2. Higgins RD, et al. Defining Bronchopulmonary Dysplasia (NICHD Workshop). Pediatrics 2018;142(6).\n"
        "3. Chong J, et al. MetaboAnalyst 5.0. Nucleic Acids Res 2021;49(W1):W388–W396.\n"
        "4. Wishart DS, et al. HMDB 5.0. Nucleic Acids Res 2022;50(D1):D622–D631."
    )
    tb, tf = txbox(slide, ref_x, ref_y + Cm(1.1), Cm(55), Cm(10))
    para(tf, refs, size=17, color=WHITE, first=True)

    # Contact (right ~25cm)
    con2_x = Cm(60)
    tb, tf = txbox(slide, con2_x, ref_y, Cm(22), Cm(1.0))
    para(tf, "CONTACT", size=22, bold=True, color=WHITE, first=True)
    rect(slide, con2_x, ref_y + Cm(0.9), Cm(22), Cm(0.1), fill=TEAL)

    contact = (
        "Chih-Yung Chiu, MD PhD\n"
        "Multi-omics Core Laboratory\n"
        "Chang Gung Memorial Hospital at Linkou\n"
        "Taoyuan, Taiwan\n"
        "pedchestIC@gmail.com"
    )
    tb, tf = txbox(slide, con2_x, ref_y + Cm(1.1), Cm(22), Cm(10))
    para(tf, contact, size=18, color=WHITE, first=True)

    # ──────────────────────────────────────────────────────────────────────────
    prs.save(OUTPUT)
    print(f"Saved: {OUTPUT}")


if __name__ == "__main__":
    main()
