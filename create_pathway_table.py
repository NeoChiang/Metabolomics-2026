#!/usr/bin/env python3
"""
Generate Table 4 Word document from pathway analysis results,
matching the style of Table 4 _Pathways_20250401_R.docx:
  - 7 columns: Metabolites | Pathway Name | Total | Hits | Raw P | FDR | Function
  - Section-divider rows (bold, top border only)
  - Header row (top + bottom border)
  - Last data row (bottom border)
  - No vertical borders
  - Note rows (italic, merged across all columns) for empty sections
"""

from docx import Document
from docx.shared import Pt, Cm
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

OUTPUT = "pathway_results/Table4_Pathways_preterm_infants.docx"
FS = 9   # table body font size (pt)
FS_CAP = 10  # caption / footnote font size

# ── Helpers ───────────────────────────────────────────────────────────────────

def fmt_p(v):
    """Format p-value / FDR for display."""
    if v < 0.001:
        return "<0.001"
    if v < 0.01:
        return f"{v:.4f}"
    return f"{v:.3f}"


def add_run(para, text, bold=False, italic=False, size=None):
    r = para.add_run(text)
    r.bold = bold
    r.italic = italic
    if size:
        r.font.size = Pt(size)
    return r


def _border_xml(name, style):
    el = OxmlElement(f"w:{name}")
    if style == "single":
        el.set(qn("w:val"), "single")
        el.set(qn("w:sz"), "4")
        el.set(qn("w:space"), "0")
        el.set(qn("w:color"), "000000")
    else:
        el.set(qn("w:val"), "none")
    return el


def set_cell_borders(cell, top="none", bottom="none"):
    tc = cell._tc
    tcPr = tc.get_or_add_tcPr()
    for old in tcPr.findall(qn("w:tcBorders")):
        tcPr.remove(old)
    b = OxmlElement("w:tcBorders")
    for name, style in [("top", top), ("bottom", bottom),
                         ("left", "none"), ("right", "none"),
                         ("insideH", "none"), ("insideV", "none")]:
        b.append(_border_xml(name, style))
    tcPr.append(b)


def set_row_borders(row, top=False, bottom=False):
    for cell in row.cells:
        set_cell_borders(cell,
                         top="single" if top else "none",
                         bottom="single" if bottom else "none")


def clear_table_borders(table):
    """Remove all automatic table-level borders."""
    tbl = table._tbl
    tblPr = tbl.tblPr
    for old in tblPr.findall(qn("w:tblBorders")):
        tblPr.remove(old)
    tblBorders = OxmlElement("w:tblBorders")
    for edge in ["top", "left", "bottom", "right", "insideH", "insideV"]:
        el = OxmlElement(f"w:{edge}")
        el.set(qn("w:val"), "none")
        tblBorders.append(el)
    tblPr.append(tblBorders)


def fill_cell(cell, text, bold=False, italic=False, size=FS):
    p = cell.paragraphs[0]
    add_run(p, text, bold=bold, italic=italic, size=size)


# ── Table data ────────────────────────────────────────────────────────────────

HEADERS = ["Metabolites", "Pathway Name", "Total", "Hits", "Raw P", "FDR", "Function"]

# Column widths sum to 16.0 cm (fits A4 with 2.54 cm margins each side)
COL_WIDTHS = [Cm(3.5), Cm(4.0), Cm(1.0), Cm(1.0), Cm(1.7), Cm(1.7), Cm(3.1)]

SECTIONS = [
    {
        "label": "6M GA only",
        "note": ("All GA-significant metabolites at 6M were also significant "
                 "in BPD group"),
        "rows": [],
    },
    {
        "label": "6M BPD only",
        "note": None,
        "rows": [
            ("N-Acetylaspartic acid, Glutamine, Citric acid, Succinic acid",
             "Alanine, aspartate and glutamate metabolism",
             28, 4, 9.1045e-06, 0.00072836, "Amino acid metabolism"),
            ("Acetic acid, Glutamine, Citric acid",
             "Glyoxylate and dicarboxylate metabolism",
             32, 3, 0.00057151, 0.022861, "Carbohydrate metabolism"),
            ("Citric acid, Succinic acid",
             "Citrate cycle (TCA cycle)",
             20, 2, 0.0051224, 0.1366, "Carbohydrate metabolism"),
            ("Glutamine",
             "Nitrogen metabolism",
             6, 1, 0.033496, 0.66991, "Energy metabolism"),
        ],
    },
    {
        "label": "6M Both GA+BPD",
        "note": None,
        "rows": [
            ("3-Methyl-2-oxopentanoic acid, Valine",
             "Valine, leucine and isoleucine biosynthesis",
             8, 2, 0.00078202, 0.062561, "Amino acid metabolism"),
            ("Pantothenic acid, Valine",
             "Pantothenate and CoA biosynthesis",
             20, 2, 0.0051224, 0.2049, "Metabolism of cofactors and vitamins"),
            ("N,N-Dimethylglycine, Creatine",
             "Glycine, serine and threonine metabolism",
             33, 2, 0.0137, 0.29924, "Amino acid metabolism"),
            ("3-Methyl-2-oxopentanoic acid, Valine",
             "Valine, leucine and isoleucine degradation",
             40, 2, 0.019825, 0.29924, "Amino acid metabolism"),
            ("4-Hydroxyphenylacetic acid, Tyrosine",
             "Tyrosine metabolism",
             42, 2, 0.021755, 0.29924, "Amino acid metabolism"),
            ("Tyrosine",
             "Phenylalanine, tyrosine and tryptophan biosynthesis",
             4, 1, 0.022443, 0.29924, "Amino acid metabolism"),
            ("Tyrosine",
             "Phenylalanine metabolism",
             8, 1, 0.044437, 0.50785, "Amino acid metabolism"),
        ],
    },
    {
        "label": "2Y GA only",
        "note": None,
        "rows": [
            ("Tyrosine, 4-Hydroxyphenylacetate",
             "Tyrosine metabolism",
             42, 2, 0.017206, 0.79173, "Amino acid metabolism"),
            ("Tyrosine",
             "Phenylalanine, tyrosine and tryptophan biosynthesis",
             4, 1, 0.019968, 0.79173, "Amino acid metabolism"),
            ("Tyrosine",
             "Phenylalanine metabolism",
             8, 1, 0.039587, 0.79173, "Amino acid metabolism"),
            ("Taurine",
             "Taurine and hypotaurine metabolism",
             8, 1, 0.039587, 0.79173, "Metabolism of other amino acids"),
        ],
    },
    {
        "label": "2Y BPD only",
        "note": ("Too few metabolites to perform meaningful enrichment analysis "
                 "(2-Aminobutyric acid, Citric acid, Propylene glycol)"),
        "rows": [],
    },
    {
        "label": "2Y Both GA+BPD",
        "note": None,
        "rows": [
            ("Glutamine, Succinic acid",
             "Alanine, aspartate and glutamate metabolism",
             28, 2, 0.0042853, 0.22341, "Amino acid metabolism"),
            ("Glycine, Glutamine",
             "Glyoxylate and dicarboxylate metabolism",
             32, 2, 0.0055852, 0.22341, "Carbohydrate metabolism"),
            ("Glutamine",
             "Nitrogen metabolism",
             6, 1, 0.022436, 0.57798, "Energy metabolism"),
            ("3-Methyl-2-oxopentanoic acid",
             "Valine, leucine and isoleucine biosynthesis",
             8, 1, 0.029821, 0.57798, "Amino acid metabolism"),
        ],
    },
]

# ── Build document ────────────────────────────────────────────────────────────

def main():
    doc = Document()

    # Page margins
    for sec in doc.sections:
        sec.left_margin = Cm(2.54)
        sec.right_margin = Cm(2.54)
        sec.top_margin = Cm(2.54)
        sec.bottom_margin = Cm(2.54)

    # ── Caption ───────────────────────────────────────────────────────────────
    cap = doc.add_paragraph()
    add_run(cap, "Table 4.", bold=True, size=FS_CAP)
    add_run(cap, (
        " Metabolic pathway enrichment analysis of urinary NMR metabolites "
        "significantly expressed in preterm infants stratified by gestational "
        "age (GA) and bronchopulmonary dysplasia (BPD) severity at 6 months "
        "(6M) and 2 years (2Y) corrected age."
    ), size=FS_CAP)

    # ── Count rows ────────────────────────────────────────────────────────────
    n_rows = 1  # header
    for s in SECTIONS:
        n_rows += 1                              # section-divider row
        if s["note"] and not s["rows"]:
            n_rows += 1                          # note row
        n_rows += len(s["rows"])                 # data rows

    # ── Create table ──────────────────────────────────────────────────────────
    table = doc.add_table(rows=n_rows, cols=7)
    clear_table_borders(table)

    # Set column widths for every row
    for row in table.rows:
        for i, w in enumerate(COL_WIDTHS):
            row.cells[i].width = w

    # ── Header row ────────────────────────────────────────────────────────────
    hdr = table.rows[0]
    for i, h in enumerate(HEADERS):
        p = hdr.cells[i].paragraphs[0]
        if h == "Raw P":
            add_run(p, "Raw ", bold=True, size=FS)
            add_run(p, "P", bold=True, italic=True, size=FS)
        else:
            add_run(p, h, bold=True, size=FS)
    set_row_borders(hdr, top=True, bottom=True)

    # ── Section + data rows ───────────────────────────────────────────────────
    ri = 1
    for si, sec in enumerate(SECTIONS):
        is_last_sec = si == len(SECTIONS) - 1

        # Section-divider row
        sec_row = table.rows[ri]
        fill_cell(sec_row.cells[0], sec["label"], bold=True)
        set_row_borders(sec_row, top=True, bottom=False)
        ri += 1

        # Note row for empty sections (merged across all 7 columns)
        if sec["note"] and not sec["rows"]:
            note_row = table.rows[ri]
            merged = note_row.cells[0].merge(note_row.cells[6])
            fill_cell(merged, sec["note"], italic=True)
            set_cell_borders(merged, top="none", bottom="none")
            ri += 1

        # Data rows
        for di, row_data in enumerate(sec["rows"]):
            met, pwy, total, hits, rawp, fdr, func = row_data
            is_last_row = is_last_sec and di == len(sec["rows"]) - 1
            dr = table.rows[ri]
            for i, v in enumerate([met, pwy, str(total), str(hits),
                                    fmt_p(rawp), fmt_p(fdr), func]):
                fill_cell(dr.cells[i], v)
            set_row_borders(dr, top=False, bottom=is_last_row)
            ri += 1

    # ── Footnote ──────────────────────────────────────────────────────────────
    fn = doc.add_paragraph()
    add_run(fn, "Total, total number of compounds in the pathway; "
                "Hits, matched number from the uploaded data; Raw ", size=FS_CAP)
    add_run(fn, "P", italic=True, size=FS_CAP)
    add_run(fn, ", original ", size=FS_CAP)
    add_run(fn, "P", italic=True, size=FS_CAP)
    add_run(fn, (" value from the enrichment analysis; FDR, false discovery rate "
                 "(Benjamini-Hochberg correction)."), size=FS_CAP)

    doc.save(OUTPUT)
    print(f"Saved: {OUTPUT}")


if __name__ == "__main__":
    main()
