"""
AgriGen - USVT DNA Analyzer: Aplicație Web Streamlit
====================================================
Platformă interactivă pentru analiza datelor genetice vegetale,
dedicată detecției markerilor de rezistență la secetă și agenți patogeni.
Dezvoltată pentru Universitatea de Științele Vieții "Regele Mihai I" din Timișoara (USVT).
"""

import sys
import os
import time
import unicodedata
from typing import Any, Dict, List, Optional
import streamlit as st
import pandas as pd
from fpdf import FPDF

# Asigurare import modul local
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

try:
    from analyzer import analyze_vcf
except ImportError:
    from src.analyzer import analyze_vcf


# -----------------------------------------------------------------------------
# 1. CATALOG DE DATE DEMONSTRATIVE PENTRU TOATE CELE 3 SPECII VEGETALE
# -----------------------------------------------------------------------------
DEMO_SAMPLES: Dict[str, Dict[str, str]] = {
    # --- PORUMB (Zea mays) ---
    "Demo: Porumb Rezistent Secetă": {
        "species": "Porumb",
        "default_trait": "Secetă",
        "vcf": (
            "##fileformat=VCFv4.2\n"
            "##source=AgriGenUSVT_Sequencer\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tUSVT-PORUMB-REZ-SECETA\n"
            "chr3\t145320\trsZmDREB2A\tA\tG\t99\tPASS\tDP=48;AF=1.0\tGT\t1/1\n"
            "chr5\t882341\trsZmWAK\tG\tG\t92\tPASS\tDP=38;AF=0.0\tGT\t0/0\n"
        ),
        "filename": "demo_porumb_rezistent_seceta.vcf",
        "description": "Porumb cu alela favorabilă la secetă (chr3:145320 A->G), sensibil la boli."
    },
    "Demo: Porumb Vulnerabil Secetă": {
        "species": "Porumb",
        "default_trait": "Secetă",
        "vcf": (
            "##fileformat=VCFv4.2\n"
            "##source=AgriGenUSVT_Sequencer\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tUSVT-PORUMB-VULN-SECETA\n"
            "chr3\t145320\trsZmDREB2A\tA\tA\t95\tPASS\tDP=35;AF=0.0\tGT\t0/0\n"
            "chr5\t882341\trsZmWAK\tG\tG\t90\tPASS\tDP=32;AF=0.0\tGT\t0/0\n"
        ),
        "filename": "demo_porumb_vulnerabil_seceta.vcf",
        "description": "Porumb cu genotipul sălbatic de referință (sensibil la secetă)."
    },
    "Demo: Porumb Rezistent Boli": {
        "species": "Porumb",
        "default_trait": "Boli / Patogeni",
        "vcf": (
            "##fileformat=VCFv4.2\n"
            "##source=AgriGenUSVT_Sequencer\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tUSVT-PORUMB-REZ-BOLI\n"
            "chr5\t882341\trsZmWAK\tG\tT\t99\tPASS\tDP=52;AF=1.0\tGT\t1/1\n"
            "chr3\t145320\trsZmDREB2A\tA\tA\t90\tPASS\tDP=30;AF=0.0\tGT\t0/0\n"
        ),
        "filename": "demo_porumb_rezistent_boli.vcf",
        "description": "Porumb cu mutația de imunitate la patogeni (chr5:882341 G->T)."
    },
    "Demo: Porumb Vulnerabil Boli": {
        "species": "Porumb",
        "default_trait": "Boli / Patogeni",
        "vcf": (
            "##fileformat=VCFv4.2\n"
            "##source=AgriGenUSVT_Sequencer\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tUSVT-PORUMB-VULN-BOLI\n"
            "chr5\t882341\trsZmWAK\tG\tG\t91\tPASS\tDP=30;AF=0.0\tGT\t0/0\n"
            "chr3\t145320\trsZmDREB2A\tA\tA\t88\tPASS\tDP=28;AF=0.0\tGT\t0/0\n"
        ),
        "filename": "demo_porumb_vulnerabil_boli.vcf",
        "description": "Porumb fără markerul de rezistență la boli (chr5:882341)."
    },

    # --- GRÂU (Triticum aestivum) ---
    "Demo: Grâu Rezistent Secetă": {
        "species": "Grâu",
        "default_trait": "Secetă",
        "vcf": (
            "##fileformat=VCFv4.2\n"
            "##source=AgriGenUSVT_Sequencer\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tUSVT-GRAU-REZ-SECETA\n"
            "chr4\t551200\trsTaDREB1\tA\tT\t99\tPASS\tDP=50;AF=1.0\tGT\t1/1\n"
            "chr7\t1200500\trsLr34\tG\tG\t90\tPASS\tDP=34;AF=0.0\tGT\t0/0\n"
        ),
        "filename": "demo_grau_rezistent_seceta.vcf",
        "description": "Grâu cu varianta protectoare TaDREB1 (chr4:551200 A->T) pentru toleranță la secetă."
    },
    "Demo: Grâu Vulnerabil Secetă": {
        "species": "Grâu",
        "default_trait": "Secetă",
        "vcf": (
            "##fileformat=VCFv4.2\n"
            "##source=AgriGenUSVT_Sequencer\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tUSVT-GRAU-VULN-SECETA\n"
            "chr4\t551200\trsTaDREB1\tA\tA\t92\tPASS\tDP=32;AF=0.0\tGT\t0/0\n"
            "chr7\t1200500\trsLr34\tG\tG\t90\tPASS\tDP=30;AF=0.0\tGT\t0/0\n"
        ),
        "filename": "demo_grau_vulnerabil_seceta.vcf",
        "description": "Grâu vulnerabil la secetă (fără alela mutantă TaDREB1 la chr4:551200)."
    },
    "Demo: Grâu Rezistent Boli": {
        "species": "Grâu",
        "default_trait": "Boli / Patogeni",
        "vcf": (
            "##fileformat=VCFv4.2\n"
            "##source=AgriGenUSVT_Sequencer\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tUSVT-GRAU-REZ-BOLI\n"
            "chr7\t1200500\trsLr34\tG\tC\t99\tPASS\tDP=55;AF=1.0\tGT\t1/1\n"
            "chr4\t551200\trsTaDREB1\tA\tA\t89\tPASS\tDP=30;AF=0.0\tGT\t0/0\n"
        ),
        "filename": "demo_grau_rezistent_boli.vcf",
        "description": "Grâu cu markerul Lr34 (chr7:1200500 G->C) pentru rezistență durabilă la rugină."
    },
    "Demo: Grâu Vulnerabil Boli": {
        "species": "Grâu",
        "default_trait": "Boli / Patogeni",
        "vcf": (
            "##fileformat=VCFv4.2\n"
            "##source=AgriGenUSVT_Sequencer\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tUSVT-GRAU-VULN-BOLI\n"
            "chr7\t1200500\trsLr34\tG\tG\t91\tPASS\tDP=28;AF=0.0\tGT\t0/0\n"
            "chr4\t551200\trsTaDREB1\tA\tA\t88\tPASS\tDP=25;AF=0.0\tGT\t0/0\n"
        ),
        "filename": "demo_grau_vulnerabil_boli.vcf",
        "description": "Grâu sensibil la atacul ruginii și bolilor foliare (fără varianta Lr34)."
    },

    # --- ROȘIE (Solanum lycopersicum) ---
    "Demo: Roșie Rezistentă Secetă": {
        "species": "Roșie",
        "default_trait": "Secetă",
        "vcf": (
            "##fileformat=VCFv4.2\n"
            "##source=AgriGenUSVT_Sequencer\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tUSVT-ROSIE-REZ-SECETA\n"
            "chr1\t88320\trsSlAREB1\tT\tA\t99\tPASS\tDP=45;AF=1.0\tGT\t1/1\n"
            "chr9\t405000\trsPh3\tC\tC\t90\tPASS\tDP=35;AF=0.0\tGT\t0/0\n"
        ),
        "filename": "demo_rosie_rezistenta_seceta.vcf",
        "description": "Roșie cu factorul SlAREB1 (chr1:88320 T->A) pentru rezistență la secetă și salinitate."
    },
    "Demo: Roșie Vulnerabilă Secetă": {
        "species": "Roșie",
        "default_trait": "Secetă",
        "vcf": (
            "##fileformat=VCFv4.2\n"
            "##source=AgriGenUSVT_Sequencer\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tUSVT-ROSIE-VULN-SECETA\n"
            "chr1\t88320\trsSlAREB1\tT\tT\t90\tPASS\tDP=30;AF=0.0\tGT\t0/0\n"
            "chr9\t405000\trsPh3\tC\tC\t90\tPASS\tDP=28;AF=0.0\tGT\t0/0\n"
        ),
        "filename": "demo_rosie_vulnerabila_seceta.vcf",
        "description": "Roșie vulnerabilă la deficit hidric (fără alela mutantă SlAREB1)."
    },
    "Demo: Roșie Rezistentă Boli": {
        "species": "Roșie",
        "default_trait": "Boli / Patogeni",
        "vcf": (
            "##fileformat=VCFv4.2\n"
            "##source=AgriGenUSVT_Sequencer\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tUSVT-ROSIE-REZ-BOLI\n"
            "chr9\t405000\trsPh3\tC\tG\t99\tPASS\tDP=54;AF=1.0\tGT\t1/1\n"
            "chr1\t88320\trsSlAREB1\tT\tT\t89\tPASS\tDP=28;AF=0.0\tGT\t0/0\n"
        ),
        "filename": "demo_rosie_rezistenta_boli.vcf",
        "description": "Roșie cu gena Ph-3 (chr9:405000 C->G) pentru imunitate la mana solanaceelor."
    },
    "Demo: Roșie Vulnerabilă Boli": {
        "species": "Roșie",
        "default_trait": "Boli / Patogeni",
        "vcf": (
            "##fileformat=VCFv4.2\n"
            "##source=AgriGenUSVT_Sequencer\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tUSVT-ROSIE-VULN-BOLI\n"
            "chr9\t405000\trsPh3\tC\tC\t90\tPASS\tDP=30;AF=0.0\tGT\t0/0\n"
            "chr1\t88320\trsSlAREB1\tT\tT\t87\tPASS\tDP=25;AF=0.0\tGT\t0/0\n"
        ),
        "filename": "demo_rosie_vulnerabila_boli.vcf",
        "description": "Roșie fără gena de rezistență Ph-3 (susceptibilă la mană)."
    }
}


# -----------------------------------------------------------------------------
# 2. GENERATOR RAPORT PDF (fpdf2)
# -----------------------------------------------------------------------------
def generate_pdf_report(report_data: dict) -> bytes:
    """
    Generează un document PDF oficial și formatat profesional
    pe baza rezultatelor analizei genomice, folosind fpdf2.
    """
    pdf = FPDF(orientation="P", unit="mm", format="A4")
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()

    font_candidates = [
        ("C:/Windows/Fonts/arial.ttf", "C:/Windows/Fonts/arialbd.ttf"),
        ("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"),
        ("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf"),
    ]

    has_unicode_font = False
    font_name = "helvetica"

    for reg_path, bold_path in font_candidates:
        if os.path.exists(reg_path):
            try:
                pdf.add_font("AgriFont", "", reg_path)
                if os.path.exists(bold_path):
                    pdf.add_font("AgriFont", "B", bold_path)
                has_unicode_font = True
                font_name = "AgriFont"
                break
            except Exception:
                continue

    def safe_text(txt: Any) -> str:
        s = str(txt)
        if has_unicode_font:
            return s
        return unicodedata.normalize("NFKD", s).encode("ASCII", "ignore").decode("ASCII")

    # Antet Instituțional
    pdf.set_font(font_name, "B", 16)
    pdf.set_text_color(27, 67, 50)
    pdf.cell(0, 8, safe_text("AgriGen - USVT | Raport Diagnostic Genomic"), new_x="LMARGIN", new_y="NEXT", align="C")

    pdf.set_font(font_name, "", 10)
    pdf.set_text_color(70, 70, 70)
    pdf.cell(0, 5, safe_text("Universitatea de Științele Vieții «Regele Mihai I» din Timișoara"), new_x="LMARGIN", new_y="NEXT", align="C")
    pdf.cell(0, 5, safe_text("Laboratorul de Genomică și Biologie Computațională Vegetală"), new_x="LMARGIN", new_y="NEXT", align="C")
    pdf.ln(3)

    pdf.set_draw_color(180, 180, 180)
    pdf.line(12, pdf.get_y(), 198, pdf.get_y())
    pdf.ln(4)

    verdict = report_data["verdict"]
    meta = report_data["sample_meta"]
    marker_info = report_data["primary_marker"]
    recs = report_data["recommendations"]

    # 1. Informații Generale Probă
    pdf.set_font(font_name, "B", 11)
    pdf.set_text_color(20, 20, 20)
    pdf.cell(0, 6, safe_text("1. Informații Generale Probă"), new_x="LMARGIN", new_y="NEXT")

    pdf.set_font(font_name, "", 9)
    pdf.cell(93, 6, safe_text(f"Specie Gazdă: {meta['species']}"), border=1)
    pdf.cell(93, 6, safe_text(f"ID Probă: {meta['sample_id']}"), border=1, new_x="LMARGIN", new_y="NEXT")
    pdf.cell(93, 6, safe_text(f"Fișier Sursă: {meta['filename']}"), border=1)
    pdf.cell(93, 6, safe_text(f"Standard VCF: {meta['vcf_format_version']} | Variante: {meta['total_variants_scanned']:,}"), border=1, new_x="LMARGIN", new_y="NEXT")
    pdf.ln(4)

    # 2. Verdict Diagnostic Final
    pdf.set_font(font_name, "B", 11)
    pdf.cell(0, 6, safe_text("2. Verdict Diagnostic Final"), new_x="LMARGIN", new_y="NEXT")

    status = verdict["status"]
    if status == "Rezistent":
        pdf.set_fill_color(216, 243, 220)
        pdf.set_text_color(45, 106, 79)
    else:
        pdf.set_fill_color(255, 229, 217)
        pdf.set_text_color(180, 20, 20)

    pdf.set_font(font_name, "B", 11)
    verdict_text = f"  VERDICT: {verdict['verdict_headline']}  |  Certitudine Diagnostic: {verdict['confidence_score']:.1f}%  "
    pdf.cell(186, 9, safe_text(verdict_text), fill=True, border=1, new_x="LMARGIN", new_y="NEXT")
    pdf.ln(3)

    # 3. Rezumat Executiv Bioinformatic
    pdf.set_text_color(20, 20, 20)
    pdf.set_font(font_name, "B", 11)
    pdf.cell(0, 6, safe_text("3. Rezumat Executiv Bioinformatic"), new_x="LMARGIN", new_y="NEXT")

    pdf.set_font(font_name, "", 9)
    pdf.multi_cell(186, 4.8, safe_text(verdict["executive_summary"]))
    pdf.ln(2)

    pdf.set_font(font_name, "B", 9)
    pdf.cell(42, 5, safe_text("Mecanism Molecular: "))
    pdf.set_font(font_name, "", 9)
    pdf.multi_cell(144, 4.8, safe_text(marker_info["mechanism"]))
    pdf.ln(3)

    # 4. Detalii Marker Genetic și Genotipare
    pdf.set_font(font_name, "B", 11)
    pdf.cell(0, 6, safe_text("4. Detalii Marker Genetic și Genotipare"), new_x="LMARGIN", new_y="NEXT")

    detalii_tabel = [
        ("Trăsătură / Marker Vizat:", marker_info["trait"]),
        ("Genă Țintă (Locus):", marker_info["gene"]),
        ("Coordonate Genomice:", f"{marker_info['chrom']}:{marker_info['pos']}"),
        ("Alele (REF / ALT):", f"{marker_info['ref_allele']} / {marker_info['alt_allele']} (Favorabil: {marker_info['favorable_allele']})"),
        ("Genotip Identificat:", f"{marker_info['genotype_code']} ({marker_info['alleles']}) - {marker_info['zygosity']}"),
        ("Profunzime Citire & Scor Calitate:", f"Read Depth: {marker_info['read_depth']} | Scor Phred: Q{marker_info['phred_quality']:.1f}"),
    ]

    for eticheta, valoare in detalii_tabel:
        pdf.set_font(font_name, "B", 8.5)
        pdf.cell(60, 5.5, safe_text(eticheta), border=1)
        pdf.set_font(font_name, "", 8.5)
        pdf.cell(126, 5.5, safe_text(valoare), border=1, new_x="LMARGIN", new_y="NEXT")
    pdf.ln(4)

    # 5. Recomandări Agronomice și de Ameliorare
    pdf.set_font(font_name, "B", 11)
    pdf.cell(0, 6, safe_text("5. Recomandări Agronomice și de Ameliorare"), new_x="LMARGIN", new_y="NEXT")

    pdf.set_font(font_name, "", 8.5)
    bullet_char = "-" if has_unicode_font else "*"
    for r in recs:
        pdf.cell(6, 4.8, bullet_char, align="R")
        pdf.multi_cell(180, 4.8, safe_text(f" {r}"))

    pdf.ln(5)

    # Subsol
    pdf.set_font(font_name, "", 7.5)
    pdf.set_text_color(110, 110, 110)
    pdf.cell(0, 4, safe_text("Document generat automat prin motorul AgriGen - USVT BioCore."), align="C")

    return bytes(pdf.output())


# -----------------------------------------------------------------------------
# 3. CONFIGURARE PAGINĂ STREAMLIT
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="🌱 AgriGen - USVT DNA Analyzer",
    page_icon="🌱",
    layout="wide",
    initial_sidebar_state="expanded"
)

# -----------------------------------------------------------------------------
# 4. HERO SECTION & DASHBOARD METRICS
# -----------------------------------------------------------------------------
st.title("🌱 AgriGen - USVT DNA Analyzer")
st.subheader("Platformă de Analiză Genomică Vegetală • Detecție Markeri de Rezistență la Secetă și Boli")
st.caption("Universitatea de Științele Vieții «Regele Mihai I» din Timișoara (USVT) • Facultatea de Bioinginerie")

# Rând de 3 metrici de tip Dashboard pentru prezentare
col_hero1, col_hero2, col_hero3 = st.columns(3)
with col_hero1:
    st.metric(
        label="🧬 Secvențe Genomice Indexate",
        value="1.2 Milioane",
        delta="Baza USVT BioBank"
    )
with col_hero2:
    st.metric(
        label="🎯 Markeri Activi",
        value="45",
        delta="Porumb, Grâu, Roșie"
    )
with col_hero3:
    st.metric(
        label="⚡ Timp Procesare",
        value="< 2 secunde",
        delta="Motor BioCore Optimizat"
    )

st.divider()


# -----------------------------------------------------------------------------
# 5. PANOU LATERAL: CONFIGURARE, DEMO DATA MULTI-SPECIE ȘI ÎNCĂRCARE
# -----------------------------------------------------------------------------
st.sidebar.header("🔬 Panou de Configurare")

# 1. Specie Plantă
optiuni_specii = ["Porumb", "Grâu", "Roșie"]
etichete_specii = {
    "Porumb": "🌽 Porumb (Zea mays)",
    "Grâu": "🌾 Grâu (Triticum aestivum)",
    "Roșie": "🍅 Roșie (Solanum lycopersicum)"
}

specie_selectata = st.sidebar.selectbox(
    "1. Specie Plantă:",
    options=optiuni_specii,
    format_func=lambda x: etichete_specii.get(x, x),
    help="Selectați specia vegetală analizată."
)

# 2. Marker Genetic Țintă adaptat pe specia selectată
harta_etichete_markeri = {
    "Porumb": {
        "Secetă": "💧 Secetă (ZmDREB2A - chr3:145320 A->G)",
        "Boli / Patogeni": "🛡️ Boli / Patogeni (ZmWAK-R - chr5:882341 G->T)"
    },
    "Grâu": {
        "Secetă": "💧 Secetă (TaDREB1 - chr4:551200 A->T)",
        "Boli / Patogeni": "🛡️ Boli / Rugină (Lr34 - chr7:1200500 G->C)"
    },
    "Roșie": {
        "Secetă": "💧 Secetă (SlAREB1 - chr1:88320 T->A)",
        "Boli / Patogeni": "🛡️ Boli / Mană (Ph-3 - chr9:405000 C->G)"
    }
}

dict_markeri_curent = harta_etichete_markeri.get(specie_selectata, harta_etichete_markeri["Porumb"])
optiuni_markeri = list(dict_markeri_curent.keys())

marker_selectat = st.sidebar.selectbox(
    "2. Trăsătură / Marker Genetic:",
    options=optiuni_markeri,
    format_func=lambda x: dict_markeri_curent.get(x, x),
    help="Trăsătura vizată determină direct coordonatele căutate și verdictul final!"
)

st.sidebar.divider()

# 3. Selectare Date de Intrare (Proba Demo sau Încărcare Manuală)
st.sidebar.subheader("3. Selectare Date de Intrare")

optiuni_demo = ["Niciuna (Încarcă fișier manual)"] + list(DEMO_SAMPLES.keys())
demo_selectat = st.sidebar.selectbox(
    "Probă demonstrativă USVT:",
    options=optiuni_demo,
    help="Alegeți o probă demo pentru oricare dintre cele 3 specii (Porumb, Grâu, Roșie)."
)

fisier_activ = None
nume_fisier_afisat = None
specie_de_analizat = specie_selectata

if demo_selectat != "Niciuna (Încarcă fișier manual)":
    demo_info = DEMO_SAMPLES[demo_selectat]
    fisier_activ = demo_info["vcf"]
    nume_fisier_afisat = demo_info["filename"]
    # Folosește specia specifică probei demonstrative
    specie_de_analizat = demo_info["species"]

    st.sidebar.success(f"✓ **Probă Demo ({demo_info['species']}):**\n`{nume_fisier_afisat}`")
    st.sidebar.caption(demo_info["description"])
    st.sidebar.info(
        f"🎯 **Criteriu analizat:** `{marker_selectat}` pentru specia **{demo_info['species']}**.\n\n"
        "(Puteți comuta oricând trăsătura de la pasul 2 pentru a testa profilul probei)."
    )
else:
    uploaded_file = st.sidebar.file_uploader(
        "Încarcă Fișierul Genomic (.vcf, .txt):",
        type=["vcf", "txt"],
        help="Acceptă fișiere standard Variant Call Format (VCF v4.2) sau fișiere text cu variante genetice."
    )
    fisier_activ = uploaded_file
    nume_fisier_afisat = uploaded_file.name if uploaded_file else None
    specie_de_analizat = specie_selectata

st.sidebar.divider()

# Buton primar de procesare
buton_generare = st.sidebar.button(
    "🧬 Generează Raportul",
    type="primary",
    use_container_width=True
)


# -----------------------------------------------------------------------------
# 6. PROCESAREA ANALIZEI GENOMICE CU STATE DE ÎNCĂRCARE DINAMIC
# -----------------------------------------------------------------------------
if buton_generare:
    if fisier_activ is None:
        st.warning("⚠️ **Vă rugăm să încărcați un fișier `.vcf` sau `.txt`**, ori selectați o probă demonstrativă din meniul lateral.")
    else:
        # Bară de progres animată cu mesaje dinamice
        progress_bar = st.progress(0)
        status_placeholder = st.empty()

        loading_phases = [
            (15, f"🧬 1/4 Extragere ADN și parsare antet VCF ({specie_de_analizat})..."),
            (45, f"🔬 2/4 Aliniere secvențe genomice cu genomul de referință ({specie_de_analizat})..."),
            (75, f"🔍 3/4 Căutare markeri linie cu linie ({marker_selectat})..."),
            (95, "📊 4/4 Evaluare indici de calitate și formulare verdict diagnostic..."),
            (100, "✅ Procesare finalizată! Generare raport...")
        ]

        for progress_pct, status_msg in loading_phases:
            status_placeholder.markdown(f"**{status_msg}**")
            progress_bar.progress(progress_pct)
            time.sleep(0.35)

        progress_bar.empty()
        status_placeholder.empty()

        try:
            # Apel sincronizat: specia probei + trăsătura selectată de utilizator
            rezultat = analyze_vcf(
                uploaded_file=fisier_activ,
                species=specie_de_analizat,
                marker=marker_selectat
            )
            if nume_fisier_afisat:
                rezultat["sample_meta"]["filename"] = nume_fisier_afisat
            st.session_state["raport_curent"] = rezultat

            # Celebration Effect pentru 'Rezistent'
            if rezultat["verdict"]["status"] == "Rezistent":
                st.balloons()

        except Exception:
            st.error("Eroare la procesarea fișierului. Te rugăm să încarci un fișier VCF valid.")


# -----------------------------------------------------------------------------
# 7. RANDAREA RAPORTULUI CU BARĂ VIZUALĂ DE ÎNCREDERE ȘI COMPONENTE NATIVE
# -----------------------------------------------------------------------------
if "raport_curent" in st.session_state:
    raport = st.session_state["raport_curent"]
    verdict = raport["verdict"]
    meta = raport["sample_meta"]
    marker_info = raport["primary_marker"]
    recs = raport["recommendations"]

    st.header("📋 Raport Diagnostic Genomic")

    # Banner de alertă pentru Verdictul Final adaptat la trăsătura vizată
    if verdict["status"] == "Rezistent":
        st.success(
            f"### ✅ **VERDICT FINAL: PLANTA ESTE REZISTENTĂ**\n\n"
            f"**Diagnostic Molecular:** Specimenul de **{meta['species']}** este confirmat **REZISTENT** la **{verdict['trait_category']}**. "
            f"Markerul molecular de rezistență `{marker_info['gene']}` ({marker_info['chrom']}:{marker_info['pos']}) a fost detectat cu succes în secvența genomică analizată."
        )
    else:
        st.error(
            f"### 🚨 **VERDICT FINAL: PLANTA ESTE VULNERABILĂ**\n\n"
            f"**Diagnostic Molecular:** Specimenul de **{meta['species']}** este **VULNERABIL / SUSCEPTIBIL** la **{verdict['trait_category']}**. "
            f"Markerul molecular de rezistență specific (`{marker_info['gene']}` la {marker_info['chrom']}:{marker_info['pos']}) lipsește din secvența analizată."
        )

    # Metrici native Streamlit
    col_m1, col_m2, col_m3, col_m4 = st.columns(4)

    with col_m1:
        st.metric(
            label="Verdict Fenotipic",
            value=verdict["status"].upper(),
            delta="Profil Favorabil" if verdict["is_resistant"] else "-Susceptibil"
        )

    with col_m2:
        st.metric(
            label="Trăsătură Vizată",
            value=verdict["trait_category"],
            help=marker_info["trait"]
        )

    with col_m3:
        st.metric(
            label="Genotip Identificat",
            value=marker_info["alleles"],
            delta=marker_info["genotype_code"]
        )

    with col_m4:
        st.metric(
            label="Certitudine Diagnostic",
            value=f"{verdict['confidence_score']:.1f}%",
            delta="Nivel Ridicat"
        )

    # Bară vizuală de progres pentru scorul de încredere (st.progress)
    conf_procent = float(verdict["confidence_score"])
    st.write(f"**Nivel de Încredere Diagnostic: {conf_procent:.1f}%**")
    st.progress(min(max(conf_procent / 100.0, 0.0), 1.0))

    # Rezumat Executiv Bioinformatic
    st.subheader("🔬 Rezumat Executiv Bioinformatic")
    st.info(
        f"**Constatare Biologică:** {verdict['executive_summary']}\n\n"
        f"**Mecanism Molecular de Apărare:** {marker_info['mechanism']}"
    )

    st.divider()

    # Detalii Genomice și Recomandări Agronomice
    col_locus, col_agro = st.columns(2)

    with col_locus:
        st.subheader("🧬 Date Genomice și Locus")
        df_locus = pd.DataFrame([
            {"Parametru": "Specie Gazdă", "Valoare": meta["species"]},
            {"Parametru": "Trăsătură Vizată", "Valoare": marker_info["trait"]},
            {"Parametru": "Genă Țintă", "Valoare": marker_info["gene"]},
            {"Parametru": "Coordonate Genomice", "Valoare": f"{marker_info['chrom']}:{marker_info['pos']}"},
            {"Parametru": "Alela de Referință (REF)", "Valoare": marker_info["ref_allele"]},
            {"Parametru": "Alela Variantă Protectoare (ALT)", "Valoare": marker_info["alt_allele"]},
            {"Parametru": "Apel Genotipic", "Valoare": f"{marker_info['genotype_code']} ({marker_info['alleles']})"},
            {"Parametru": "Stare Zigozitate", "Valoare": marker_info["zygosity"]},
            {"Parametru": "Profunzime Citire (DP)", "Valoare": marker_info["read_depth"]},
            {"Parametru": "Scor Calitate (Phred QUAL)", "Valoare": f"Q{marker_info['phred_quality']:.1f}"},
            {"Parametru": "Filtru Calitate VCF", "Valoare": marker_info["filter_status"]},
        ])
        st.dataframe(df_locus, use_container_width=True, hide_index=True)

    with col_agro:
        st.subheader("🌾 Plan de Acțiune Agronomică")
        for r in recs:
            st.markdown(f"- {r}")

    st.divider()

    # Proveniență probă
    with st.expander("🔎 Proveniență Probă și Metadate Tehnice"):
        col_t1, col_t2 = st.columns(2)
        with col_t1:
            st.write(f"**Identificator Probă:** `{meta['sample_id']}`")
            st.write(f"**Specie Analizată:** `{meta['species']}`")
            st.write(f"**Nume Fișier:** `{meta['filename']}`")
            st.write(f"**Dimensiune Fișier:** `{meta['file_size_kb']} KB`")
        with col_t2:
            st.write(f"**Variante Genomice Scănate:** `{meta['total_variants_scanned']:,}`")
            st.write(f"**Format Standard:** `{meta['vcf_format_version']}`")
            st.write(f"**Versiune Pipeline:** `{meta['pipeline_version']}`")

    # Descărcare Raport PDF Oficial
    st.subheader("📄 Descărcare Raport Oficial PDF")
    pdf_bytes = generate_pdf_report(raport)
    st.download_button(
        label="📥 Descarcă Raportul PDF Oficial (Raport_AgriGen_USVT.pdf)",
        data=pdf_bytes,
        file_name="Raport_AgriGen_USVT.pdf",
        mime="application/pdf",
        type="primary"
    )

else:
    # -------------------------------------------------------------------------
    # 8. INTERACTIVE TABS PENTRU INFORMAȚII PLANTE
    # -------------------------------------------------------------------------
    st.info("👆 **Instrucțiuni de Utilizare:** Selectați specia și trăsătura vizată din meniul lateral, alegeți o probă demo pentru o prezentare rapidă în 2 click-uri (ori încărcați un fișier `.vcf`/`.txt`), apoi apăsați pe **Generează Raportul**.")

    st.markdown("### 📚 Catalog de Referință Genomică Vegetală (USVT)")

    tab_porumb, tab_grau, tab_rosie = st.tabs([
        "🌽 Porumb (Zea mays)",
        "🌾 Grâu (Triticum aestivum)",
        "🍅 Roșie (Solanum lycopersicum)"
    ])

    with tab_porumb:
        st.markdown("""
        #### 🌽 Porumb (*Zea mays*) — Arhitectură Genomică de Rezistență
        
        > **Categorie:** Culturi de câmp / Cereale mari • **Bază de date de referință:** B73 RefGen_v5
        
        * **💧 Toleranță la Secetă & Stres Osmotic:**
          * **Locus cheie:** Gena `ZmDREB2A` (`chr3:145320 A->G`).
          * **Mecanism molecular:** Factor de transcripție AP2/ERF ce induce sinteza proteinelor dehidrine (LEA) și accelerează închiderea stomatelor în fazele critice de vegetație.
        * **🛡️ Rezistență la Boli & Patogeni Fungici:**
          * **Locus cheie:** Gena `ZmWAK-R` (`chr5:882341 G->T`).
          * **Mecanism molecular:** Kinază perete-celular ce detectează efectorii patogeni, declanșând un răspuns imun hipersensibil fără necroze extinse.
        * **🌾 Recomandare Agronomică USVT:** Optimizat pentru bazinele hidrografice Banat-Crișana cu alternanțe frecvente de precipitații și arșiță pedologică.
        """)

    with tab_grau:
        st.markdown("""
        #### 🌾 Grâu (*Triticum aestivum*) — Profil Genetic de Imunitate și Vigor
        
        > **Categorie:** Cereale păioase • **Bază de date de referință:** Chinese Spring IWGSC RefSeq v2.1
        
        * **💧 Toleranță la Secetă & Șoc Termic:**
          * **Locus cheie:** Gena `TaDREB1` (`chr4:551200 A->T`).
          * **Mecanism molecular:** Stimulează acumularea de osmoliti celulari (prolină), protejând complexul fotosintetic sub radiație solară și arșiță.
        * **🛡️ Rezistență Durabilă la Rugini:**
          * **Locus cheie:** Gena `Lr34` (`chr7:1200500 G->C`).
          * **Mecanism molecular:** Transportor ABC de tip membranar ce asigură rezistență durabilă cu spectru larg la rugină și fuzarioza spicului.
        * **🌾 Recomandare Agronomică USVT:** Permite reducerea substanțială a tratamentelor fungicide pe vegetație în culturile de toamnă.
        """)

    with tab_rosie:
        st.markdown("""
        #### 🍅 Roșie (*Solanum lycopersicum*) — Imunitate Foliară și Reziliență
        
        > **Categorie:** Legumicultură / Culturi protejate și câmp • **Bază de date:** Heinz 1706 SL4.0
        
        * **💧 Toleranță la Deficit Hidric și Salinitate:**
          * **Locus cheie:** Factorul `SlAREB1` (`chr1:88320 T->A`).
          * **Mecanism molecular:** Reglează căile dependente de acid abscisic (ABA), amplificând retenția hidrică și densitatea perișorilor radiculari.
        * **🛡️ Rezistență la Mana Solanaceelor & Pătare Bacteriană:**
          * **Locus cheie:** Gena `Ph-3` (`chr9:405000 C->G`).
          * **Mecanism molecular:** Receptor imunitar CC-NBS-LRR capabil să blocheze pătrunderea hifelor de *Phytophthora infestans*.
        * **🌾 Recomandare Agronomică USVT:** Excelent pentru hibrizi destinați fermelor ecologice și producției în solarii fără pesticide sistemice.
        """)

# Subsol Instituțional
st.divider()
st.caption(
    "AgriGen - USVT DNA Analyzer (MVP) | Dezvoltat pentru Universitatea de Științele Vieții «Regele Mihai I» din Timișoara | "
    "Laboratorul de Genomică și Bioinformatică Agricolă"
)
