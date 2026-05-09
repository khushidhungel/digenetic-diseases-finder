# digenetic-diseases-finder
"BBS Digenic Explorer is a cross-database variant prioritization pipeline that identifies the most likely digenic gene pairs in Bardet-Biedl Syndrome. It integrates Open Targets, gnomAD, STRING, Reactome, and ClinVar through a DiVaS-inspired scoring formula to rank candidate pairs by biological and clinical evidence."
# 🧬 Digenic Disease Architecture Explorer

> **Biohackathon 2026 · Kathmandu University Bioinformatics**
> Built by **Khushi Dhungel** and **Prapti Poudel**

A cross-database variant prioritization pipeline that identifies the most likely **digenic gene pairs** in rare diseases — demonstrated on **Bardet-Biedl Syndrome (BBS)**.

---

## 🔬 What is Digenic Inheritance?

Most rare disease tools look for one broken gene. But some diseases — like BBS — are caused by **two genes mutating simultaneously**. This is called digenic inheritance.

> If BBS9 is mutated alone → mild symptoms or carrier
> If BBS9 **and** BBS4 are both mutated → full Bardet-Biedl Syndrome

Our tool identifies which gene pairs, when both mutated, are most likely to cause disease — and ranks them by evidence strength.

---

## 🏥 Clinical Use Case

A patient presents with BBS symptoms but only **one mutation is found**. The doctor asks:
*"Which second gene should I sequence next?"*

This tool answers that question — ranked by a multi-database digenic score.

---

## ⚙️ The 6-Layer Pipeline

```
Layer 1  →  Gene Harvesting       Open Targets GraphQL API
Layer 2  →  Constraint Filtering  gnomAD v4 (pLI > 0.9)
Layer 3  →  PPI Network           STRING DB (confidence > 700)
Layer 4  →  Pathway Analysis      Reactome co-membership
Layer 4b →  Clinical Evidence     ClinVar pathogenicity
Layer 5  →  Digenic Scoring       DiVaS-inspired formula
Layer 6  →  AI Interpretation     Gemini AI + MyGene.info + ClinVar
```

### Digenic Score Formula (DiVaS-inspired)

```
DS = (string_conf × pLI_A × pLI_B × pathway_overlap × clinvar_pair)
     ÷ (AF_A × AF_B × 1e6)

Normalized to 0–100 scale
Reference: Gazzo et al. 2016
```

---

## 🗄️ Databases Used

| Database | Purpose |
|----------|---------|
| Open Targets | Disease-associated gene harvesting |
| gnomAD v4 | Population constraint scoring (pLI) |
| STRING DB | Protein-protein interaction network |
| Reactome | Biological pathway co-membership |
| ClinVar | Clinical pathogenicity evidence |
| Gemini AI | Clinical interpretation briefs |

All APIs are **free and open** — no API keys required except for Gemini AI (Layer 6).

---

## 🏆 Top Findings (BBS)

| Rank | Gene Pair | Score | Pathways | STRING |
|------|-----------|-------|----------|--------|
| 🥇 | BBS9 ↔ BBS4 | 100/100 | 6 | 999 |
| 🥈 | BBS2 ↔ BBS4 | 97.5/100 | 5 | 999 |
| 🥉 | BBS7 ↔ BBS4 | 92.2/100 | 5 | 999 |

**BBS4** dominates the top results — a genuine biological finding confirmed by real database queries showing BBS4's central role in the BBSome complex.

---

## 🌍 Disease Agnostic

The pipeline works for **any rare disease** with GWAS data. The Disease Explorer page lets you type any disease name and instantly fetch its gene network:

- Joubert syndrome
- Usher syndrome
- Alström syndrome
- Meckel syndrome
- NPHP
- Any ciliopathy or rare Mendelian disease

---

## 🚀 Running Locally

```bash
# Clone the repo
git clone https://github.com/khushidhungel/digenetic-diseases-finder
cd digenetic-diseases-finder

# Create virtual environment
python -m venv venv
venv\Scripts\activate  # Windows
source venv/bin/activate  # Mac/Linux

# Install dependencies
pip install -r requirements.txt

# Run the app
streamlit run app.py
```

### Optional — Gemini AI (Layer 6)

Create a `.env` file:
```
GEMINI_API_KEY=your-key-here
```

Get a free key at [aistudio.google.com](https://aistudio.google.com)

---

## 📁 Project Structure

```
digenetic-diseases-finder/
├── app.py                        ← Streamlit dashboard
├── layer1_genes.py               ← Open Targets gene harvesting
├── layer2_gnomad.py              ← gnomAD constraint filtering
├── layer3_string.py              ← STRING PPI network
├── layer4_reactome.py            ← Reactome pathway analysis
├── layer4b_clinvar.py            ← ClinVar clinical evidence
├── layer5_scoring.py             ← Digenic score calculation
├── layer6_ai.py                  ← AI interpretation
├── requirements.txt              ← Python dependencies
└── outputs/                      ← Pipeline output files
    ├── layer5_digenic_scores.csv
    ├── layer6_ai_interpretations.txt
    └── ...
```

---

## 📦 Requirements

```
streamlit
pandas
numpy
networkx
plotly
requests
seaborn
matplotlib
python-dotenv
google-generativeai
```

---

## 🔬 References

- Gazzo AM et al. (2016) DiVaS: A method for digenic variant scoring. *Human Mutation*
- Chiang AP et al. (2006) Homozygosity mapping with SNP arrays identifies TRIM32, an E3 ubiquitin ligase, as a Bardet-Biedl syndrome gene (BBS11). *PNAS*
- Forsythe E, Beales PL. (2013) Bardet-Biedl syndrome. *European Journal of Human Genetics*

---

## 👩‍💻 Authors

**Khushi Dhungel** · **Prapti Poudel**
BSc Bioinformatics · Kathmandu University
Biohackathon 2026

---

*Built with Python, Streamlit, Plotly, and five open genomic databases.*
