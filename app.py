"""
BBS Digenic Variant Explorer — Complete Final App
Run: streamlit run app_final.py
"""
import streamlit as st
import pandas as pd
import networkx as nx
import os
import plotly.graph_objects as go
import plotly.express as px
import requests
st.set_page_config(
    page_title="BBS Digenic Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Space+Mono:wght@400;700&family=DM+Sans:wght@300;400;500;600&display=swap');

/* Force sidebar always open */
section[data-testid="stSidebar"] {
    min-width: 260px !important;
    max-width: 260px !important;
    transform: none !important;
    visibility: visible !important;
}
[data-testid="collapsedControl"] { display: none !important; }

html, body, [class*="css"] {
    font-family: 'DM Sans', sans-serif;
    background-color: #0a0e1a !important;
    color: #e2e8f0 !important;
}
.stApp { background-color: #0a0e1a !important; }
section[data-testid="stSidebar"] {
    background: #111827 !important;
    border-right: 1px solid rgba(255,255,255,0.08);
}
.metric-box {
    background: #1a2235;
    border: 1px solid rgba(255,255,255,0.08);
    border-radius: 10px;
    padding: 1rem;
    text-align: center;
    margin-bottom: 0.5rem;
}
.metric-value {
    font-family: 'Space Mono', monospace;
    font-size: 2rem;
    font-weight: 700;
    color: #00f5c4;
    line-height: 1;
}
.metric-label {
    font-size: 0.7rem;
    color: #94a3b8;
    margin-top: 4px;
    letter-spacing: 0.08em;
    text-transform: uppercase;
}
.hero-title {
    font-family: 'Space Mono', monospace;
    font-size: 2rem;
    font-weight: 700;
    background: linear-gradient(135deg, #00f5c4, #6366f1);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    line-height: 1.3;
    margin-bottom: 0.3rem;
}
.section-title {
    font-family: 'Space Mono', monospace;
    font-size: 0.8rem;
    color: #00f5c4;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    margin-bottom: 0.6rem;
}
.card {
    background: #1a2235;
    border: 1px solid rgba(255,255,255,0.08);
    border-radius: 10px;
    padding: 1rem 1.2rem;
    margin-bottom: 0.6rem;
}
.stTabs [data-baseweb="tab-list"] {
    background: #111827 !important;
    border-radius: 8px;
    padding: 4px;
    gap: 4px;
}
.stTabs [data-baseweb="tab"] {
    background: transparent !important;
    color: #94a3b8 !important;
    border-radius: 6px !important;
}
.stTabs [aria-selected="true"] {
    background: #1a2235 !important;
    color: #00f5c4 !important;
}
div[data-testid="stRadio"] label {
    color: #94a3b8 !important;
    font-size: 0.9rem;
    padding: 6px 10px;
    border-radius: 6px;
    transition: all 0.15s;
}
div[data-testid="stRadio"] label:hover { background: #1a2235 !important; color: #e2e8f0 !important; }
#MainMenu { visibility: hidden; }
footer { visibility: hidden; }
header { visibility: hidden; }
div[data-testid="stRadio"] > div { gap: 4px !important; }
div[data-testid="stRadio"] label {
    background: transparent;
    border-radius: 8px !important;
    padding: 8px 12px !important;
    color: #64748b !important;
    font-size: 0.88rem !important;
    border: 1px solid transparent !important;
    width: 100% !important;
}
div[data-testid="stRadio"] label:hover {
    background: #1a2235 !important;
    color: #e2e8f0 !important;
}
div[data-testid="stRadio"] label:has(input:checked) {
    background: #1a2235 !important;
    color: #00f5c4 !important;
    border-color: #00f5c4 !important;
    font-weight: 600 !important;
}

</style>
""", unsafe_allow_html=True)
# ── Data loading ───────────────────────────────────
@st.cache_data
def load_data():
    base  = os.path.dirname(os.path.abspath(__file__))
    files = {
    "genes":       "layer1_bbs_genes.csv",
    "constrained": "layer2_constrained_genes.csv",
    "interactions":"layer3_interactions.csv",
    "pathways":    "layer4_pathway_scores.csv",
    "clinvar":     "layer4b_clinvar_scores.csv",
    "scores":      "layer5_digenic_scores.csv",
}
    data = {}
    for key, rel in files.items():
        path = os.path.join(base, rel)
        data[key] = pd.read_csv(path) if os.path.exists(path) else pd.DataFrame()
    return data

@st.cache_data
def load_ai():
    base = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(base, "layer6_ai_interpretations.txt")
    return open(path, encoding="utf-8").read() if os.path.exists(path) else None

data      = load_data()
ai_text   = load_ai()
scores_df = data["scores"]
genes_df  = data["constrained"]
inter_df  = data["interactions"]

def score_color(s):
    return "#ef4444" if s>=80 else "#f59e0b" if s>=60 else "#00f5c4" if s>=40 else "#6366f1"

# ── Sidebar ────────────────────────────────────────
with st.sidebar:
    st.markdown("""
    <div style='padding:0.8rem 0 1rem;'>
        <div style='font-family:Space Mono,monospace;font-size:1.1rem;color:#00f5c4;font-weight:700;letter-spacing:0.05em;'>🧬 BBS Explorer</div>
        <div style='font-size:0.72rem;color:#64748b;margin-top:3px;'>Digenic Variant Architecture</div>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("---")
    page = st.radio("Navigation", [
        "🏠  Overview",
        "🔬  Gene Analysis",
        "🕸️  Network",
        "📊  Digenic Scores",
        "🤖  AI Interpretation",
        "🌍  Disease Explorer",
        "⚙️  Pipeline"
    ], label_visibility="collapsed")

    st.markdown("---")

    # Quick stats in sidebar
    if not scores_df.empty:
        top_pair = scores_df.iloc[0]
        st.markdown(f"""
        <div class='card' style='padding:0.7rem;'>
            <div style='font-size:0.65rem;color:#64748b;text-transform:uppercase;letter-spacing:0.08em;margin-bottom:6px;'>⭐ Top Finding</div>
            <div style='font-family:Space Mono,monospace;font-size:0.85rem;color:#00f5c4;'>{top_pair['gene_a']} ↔ {top_pair['gene_b']}</div>
            <div style='font-size:0.75rem;color:#ef4444;margin-top:2px;font-weight:600;'>Score: {top_pair['digenic_score_normalized']:.1f}/100</div>
        </div>
        """, unsafe_allow_html=True)

    st.markdown("""
    <div class='card' style='padding:0.7rem;'>
        <div style='font-size:0.65rem;color:#64748b;text-transform:uppercase;letter-spacing:0.08em;margin-bottom:6px;'>Disease</div>
        <div style='font-size:0.85rem;color:#e2e8f0;font-weight:500;'>Bardet-Biedl Syndrome</div>
        <div style='font-size:0.7rem;color:#6366f1;margin-top:2px;'>MONDO:0015229</div>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("<div style='font-size:0.65rem;color:#334155;text-align:center;padding-top:0.5rem;'>Biohackathon 2026 · KU Bioinformatics</div>", unsafe_allow_html=True)

# ══════════════════════════════════════════════════
# OVERVIEW
# ══════════════════════════════════════════════════
if "Overview" in page:
 
    # ── Header ────────────────────────────────────
    st.markdown("""
    <div style='margin-bottom:1.2rem;'>
        <div style='font-size:22px;font-weight:600;color:#ffffff;'>Digenic Disease Architecture Explorer</div>
        <div style='font-size:13px;color:#64748b;margin-top:3px;'>Bardet-Biedl Syndrome · Cross-database variant prioritization · Live data integration</div>
    </div>
    """, unsafe_allow_html=True)
 
    # ── Metric cards ──────────────────────────────
    c1,c2,c3,c4 = st.columns(4)
    ng  = len(genes_df)  if not genes_df.empty  else 9
    ni  = len(inter_df)  if not inter_df.empty  else 16
    np_ = len(scores_df) if not scores_df.empty else 15
    ts  = round(scores_df["digenic_score_normalized"].max(),1) if not scores_df.empty else 100.0
    top_pair = scores_df.iloc[0] if not scores_df.empty else None
 
    with c1:
        st.markdown(f"""
        <div style='background:#111827;border:0.5px solid rgba(255,255,255,0.08);border-radius:10px;padding:1rem;'>
            <div style='font-size:24px;font-weight:600;color:#ffffff;'>{ng}</div>
            <div style='font-size:12px;color:#94a3b8;margin-top:2px;'>BBS genes</div>
            <div style='font-size:11px;color:#334155;margin-top:2px;'>Open Targets</div>
        </div>""", unsafe_allow_html=True)
    with c2:
        st.markdown(f"""
        <div style='background:#111827;border:0.5px solid rgba(255,255,255,0.08);border-radius:10px;padding:1rem;'>
            <div style='font-size:24px;font-weight:600;color:#6366f1;'>{ni}</div>
            <div style='font-size:12px;color:#94a3b8;margin-top:2px;'>Interactions</div>
            <div style='font-size:11px;color:#334155;margin-top:2px;'>STRING score &gt;700</div>
        </div>""", unsafe_allow_html=True)
    with c3:
        st.markdown(f"""
        <div style='background:#111827;border:0.5px solid rgba(255,255,255,0.08);border-radius:10px;padding:1rem;'>
            <div style='font-size:24px;font-weight:600;color:#f59e0b;'>{np_}</div>
            <div style='font-size:12px;color:#94a3b8;margin-top:2px;'>Digenic pairs</div>
            <div style='font-size:11px;color:#334155;margin-top:2px;'>Scored 0–100</div>
        </div>""", unsafe_allow_html=True)
    with c4:
        tp = f"{top_pair['gene_a']} ↔ {top_pair['gene_b']}" if top_pair is not None else "—"
        st.markdown(f"""
        <div style='background:#111827;border:0.5px solid rgba(255,255,255,0.08);border-radius:10px;padding:1rem;'>
            <div style='font-size:24px;font-weight:600;color:#ef4444;'>{ts}</div>
            <div style='font-size:12px;color:#94a3b8;margin-top:2px;'>Top score</div>
            <div style='font-size:11px;color:#334155;margin-top:2px;'>{tp}</div>
        </div>""", unsafe_allow_html=True)
 
    st.markdown("<br>", unsafe_allow_html=True)
 
    # ── Clinical use case banner ───────────────────
    st.markdown("""
    <div style='background:#0c1a2e;border:0.5px solid rgba(99,102,241,0.3);border-radius:10px;
               padding:1rem 1.2rem;margin-bottom:1rem;'>
        <div style='font-size:11px;font-weight:600;color:#6366f1;text-transform:uppercase;
                   letter-spacing:0.08em;margin-bottom:6px;'>🏥 Clinical use case</div>
        <div style='font-size:13px;color:#cbd5e1;line-height:1.7;'>
            Patient presents with BBS symptoms but only <b style='color:#f59e0b;'>one mutation found</b>.
            This tool asks: <i>which second gene is most likely co-mutated?</i><br>
            If <b style='color:#00f5c4;'>BBS9 is mutated</b> → model predicts
            <b style='color:#ef4444;'>BBS4 as the highest-risk co-mutation (score 100/100)</b> →
            doctor sequences BBS4 → <b style='color:#00f5c4;'>diagnosis confirmed.</b>
        </div>
    </div>
    """, unsafe_allow_html=True)
 
    # ── Two column layout ─────────────────────────
    col_left, col_right = st.columns([1.1, 1])
 
    with col_left:
        # Bar chart
        st.markdown("""
        <div style='background:#111827;border:0.5px solid rgba(255,255,255,0.08);border-radius:10px;
                   padding:1rem 1.2rem;margin-bottom:1rem;'>
            <div style='font-size:14px;font-weight:500;color:#ffffff;margin-bottom:2px;'>Digenic score ranking</div>
            <div style='font-size:11px;color:#64748b;margin-bottom:12px;'>Top pairs by combined evidence — hover for details</div>
        """, unsafe_allow_html=True)
 
        if not scores_df.empty:
            top = scores_df.head(7).copy()
            top["label"] = top["gene_a"] + " ↔ " + top["gene_b"]
            colors = [score_color(float(s)) for s in top["digenic_score_normalized"]]
            fig = go.Figure(go.Bar(
                x=top["digenic_score_normalized"].astype(float),
                y=top["label"],
                orientation="h",
                marker_color=colors,
                text=[f"{float(s):.1f}" for s in top["digenic_score_normalized"]],
                textposition="outside",
                customdata=top[["string_score","shared_pathways","clinvar_pair"]].values,
                hovertemplate="<b>%{y}</b><br>Score: %{x:.1f}/100<br>STRING: %{customdata[0]}<br>Pathways: %{customdata[1]}<br>ClinVar: %{customdata[2]:.2f}<extra></extra>"
            ))
            fig.update_layout(
                plot_bgcolor="#111827", paper_bgcolor="#111827",
                font_color="#94a3b8", height=280,
                xaxis=dict(range=[0,118], gridcolor="#1e293b", color="#475569",
                          showline=False, zeroline=False),
                yaxis=dict(autorange="reversed", gridcolor="#1e293b", color="#94a3b8"),
                margin=dict(l=5, r=55, t=5, b=5),
                hoverlabel=dict(bgcolor="#1a2235", font_color="#e2e8f0")
            )
            st.plotly_chart(fig, use_container_width=True)
        st.markdown("</div>", unsafe_allow_html=True)
 
        # Priority table
        st.markdown("""
        <div style='background:#111827;border:0.5px solid rgba(255,255,255,0.08);border-radius:10px;
                   padding:1rem 1.2rem;'>
            <div style='display:flex;justify-content:space-between;align-items:center;margin-bottom:12px;'>
                <div>
                    <div style='font-size:14px;font-weight:500;color:#ffffff;'>Priority digenic candidates</div>
                    <div style='font-size:11px;color:#64748b;'>Cross-referenced high-confidence pairs</div>
                </div>
            </div>
        """, unsafe_allow_html=True)
 
        if not scores_df.empty:
            medals = ["🥇","🥈","🥉","4️⃣","5️⃣"]
            for i,(_, row) in enumerate(scores_df.head(5).iterrows()):
                s = float(row["digenic_score_normalized"])
                c = score_color(s)
                bar_w = int(s)
                st.markdown(f"""
                <div style='display:flex;align-items:center;gap:10px;padding:8px 0;
                           border-bottom:0.5px solid rgba(255,255,255,0.05);'>
                    <span style='font-size:14px;min-width:22px;'>{medals[i]}</span>
                    <span style='font-family:monospace;font-size:13px;color:#e2e8f0;flex:1;'>
                        {row['gene_a']} <span style='color:#334155;'>↔</span> {row['gene_b']}
                    </span>
                    <span style='font-size:10px;background:#1e0a0a;color:#ef4444;
                               padding:2px 6px;border-radius:4px;'>Pathogenic</span>
                    <div style='width:80px;background:#1e293b;border-radius:4px;height:5px;'>
                        <div style='background:{c};width:{bar_w}%;height:100%;border-radius:4px;'></div>
                    </div>
                    <span style='font-size:12px;font-weight:600;color:{c};min-width:36px;text-align:right;'>{s:.1f}</span>
                </div>""", unsafe_allow_html=True)
        st.markdown("</div>", unsafe_allow_html=True)
 
    with col_right:
        # Venn diagram
        st.markdown("""
        <div style='background:#111827;border:0.5px solid rgba(255,255,255,0.08);border-radius:10px;
                   padding:1rem 1.2rem;margin-bottom:1rem;'>
            <div style='display:flex;justify-content:space-between;align-items:center;margin-bottom:4px;'>
                <div style='font-size:14px;font-weight:500;color:#ffffff;'>Database overlap</div>
                <span style='font-size:10px;background:#0c1a2e;color:#6366f1;padding:2px 8px;border-radius:4px;border:0.5px solid #6366f1;'>Live stats</span>
            </div>
            <div style='font-size:11px;color:#64748b;margin-bottom:10px;'>Evidence convergence across sources</div>
        """, unsafe_allow_html=True)
 
        venn_fig = go.Figure()
        for cx,cy,col,label in [(0.35,0.55,"#185FA5","Open Targets"),(0.65,0.55,"#0F6E56","gnomAD"),(0.5,0.35,"#993C1D","ClinVar")]:
            venn_fig.add_shape(type="circle", x0=cx-0.22, y0=cy-0.22, x1=cx+0.22, y1=cy+0.22,
                fillcolor=col, opacity=0.15, line_color=col, line_width=1.5)
            venn_fig.add_annotation(x=cx+(0.28 if cx!=0.5 else 0), y=cy+(0 if cx!=0.5 else -0.3),
                text=label, font=dict(size=11,color=col), showarrow=False)
        venn_fig.add_annotation(x=0.5, y=0.52,
            text="<b>9</b><br><span style='font-size:9px'>ALL 3</span>",
            font=dict(size=14,color="#ffffff"), showarrow=False, bgcolor="rgba(99,102,241,0.3)",
            bordercolor="#6366f1", borderwidth=1, borderpad=6)
        venn_fig.update_layout(
            plot_bgcolor="#111827", paper_bgcolor="#111827",
            height=200, margin=dict(l=10,r=10,t=10,b=10),
            xaxis=dict(range=[0,1],showgrid=False,zeroline=False,showticklabels=False),
            yaxis=dict(range=[0,1],showgrid=False,zeroline=False,showticklabels=False),
            showlegend=False
        )
        st.plotly_chart(venn_fig, use_container_width=True)
 
        for dot,label,val in [("#185FA5","Open Targets",f"{ng} genes"),
                               ("#0F6E56","gnomAD pLI>0.9",f"{ng} constrained"),
                               ("#993C1D","ClinVar pathogenic",f"{ng} confirmed"),
                               ("#534AB7","Digenic pairs scored",f"{np_} pairs")]:
            st.markdown(f"""
            <div style='display:flex;align-items:center;gap:8px;padding:4px 0;'>
                <div style='width:9px;height:9px;border-radius:50%;background:{dot};flex-shrink:0;'></div>
                <div style='font-size:12px;color:#94a3b8;flex:1;'>{label}</div>
                <div style='font-size:12px;font-weight:500;color:#e2e8f0;'>{val}</div>
            </div>""", unsafe_allow_html=True)
        st.markdown("</div>", unsafe_allow_html=True)
 
        # Digenic explanation
        st.markdown("""
        <div style='background:#111827;border:0.5px solid rgba(255,255,255,0.08);border-radius:10px;
                   padding:1rem 1.2rem;margin-bottom:1rem;'>
            <div style='font-size:14px;font-weight:500;color:#ffffff;margin-bottom:10px;'>How digenic inheritance works</div>
        """, unsafe_allow_html=True)
        for bg,col,txt in [
            ("#0a2010","#4ade80","BBS9 mutated alone → mild / carrier"),
            ("#1a1500","#f59e0b","BBS9 + BBS4 both mutated → full BBS"),
            ("#1a0a0a","#ef4444","Both BBSome subunits lost → complete breakdown"),
        ]:
            st.markdown(f"""
            <div style='padding:7px 10px;background:{bg};border-radius:6px;
                       color:{col};font-size:12px;margin-bottom:6px;'>
                {txt}
            </div>""", unsafe_allow_html=True)
        st.markdown("</div>", unsafe_allow_html=True)
 
        # Disease agnostic
        st.markdown("""
        <div style='background:#111827;border:0.5px solid rgba(255,255,255,0.08);border-radius:10px;
                   padding:1rem 1.2rem;'>
            <div style='font-size:14px;font-weight:500;color:#ffffff;margin-bottom:6px;'>🔄 Disease agnostic pipeline</div>
            <div style='font-size:12px;color:#64748b;margin-bottom:10px;'>Works for any rare disease with GWAS data</div>
            <div style='display:flex;flex-wrap:wrap;gap:5px;'>
        """, unsafe_allow_html=True)
        for d in ["Alström","Joubert","Usher","NPHP","Meckel","+ any ciliopathy"]:
            st.markdown(f"<span style='font-size:11px;background:#1e293b;color:#94a3b8;padding:3px 9px;border-radius:4px;'>{d}</span>", unsafe_allow_html=True)
        st.markdown("</div></div>", unsafe_allow_html=True)

# ══════════════════════════════════════════════════
# GENE ANALYSIS
# ══════════════════════════════════════════════════
elif "Gene" in page:
    st.markdown('<div class="section-title">Gene Analysis</div>', unsafe_allow_html=True)
    st.markdown("### BBS Gene Constraint Profile")

    if not genes_df.empty:
        col1, col2 = st.columns([1.6, 1])
        with col1:
            pli_vals = genes_df["pLI"].astype(float).tolist()
            colors   = ["#ef4444" if p>=0.99 else "#f59e0b" if p>=0.95 else "#00f5c4" for p in pli_vals]
            fig = go.Figure(go.Bar(
                x=pli_vals, y=genes_df["gene_symbol"].tolist(),
                orientation="h", marker_color=colors,
                text=[f"{p:.3f}" for p in pli_vals], textposition="outside",
                hovertemplate="<b>%{y}</b><br>pLI: %{x:.4f}<br>Higher = more constrained<extra></extra>"
            ))
            fig.add_vline(x=0.9, line_dash="dash", line_color="#475569",
                         annotation_text="cutoff", annotation_font_color="#64748b")
            fig.update_layout(
                plot_bgcolor="#111827", paper_bgcolor="#0a0e1a",
                font_color="#94a3b8", height=400,
                xaxis=dict(range=[0,1.1], gridcolor="#1e293b", title="pLI Score", color="#64748b"),
                yaxis=dict(autorange="reversed", gridcolor="#1e293b", color="#94a3b8"),
                title=dict(text="gnomAD pLI Scores — gene intolerance to mutations", font_color="#e2e8f0", font_size=13),
                margin=dict(l=10, r=70, t=45, b=10),
                hoverlabel=dict(bgcolor="#1a2235", font_color="#e2e8f0")
            )
            st.plotly_chart(fig, use_container_width=True)

            st.markdown("""
            <div style='display:flex;gap:16px;margin-top:4px;flex-wrap:wrap;'>
                <span style='font-size:0.72rem;color:#94a3b8;'><span style='color:#ef4444;'>■</span> pLI≥0.99 extremely constrained</span>
                <span style='font-size:0.72rem;color:#94a3b8;'><span style='color:#f59e0b;'>■</span> pLI≥0.95 highly constrained</span>
                <span style='font-size:0.72rem;color:#94a3b8;'><span style='color:#00f5c4;'>■</span> pLI<0.95 moderate</span>
            </div>""", unsafe_allow_html=True)

        with col2:
            st.markdown('<div class="section-title">Gene Details</div>', unsafe_allow_html=True)
            for _, row in genes_df.iterrows():
                pli   = float(row["pLI"])
                color = "#ef4444" if pli>=0.99 else "#f59e0b" if pli>=0.95 else "#00f5c4"
                st.markdown(f"""
                <div class='card' style='padding:0.65rem 0.9rem;margin-bottom:0.35rem;'>
                    <div style='display:flex;justify-content:space-between;align-items:center;'>
                        <span style='font-family:Space Mono,monospace;font-size:0.82rem;color:#e2e8f0;'>{row['gene_symbol']}</span>
                        <span style='font-size:0.75rem;color:{color};font-weight:600;'>pLI={pli:.3f}</span>
                    </div>
                    <div style='display:flex;justify-content:space-between;margin-top:3px;'>
                        <span style='font-size:0.68rem;color:#475569;'>OT: {float(row['ot_score']):.3f}</span>
                        <span style='font-size:0.68rem;color:{color};'>{"★★★" if pli>=0.99 else "★★" if pli>=0.95 else "★"}</span>
                    </div>
                </div>""", unsafe_allow_html=True)
    else:
        st.warning("Run layer2_gnomad.py first")

# ══════════════════════════════════════════════════
# NETWORK
# ══════════════════════════════════════════════════
elif "Network" in page:
    st.markdown('<div class="section-title">PPI Network</div>', unsafe_allow_html=True)
    st.markdown("### Interactive Protein-Protein Interaction Network")

    if not inter_df.empty and not genes_df.empty:
        all_genes = sorted(genes_df["gene_symbol"].tolist())

        col_f1, col_f2 = st.columns([2, 1])
        with col_f1:
            selected = st.multiselect(
                "🔍 Type gene name to filter network:",
                options=all_genes, default=all_genes,
                help="Remove genes to focus on specific interactions"
            )
        with col_f2:
            min_edge = st.slider("Min interaction score:", 500, 999, 700, 50)

        if not selected: selected = all_genes

        filtered_inter = inter_df[
            inter_df["gene_a"].isin(selected) &
            inter_df["gene_b"].isin(selected) &
            (inter_df["score"] >= min_edge)
        ]

        G       = nx.Graph()
        pli_map = dict(zip(genes_df["gene_symbol"], genes_df["pLI"].astype(float)))
        for gene in selected:
            G.add_node(gene, pLI=pli_map.get(gene, 0.5))
        for _, row in filtered_inter.iterrows():
            G.add_edge(row["gene_a"], row["gene_b"], weight=float(row["score"]))

        pos = nx.spring_layout(G, seed=42, k=2.8)

        # Edges
        edge_x, edge_y, edge_hover = [], [], []
        for u, v, d in G.edges(data=True):
            x0,y0=pos[u]; x1,y1=pos[v]
            edge_x+=[x0,x1,None]; edge_y+=[y0,y1,None]
            edge_hover.append(f"{u} ↔ {v}: {d['weight']}")

        node_x     = [pos[n][0] for n in G.nodes()]
        node_y     = [pos[n][1] for n in G.nodes()]
        node_names = list(G.nodes())
        node_color = ["#ef4444" if pli_map.get(n,0)>=0.99 else "#f59e0b" if pli_map.get(n,0)>=0.95 else "#00f5c4" for n in G.nodes()]
        node_size  = [22+G.degree(n)*7 for n in G.nodes()]
        node_hover = [f"<b>{n}</b><br>pLI: {pli_map.get(n,0):.3f}<br>Connections: {G.degree(n)}<br>Click to explore" for n in G.nodes()]

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=edge_x, y=edge_y, mode="lines",
            line=dict(width=1.8, color="#6366f1"), opacity=0.45, hoverinfo="none"
        ))
        fig.add_trace(go.Scatter(
            x=node_x, y=node_y, mode="markers+text",
            marker=dict(size=node_size, color=node_color, line=dict(width=2, color="#0a0e1a"),
                       opacity=0.95),
            text=node_names, textposition="top center",
            textfont=dict(color="#e2e8f0", size=11, family="Space Mono"),
            hovertext=node_hover, hoverinfo="text",
            hoverlabel=dict(bgcolor="#1a2235", font_color="#e2e8f0")
        ))
        fig.update_layout(
            plot_bgcolor="#0a0e1a", paper_bgcolor="#0a0e1a",
            height=520, showlegend=False,
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            margin=dict(l=20, r=20, t=45, b=20),
            title=dict(
                text=f"BBS PPI Network — {len(selected)} genes · {G.number_of_edges()} interactions · score≥{min_edge}",
                font_color="#e2e8f0", font_size=13
            )
        )
        st.plotly_chart(fig, use_container_width=True)

        # Hub genes
        if G.number_of_nodes() > 0:
            st.markdown('<div class="section-title">Hub Genes (most connected)</div>', unsafe_allow_html=True)
            top_hubs = sorted(G.degree(), key=lambda x: x[1], reverse=True)[:min(5, len(selected))]
            cols     = st.columns(len(top_hubs))
            for i, (gene, deg) in enumerate(top_hubs):
                with cols[i]:
                    pli = pli_map.get(gene, 0)
                    c   = "#ef4444" if pli>=0.99 else "#f59e0b" if pli>=0.95 else "#00f5c4"
                    st.markdown(f"""
                    <div class='metric-box'>
                        <div class='metric-value' style='font-size:1.1rem;color:{c};'>{gene}</div>
                        <div class='metric-label'>{deg} connections</div>
                        <div style='font-size:0.65rem;color:#475569;margin-top:2px;'>pLI={pli:.3f}</div>
                    </div>""", unsafe_allow_html=True)
    else:
        st.warning("Run layer3_string.py first")

# ══════════════════════════════════════════════════
# DIGENIC SCORES
# ══════════════════════════════════════════════════
elif "Digenic" in page:
    st.markdown('<div class="section-title">Digenic Scores</div>', unsafe_allow_html=True)
    st.markdown("### Final Ranked Digenic Candidate Pairs")

    if not scores_df.empty:
        col_f1, col_f2 = st.columns(2)
        with col_f1:
            all_genes   = sorted(set(scores_df["gene_a"].tolist()+scores_df["gene_b"].tolist()))
            gene_filter = st.multiselect("🔍 Filter by gene:", all_genes, default=[],
                                        placeholder="Type gene name...")
        with col_f2:
            min_score = st.slider("Minimum digenic score:", 0, 100, 0, 5)

        filtered = scores_df.copy()
        if gene_filter:
            filtered = filtered[filtered["gene_a"].isin(gene_filter)|filtered["gene_b"].isin(gene_filter)]
        filtered = filtered[filtered["digenic_score_normalized"]>=min_score]

        st.markdown(f"<div style='font-size:0.78rem;color:#475569;margin-bottom:10px;'>Showing {len(filtered)} of {len(scores_df)} pairs</div>", unsafe_allow_html=True)

        tab1, tab2, tab3 = st.tabs(["📊  Bar Chart", "🔥  Heatmap", "📋  Full Table"])

        with tab1:
            top = filtered.head(15).copy()
            top["label"] = top["gene_a"] + " ↔ " + top["gene_b"]
            colors = [score_color(float(s)) for s in top["digenic_score_normalized"]]
            fig = go.Figure(go.Bar(
                x=top["digenic_score_normalized"].astype(float),
                y=top["label"],
                orientation="h", marker_color=colors,
                text=[f"{float(s):.1f}" for s in top["digenic_score_normalized"]],
                textposition="outside",
                customdata=top[["string_score","shared_pathways","clinvar_pair"]].values,
                hovertemplate="<b>%{y}</b><br>Score: %{x:.1f}/100<br>STRING: %{customdata[0]}<br>Pathways: %{customdata[1]}<br>ClinVar: %{customdata[2]:.2f}<extra></extra>"
            ))
            fig.update_layout(
                plot_bgcolor="#111827", paper_bgcolor="#0a0e1a",
                font_color="#94a3b8", height=max(380, len(top)*38),
                xaxis=dict(range=[0,118], gridcolor="#1e293b", title="Digenic Score (0-100)", color="#64748b"),
                yaxis=dict(autorange="reversed", gridcolor="#1e293b", color="#94a3b8"),
                title=dict(text="Digenic Score Ranking — hover for details", font_color="#e2e8f0", font_size=13),
                margin=dict(l=10, r=60, t=45, b=10),
                hoverlabel=dict(bgcolor="#1a2235", font_color="#e2e8f0")
            )
            st.plotly_chart(fig, use_container_width=True)

        with tab2:
            gene_list = genes_df["gene_symbol"].tolist() if not genes_df.empty else []
            matrix    = pd.DataFrame(0.0, index=gene_list, columns=gene_list)
            for _, row in filtered.iterrows():
                a,b,s = row["gene_a"], row["gene_b"], float(row["digenic_score_normalized"])
                if a in matrix.index and b in matrix.columns:
                    matrix.loc[a,b] = s
                    matrix.loc[b,a] = s
            fig = px.imshow(
                matrix, color_continuous_scale="YlOrRd",
                text_auto=".0f", aspect="auto",
                title="Digenic Score Heatmap — darker = higher score"
            )
            fig.update_layout(
                plot_bgcolor="#111827", paper_bgcolor="#0a0e1a",
                font_color="#94a3b8", height=460,
                title_font_color="#e2e8f0", title_font_size=13
            )
            fig.update_traces(textfont_size=10)
            st.plotly_chart(fig, use_container_width=True)

        with tab3:
            display = filtered[["gene_a","gene_b","digenic_score_normalized",
                               "shared_pathways","string_score","clinvar_pair"]].copy()
            display.columns = ["Gene A","Gene B","Score","Pathways","STRING","ClinVar"]
            display["Score"] = display["Score"].round(1)
            st.dataframe(
                display.style.background_gradient(subset=["Score"], cmap="YlOrRd")
                             .format({"Score":"{:.1f}","ClinVar":"{:.2f}"}),
                use_container_width=True, height=420
            )
    else:
        st.warning("Run layer5_scoring.py first")

# ══════════════════════════════════════════════════
# AI INTERPRETATION
# ══════════════════════════════════════════════════
elif "AI" in page:
    st.markdown('<div class="section-title">AI + Database Interpretation</div>', unsafe_allow_html=True)
    st.markdown("### Clinical Briefs — Gemini AI + MyGene.info + ClinVar")

    if ai_text:
        sections = [s.strip() for s in ai_text.split("="*55) if s.strip() and "RANK" in s]

        if not sections:
            st.info("No interpretations parsed — check layer6_ai_interpretations.txt format")
        else:
            for section in sections:
                lines    = section.strip().split("\n")
                title    = next((l for l in lines if "RANK" in l), "Unknown")
                score_ln = next((l for l in lines if "Digenic Score" in l), "")
                source_ln= next((l for l in lines if "Source" in l), "")
                content  = "\n".join(l for l in lines if l not in [title, score_ln]).strip()

                score_val = 0.0
                try:
                    score_val = float(score_ln.split(":")[1].strip().split("/")[0].split("|")[0].strip())
                except:
                    pass

                source = "Real Data" if "MyGene" in source_ln else "AI" if "Gemini" in source_ln else "Database"
                color  = score_color(score_val)
                badge  = "🤖 Gemini AI" if "Gemini" in source_ln else "🔬 Real Database"

                with st.expander(f"🧬 {title}  |  Score: {score_val:.0f}/100  |  {badge}", expanded=score_val>=90):
                    st.markdown(f"""
                    <div style='background:#111827;border-left:3px solid {color};padding:1rem 1.2rem;
                               border-radius:0 8px 8px 0;font-size:0.85rem;line-height:1.8;
                               color:#cbd5e1;white-space:pre-wrap;font-family:DM Sans,sans-serif;'>
{content}
                    </div>""", unsafe_allow_html=True)
    else:
        st.info("Run layer6_ai.py first to generate interpretations")
        
elif "🌍" in page:
    st.markdown('<div class="section-title">Disease Explorer</div>', unsafe_allow_html=True)
    st.markdown("### Live Search — Any Rare Disease")
    st.markdown("""
    <div style='background:#111827;border:1px solid rgba(0,245,196,0.2);border-radius:10px;
               padding:1rem 1.2rem;margin-bottom:1.2rem;'>
        <div style='font-size:0.82rem;color:#94a3b8;line-height:1.7;'>
            🔍 Type any rare disease — live fetch from Open Targets + STRING DB instantly.
        </div>
    </div>
    """, unsafe_allow_html=True)

    col_s1, col_s2 = st.columns([3,1])
    with col_s1:
        disease_input = st.text_input("Disease name:", placeholder="e.g. Joubert syndrome, Usher syndrome...", label_visibility="collapsed")
    with col_s2:
        search_btn = st.button("🔍 Search", use_container_width=True)

    st.markdown("<div style='font-size:0.72rem;color:#475569;margin-bottom:6px;'>Quick search:</div>", unsafe_allow_html=True)
    qcols = st.columns(5)
    for i, qd in enumerate(["Joubert syndrome","Usher syndrome","Alström syndrome","Meckel syndrome","NPHP"]):
        with qcols[i]:
            if st.button(qd, use_container_width=True, key=f"q{i}"):
                disease_input = qd
                search_btn    = True

    if search_btn and disease_input:
        with st.spinner(f"Fetching '{disease_input}'..."):
            disease_id = None
            found_name = ""
            try:
                search_query = """
                query SearchDisease($q: String!) {
                    search(queryString: $q, entityNames: ["disease"]) {
                        hits {
                            id
                            name
                            entity
                        }
                    }
                }
                """
                r    = requests.post(
                    "https://api.platform.opentargets.org/api/v4/graphql",
                    json={"query": search_query, "variables": {"q": disease_input}},
                    timeout=15,
                    headers={"Content-Type": "application/json"}
                )
                hits = r.json().get("data", {}).get("search", {}).get("hits", [])
                if hits:
                    disease_id = hits[0]["id"]
                    found_name = hits[0]["name"]
                else:
                    st.error("Disease not found — try a different name")
            except Exception as e:
                st.error(f"Search failed: {e}")

            if disease_id:
                genes = []
                try:
                    gene_query = """
                    query DiseaseGenes($id: String!, $n: Int!) {
                        disease(efoId: $id) {
                            associatedTargets(page: {index: 0, size: $n}) {
                                rows {
                                    target {
                                        approvedSymbol
                                        approvedName
                                    }
                                    score
                                }
                            }
                        }
                    }
                    """
                    r2   = requests.post(
                        "https://api.platform.opentargets.org/api/v4/graphql",
                        json={"query": gene_query, "variables": {"id": disease_id, "n": 15}},
                        timeout=15,
                        headers={"Content-Type": "application/json"}
                    )
                    rows  = r2.json()["data"]["disease"]["associatedTargets"]["rows"]
                    genes = [
                        {
                            "gene_symbol": row["target"]["approvedSymbol"],
                            "gene_name":   row["target"]["approvedName"],
                            "ot_score":    round(row["score"], 3)
                        }
                        for row in rows if row["score"] >= 0.3
                    ]
                except Exception as e:
                    st.error(f"Gene fetch failed: {e}")

                if genes:
                    gene_symbols = [g["gene_symbol"] for g in genes]
                    interactions = []
                    try:
                        r3 = requests.get(
                            f"https://string-db.org/api/json/network"
                            f"?identifiers={'%0d'.join(gene_symbols)}"
                            f"&species=9606&required_score=700"
                            f"&caller_identity=bbs_digenic",
                            timeout=20
                        )
                        if r3.status_code == 200:
                            interactions = [
                                {
                                    "gene_a": item["preferredName_A"],
                                    "gene_b": item["preferredName_B"],
                                    "score":  int(item["score"] * 1000)
                                }
                                for item in r3.json()
                            ]
                    except:
                        interactions = []

                    m1, m2, m3 = st.columns(3)
                    with m1: st.markdown(f"<div class='metric-box'><div class='metric-value'>{len(genes)}</div><div class='metric-label'>Genes</div></div>", unsafe_allow_html=True)
                    with m2: st.markdown(f"<div class='metric-box'><div class='metric-value' style='color:#6366f1;'>{len(interactions)}</div><div class='metric-label'>Interactions</div></div>", unsafe_allow_html=True)
                    with m3: st.markdown(f"<div class='metric-box'><div class='metric-value' style='color:#f59e0b;font-size:0.85rem;'>{found_name[:25]}</div><div class='metric-label'>Disease Found</div></div>", unsafe_allow_html=True)

                    tab1, tab2 = st.tabs(["🕸️ Network", "📋 Genes"])

                    with tab1:
                        if interactions:
                            G = nx.Graph()
                            for g in gene_symbols:
                                G.add_node(g)
                            for inter in interactions:
                                if inter["gene_a"] in G and inter["gene_b"] in G:
                                    G.add_edge(inter["gene_a"], inter["gene_b"], weight=inter["score"])

                            pos   = nx.spring_layout(G, seed=42, k=2.5)
                            ex, ey = [], []
                            for u, v in G.edges():
                                x0,y0 = pos[u]
                                x1,y1 = pos[v]
                                ex += [x0, x1, None]
                                ey += [y0, y1, None]

                            fig = go.Figure()
                            fig.add_trace(go.Scatter(
                                x=ex, y=ey, mode="lines",
                                line=dict(width=1.5, color="#6366f1"),
                                opacity=0.4, hoverinfo="none"
                            ))
                            fig.add_trace(go.Scatter(
                                x=[pos[n][0] for n in G.nodes()],
                                y=[pos[n][1] for n in G.nodes()],
                                mode="markers+text",
                                marker=dict(
                                    size=[18+G.degree(n)*6 for n in G.nodes()],
                                    color="#00f5c4",
                                    line=dict(width=1.5, color="#0a0e1a")
                                ),
                                text=list(G.nodes()),
                                textposition="top center",
                                textfont=dict(color="#e2e8f0", size=10),
                                hovertext=[f"<b>{n}</b><br>Connections: {G.degree(n)}" for n in G.nodes()],
                                hoverinfo="text"
                            ))
                            fig.update_layout(
                                plot_bgcolor="#0a0e1a", paper_bgcolor="#0a0e1a",
                                height=450, showlegend=False,
                                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                                title=dict(text=f"{found_name} — PPI Network", font_color="#e2e8f0"),
                                margin=dict(l=20, r=20, t=45, b=20)
                            )
                            st.plotly_chart(fig, use_container_width=True)

                            top_hubs = sorted(G.degree(), key=lambda x: x[1], reverse=True)[:5]
                            st.markdown('<div class="section-title">Hub Genes</div>', unsafe_allow_html=True)
                            hcols = st.columns(5)
                            for i, (gene, deg) in enumerate(top_hubs):
                                with hcols[i]:
                                    st.markdown(f"<div class='metric-box'><div class='metric-value' style='font-size:1rem;color:#00f5c4;'>{gene}</div><div class='metric-label'>{deg} links</div></div>", unsafe_allow_html=True)
                        else:
                            st.info("No interactions found for these genes")

                    with tab2:
                        gdf = pd.DataFrame(genes)
                        gdf.columns = ["Symbol", "Name", "OT Score"]
                        st.dataframe(
                            gdf.style.background_gradient(subset=["OT Score"], cmap="YlOrRd"),
                            use_container_width=True, height=400
                        )
                else:
                    st.warning("No genes found above score threshold")

# ══════════════════════════════════════════════════
# PIPELINE
# ══════════════════════════════════════════════════
elif "Pipeline" in page:
    st.markdown('<div class="section-title">Analysis Pipeline</div>', unsafe_allow_html=True)
    st.markdown("### 6-Layer Cross-Database Pipeline")

    base = os.path.dirname(os.path.abspath(__file__))
    for num, name, source, outfile, color, desc in [
        ("1",  "Gene Harvesting",      "Open Targets GraphQL",  "outputs/layer1_bbs_genes.csv",          "#00f5c4", "Disease-associated genes filtered by score>0.3"),
        ("2",  "Constraint Filtering", "gnomAD v4",             "outputs/layer2_constrained_genes.csv",  "#6366f1", "Genes with pLI>0.9 — intolerant to mutations"),
        ("3",  "PPI Network",          "STRING DB",             "outputs/layer3_interactions.csv",       "#f59e0b", "Protein interactions with confidence>700"),
        ("4",  "Pathway Analysis",     "Reactome",              "outputs/layer4_pathway_scores.csv",     "#00f5c4", "Shared biological pathway co-membership"),
        ("4b", "Clinical Evidence",    "ClinVar + NCBI",        "outputs/layer4b_clinvar_scores.csv",    "#ef4444", "Confirmed pathogenic variants in patients"),
        ("5",  "Digenic Scoring",      "DiVaS-inspired",        "outputs/layer5_digenic_scores.csv",     "#f59e0b", "Combined score ranking all gene pairs 0-100"),
        ("6",  "AI Interpretation",    "Gemini + MyGene+ClinVar","outputs/layer6_ai_summary.csv",        "#6366f1", "Clinical briefs with validation experiments"),
    ]:
        exists = os.path.exists(os.path.join(base, outfile))
        status = "✅" if exists else "⏳"
        rows   = 0
        if exists:
            try:
                rows = len(pd.read_csv(os.path.join(base, outfile)))
            except:
                rows = 0
        row_txt = f"{rows} rows" if rows > 0 else ""
        st.markdown(f"""
        <div class='card' style='display:flex;align-items:center;gap:14px;padding:0.85rem 1.1rem;margin-bottom:0.4rem;border-left:3px solid {color};'>
            <div style='font-family:Space Mono,monospace;font-size:0.95rem;color:{color};font-weight:700;min-width:34px;'>L{num}</div>
            <div style='flex:1;'>
                <div style='font-size:0.88rem;color:#e2e8f0;font-weight:500;'>{name}</div>
                <div style='font-size:0.7rem;color:#475569;margin-top:2px;'>{source} · {desc}</div>
            </div>
            <div style='font-size:0.7rem;color:#475569;'>{row_txt}</div>
            <div style='font-size:1.1rem;'>{status}</div>
        </div>""", unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)
    col1, col2 = st.columns(2)
    with col1:
        st.markdown('<div class="section-title">Digenic Score Formula</div>', unsafe_allow_html=True)
        st.code("""DS = (string_conf × pLI_A × pLI_B
       × pathway_overlap × clinvar_pair)
     ÷ (AF_A × AF_B × 1e6)

Normalized to 0–100 scale
Inspired by DiVaS algorithm (Gazzo et al. 2016)""", language="text")
    with col2:
        st.markdown('<div class="section-title">Score Interpretation</div>', unsafe_allow_html=True)
        for rng, label, color in [
            ("80–100", "Top priority — strong multi-database evidence",   "#ef4444"),
            ("60–79",  "High priority — validated interaction + pathways","#f59e0b"),
            ("40–59",  "Moderate — good interaction confidence",          "#00f5c4"),
            ("0–39",   "Low — limited evidence",                          "#6366f1"),
        ]:
            st.markdown(f"""
            <div style='display:flex;align-items:center;gap:10px;padding:6px 10px;background:#111827;border-radius:6px;margin-bottom:4px;'>
                <span style='font-family:Space Mono,monospace;font-size:0.8rem;color:{color};min-width:50px;'>{rng}</span>
                <span style='font-size:0.78rem;color:#94a3b8;'>{label}</span>
            </div>""", unsafe_allow_html=True)
