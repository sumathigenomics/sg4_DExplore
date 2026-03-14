import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import mygene
import gseapy as gp

# --- BRANDING & CONFIG ---
st.set_page_config(page_title="sg4_DExplore", layout="wide", page_icon="🧬")

def local_css():
    st.markdown("""
        <style>
        .main { font-family: 'Helvetica', sans-serif; }
        .stButton>button { width: 100%; border-radius: 5px; height: 3em; background-color: #007BFF; color: white; }
        </style>
        """, unsafe_allow_html=True)

local_css()

# --- HEADER ---
st.title("sg4_DExplore")
st.subheader("Interactive RNA-Seq Discovery Portal")
st.caption("Developed by Dr. Saminathan Sivaprakasham Murugesan | lab@sumathigenomics.com")

# --- DATA LOADING ---
with st.sidebar:
    st.header("📁 Data Upload")
    res_file = st.file_uploader("Upload DESeq2 CSV", type="csv")
    counts_file = st.file_uploader("Upload Normalized Counts", type=["csv", "tsv"])
    
    st.divider()
    st.header("⚙️ Thresholds")
    lfc_thresh = st.slider("Log2 Fold Change", 0.0, 5.0, 1.0, step=0.1)
    pval_thresh = st.number_input("Adj. P-value (padj)", 0.0, 1.0, 0.05, step=0.01)

if res_file:
    df = pd.read_csv(res_file).rename(columns={'Unnamed: 0': 'Ensembl_ID'})
    df['padj'] = df['padj'].fillna(1.0)
    
    # Define Significance
    df['Status'] = 'Not Sig'
    df.loc[(df['padj'] < pval_thresh) & (df['log2FoldChange'] > lfc_thresh), 'Status'] = 'UP'
    df.loc[(df['padj'] < pval_thresh) & (df['log2FoldChange'] < -lfc_thresh), 'Status'] = 'DOWN'

    # Mapping Logic
    if st.button("🧬 1. Map Gene Symbols"):
        with st.spinner("Mapping Ensembl IDs..."):
            mg = mygene.MyGeneInfo()
            res = mg.querymany(df['Ensembl_ID'].tolist(), scopes='ensembl.gene', fields='symbol', species='human', as_dataframe=True)
            df['Symbol'] = df['Ensembl_ID'].map(res['symbol'].to_dict()) if 'symbol' in res.columns else df['Ensembl_ID']
            st.session_state['data'] = df
            st.success("Mapping Complete!")

    if 'data' in st.session_state:
        df = st.session_state['data']
        
        tab1, tab2, tab3 = st.tabs(["📊 Volcano Plot", "🔥 Heatmap", "🧬 Pathway Enrichment"])

        with tab1:
            df['-log10_padj'] = -np.log10(df['padj'].replace(0, 1e-300))
            fig = px.scatter(df, x='log2FoldChange', y='-log10_padj', color='Status', 
                             hover_name='Symbol', color_discrete_map={'UP':'red', 'DOWN':'blue', 'Not Sig':'gray'})
            st.plotly_chart(fig, use_container_width=True)

        with tab2:
            if counts_file:
                counts_df = pd.read_csv(counts_file, index_col=0)
                top_genes = df.nsmallest(50, 'padj')['Ensembl_ID']
                subset = counts_df.loc[counts_df.index.isin(top_genes)]
                z_score = subset.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
                st.plotly_chart(px.imshow(z_score, color_continuous_scale='RdBu_r', title="Top 50 Differentially Expressed Genes"))
            else:
                st.warning("Please upload normalized counts in the sidebar for heatmaps.")

        with tab3:
            if st.button("🔍 Run KEGG Analysis"):
                sig_genes = df[df['Status'] != 'Not Sig']['Symbol'].dropna().tolist()
                enr = gp.enrichr(gene_list=sig_genes, gene_sets=['KEGG_2021_Human'], organism='human', outdir=None).results
                st.dataframe(enr.head(10))
                st.plotly_chart(px.bar(enr.head(10), x='Combined Score', y='Term', orientation='h', color='Adjusted P-value'))

        st.download_button("📥 Download Annotated Results", df.to_csv(index=False), "sg4_DExplore_results.csv", "text/csv")
