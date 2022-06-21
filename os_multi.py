import streamlit as st
from multiapp import MultiApp
from os_apps import gene_search, home, data, analysis, top_genes, snp_search, drug_search, data_sources, network, gene_snp # import your app modules here


app = MultiApp()

st.markdown("""
# Welcome to OmicSynth!
This multi-page app is using the [streamlit-multiapps](https://github.com/upraneelnihar/streamlit-multiapps) framework developed by [Praneel Nihar](https://medium.com/@u.praneel.nihar). Also check out his [Medium article](https://medium.com/@u.praneel.nihar/building-multi-page-web-app-using-streamlit-7a40d55fa5b4).
""")

# Add all your application here
app.add_app("Home", home.app)
app.add_app("Data", data.app)
app.add_app("Analysis", analysis.app)
app.add_app("Top Genes & SNPs", top_genes.app)
app.add_app("Gene and SNP search", gene_snp.app)
#app.add_app("Gene Search", gene_search.app)
#app.add_app("SNP Search", snp_search.app)
app.add_app("Drug Target Search", drug_search.app)
app.add_app("Gene Network", network.app)
app.add_app("Data Sources and Information", data_sources.app)


# The main app
app.run()
