import pandas as pd
import streamlit as st

@st.cache(show_spinner=False)
def create_main():
    path = '/Users/alvaradocx/Documents/omicsynth/diseases/Neuro/NDDs_SMR.csv'
    df = pd.read_csv(path, sep = ',')

    # modifications
    df.rename({'Gene_rename':'annotated_gene', 'topSNP':'topSNP_full'}, axis = 1, inplace = True) # rename one gene

    snp_list = []
    for snp in df['topSNP_full']:
        if ':' in snp:
            snp_list.append(snp.split(':')[2])
        else:
            snp_list.append(snp)
    
    df['topSNP'] = snp_list

    return df
