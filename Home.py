# Streamlit Home

# module load in
import streamlit as st
import pandas as pd

@st.cache(show_spinner=False)
def create_main():
    path = './NDD_SMR_all.csv'
    df = pd.read_csv(path, sep = ',')

    # modifications
    df = df.rename({'Gene_rename':'annotated_gene'}, axis = 1) # rename columns as desired
    df = df.drop(columns = ['gene_clean'])
    df['p_SMR'] = df['p_SMR'].map('{:.6f}'.format)
    df['b_GWAS'] = df['b_GWAS'].map('{:.6f}'.format)
    df['se_GWAS'] = df['se_GWAS'].map('{:.6f}'.format)
    df['p_GWAS'] = df['p_GWAS'].map('{:.6f}'.format)
    df['p_eQTL'] = df['p_eQTL'].map('{:.6f}'.format)
    df['p_SMR_multi'] =  df['p_SMR_multi'].map('{:.6f}'.format)
    df['p_HEIDI'] =  df['p_HEIDI'].map('{:.6f}'.format)

    return df

# Welcome message
st.title("OmicSynth Functional NDD Gene Browser")
st.markdown("""Welcome to OmicSynth's Neurodegenrative Disorders functionl analysis browser!
            This application allows you to browse SMR data from our paper and conduct customized analysis. 
            Please report any issues, provide feedback, or ask general questions to chelsea.alvarado@nih.gov""" ) 

# load in main df and keep in cache to mitigate reloading on every page
if 'main_data' not in st.session_state: # create session state
        st.session_state['main_data'] = pd.DataFrame()

with st.spinner('Loading in data ... only happens once :)'):
    main_df = create_main()
    st.session_state['main_data'] = main_df
    st.success('Done!')