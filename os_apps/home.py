import streamlit as st
import pandas as pd
import os_apps.main_df as main

def app():
    if 'main_data' not in st.session_state:
        st.session_state['main_data'] = pd.DataFrame()

    st.title('OmicSynth')

    st.write('OmicSynth aims to provide users with the ability to navigate Genome Wide Assocation Study (GWAS) summary statistics and SMR Analysis for multiple diseases. ')

    st.write('Please use the navigation bar to locate existing functionality')

    with st.spinner('Loading in data ... only happens once :)'):
        st.session_state['main_data'] = main.create_main()
        st.success('Done!')
   
   