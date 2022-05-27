import streamlit as st
import pandas as pd
import os_apps.main_df as main


def convert_df(df):
     # IMPORTANT: Cache the conversion to prevent computation on every rerun
     return df.to_csv(index = False)

def create_df(df, diseases, omics):
  # determine max possible options 
  dx_max = len(df['Disease'])
  o_max = len(df['Omic'].unique())

  # num of options in provided lists
  dx_len = len(diseases)
  o_len = len(omics)

  if dx_len == dx_max and o_max == o_len:
    return df

  elif dx_len == dx_max:
    new_df = df[df['Omic'].isin(omics)]
    return new_df

  elif o_max == o_len:
    new_df = df[df['Disease'].isin(diseases)]
    return new_df
  
  else:
    new_df = df.query("Disease in @diseases and Omic in @omics")
    return new_df

def holm_correction(df):
    adjusted_p = []
    # sort df smallest p_smr_multi to largest
    df = df.sort_values('p_SMR_multi')
    
    # total number of genes
    n = df.shape[0]
    for p in df['p_SMR_multi']:
        adjusted_p.append(p * n)
        n = n - 1
    df['adjusted_pval'] = adjusted_p
    
    return df 

def format_df(path, mtc = False):
    # format df
    df = pd.read_csv(path, sep = ',')
    df.rename({'Gene_rename':'annotated_gene'}, axis = 1, inplace = True)
    
    # order columns
    cols_start = ['Omic', 'Disease', 'annotated_gene','topSNP']
    df = df[[c for c in cols_start if c in df]
                      + [c for c in df if c not in cols_start]]
    
    if mtc:
        # run holm correction
        holm_df = holm_correction(df)
        return holm_df
    else:
        return df 

def zscore(df):
    df = df.copy()
    df['z_score'] = df.loc[:,'b_SMR']/df.loc[:,'se_SMR']

    return df

def app():
    # CSS to inject contained in a string
    hide_dataframe_row_index = """<style>.row_heading.level0 {display:none}.blank {display:none}</style>"""

    # Inject CSS with Markdown
    st.markdown(hide_dataframe_row_index, unsafe_allow_html=True)

    # session state variables for filtered dataframe
    if 'filterdf' not in st.session_state:
        st.session_state['filterdf'] = None
    if 'filter_submit' not in st.session_state:
        st.session_state['filter_submit'] = False
    if 'filter_name' not in st.session_state:
        st.session_state['filter_name'] = 'not_run'
    
    # session state variables for multiple test correction
    if 'mtc_status' not in st.session_state:
        st.session_state['mtc_status'] = 'not_run'
    if 'mtc_df' not in st.session_state:
        st.session_state['mtc_df'] = None
    if 'mtc_name' not in st.session_state:
        st.session_state['mtc_name'] = ''

    # session state variables for z-score copmutation
    if 'z_status' not in st.session_state:
        st.session_state['z_status'] = 'not_run'
    if 'z_df' not in st.session_state:
        st.session_state['z_df'] = None
    if 'z_name' not in st.session_state:
        st.session_state['z_name'] = ''

    main_df =st.session_state['main_data']

    st.title('Data Analysis')

    st.write("Here you can select data based on a selected disease and omic list. Currently only for NDDs.")

    
    with st.form("Filter_Results"):
        # diseases
        unique_dx = list(main_df['Disease'].unique())
        unique_dx.append('All')

        diseases = st.multiselect('Please select disease(s)', unique_dx)

        if "All" in diseases:
            diseases = list(main_df['Disease'].unique())

        st.write(f'**Selected Diseases**: {", ".join(diseases)}')

        # omics
        unique_omic = list(main_df['Omic'].unique())
        unique_omic.append('All')

        omics = st.multiselect('Please select omic(s)', unique_omic)

        if "All" in omics:
            omics = list(main_df['Omic'].unique())


        submitted = st.form_submit_button("Filter Results!")
        
        if submitted:
            
            result_filter_df = create_df(main_df, diseases, omics)
            st.session_state['filterdf'] = result_filter_df       
            st.session_state['filter_submit'] = 'run'     
    # You can call any Streamlit command, including custom components:
    if st.session_state['filter_submit'] == 'run':
        st.dataframe(st.session_state['filterdf'])
    
        output_name = st.text_input('Please provide an output file name if you would like to download the results', placeholder = 'example.csv')
        st.session_state['filter_name'] = output_name
        if st.session_state['filter_name']:
          st.download_button(label="Download data as CSV", data=convert_df(st.session_state['filterdf']),file_name=output_name, mime='text/csv')
    
    st.title('Analysis')
    with st.container():
      col1, col2= st.columns(2)

      with col1:
        with st.form('MTC'): # form to run holm adjustment
            st.subheader("Multiple Test Correction")
            st.write('Perform Holm (Bonferroni step-down) Correction on your filtered dataframe and add an "adjusted_pval" column')
          
            mtc_run = st.form_submit_button("Run correction!")
          
            if mtc_run: # if run selection is pressed, run adjustment on filtered df
                
                try:
                    mtc_result = holm_correction(st.session_state['filterdf'])
                    st.session_state['mtc_df'] = mtc_result
                    st.session_state['mtc_status'] = 'run'
                except AttributeError:
                    st.error('Make sure to use the options above to filter your results and then run the correction')
                    st.session_state['mtc_status'] == 'not_run'  

        if st.session_state['mtc_status'] == 'run':
            st.dataframe(st.session_state['mtc_df'])
            mtc_filename = st.text_input('Please provide an output file name for the results of the correction', placeholder = 'example.csv')
            st.session_state['mtc_name'] = mtc_filename
            
            if st.session_state['mtc_name']:
                st.download_button(label="Download data as CSV", data=convert_df(st.session_state['mtc_df']),file_name=st.session_state['mtc_name'], mime='text/csv')

      with col2:
            with st.form('zscore'):
                st.subheader("Z-score")
                st.write('Compute z-score on your filtered dataframe and add an "z-score" column')

                z_run = st.form_submit_button("Compute z-scores!")

                if z_run:
                    try:
                        z_df = zscore(st.session_state['filterdf'])
                        st.session_state['z_df'] = z_df
                        st.session_state['z_status'] = 'run'
                    except TypeError:
                        st.error('Make sure to use the options above to filter your results and then run the computation')
                        st.session_state['z_status'] == 'not_run' 

            if st.session_state['z_status'] == 'run':
                st.dataframe(st.session_state['z_df'])
                z_filename = st.text_input('Please provide an output file name for the results of the z-score computation', placeholder = 'example.csv')
                st.session_state['z_name'] = z_filename
                
                if st.session_state['z_name']:
                    st.download_button(label="Download data as CSV", data=convert_df(st.session_state['z_df']),file_name=st.session_state['z_name'], mime='text/csv')
            