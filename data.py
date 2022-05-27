from turtle import down
import streamlit as st
import pandas as pd
from st_aggrid import GridOptionsBuilder, AgGrid
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
   

def app():
    # CSS to inject contained in a string
   # hide_dataframe_row_index = """<style>.row_heading.level0 {display:none}.blank {display:none}</style>"""

    # Inject CSS with Markdown
    #st.markdown(hide_dataframe_row_index, unsafe_allow_html=True)


    if 'filterdf' not in st.session_state:
        st.session_state['filterdf'] = None
    if 'filter_submit' not in st.session_state:
        st.session_state['filter_submit'] = False

    main_df = st.session_state['main_data']

    st.title('Data')

    st.write("Here you can select data based on a selected disease and omic list. Currently only for NDDs.")

    with st.sidebar:
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
            st.session_state['filter_submit'] = submitted
            if submitted:
                result_filter_df = create_df(main_df, diseases, omics)
                st.session_state['filterdf'] = result_filter_df            


        output_name = st.text_input('Please provide an output file name if you would like to download the results', placeholder = 'example.csv')
        if output_name:
            st.download_button(label="Download data as CSV", data=convert_df(st.session_state['filterdf']),file_name=output_name, mime='text/csv')

    with st.container():
        # You can call any Streamlit command, including custom components:
        if st.session_state['filter_submit']:
            gb = GridOptionsBuilder.from_dataframe(st.session_state['filterdf'])
            gb.configure_pagination(paginationAutoPageSize=True) #Add pagination
            gb.configure_side_bar() #Add a sidebar
            gb.configure_selection('multiple', use_checkbox=True, groupSelectsChildren="Group checkbox select children") #Enable multi-row selection
            gridOptions = gb.build()

            grid_response = AgGrid(
                st.session_state['filterdf'],
                gridOptions=gridOptions,
                data_return_mode='AS_INPUT', 
                update_mode='MODEL_CHANGED', 
                fit_columns_on_grid_load=False,
                theme='dark', #Add theme color to the table
                enable_enterprise_modules=True,
                height=350, 
                reload_data=True
            )

            data = grid_response['data']
            selected = grid_response['selected_rows'] 
            df = pd.DataFrame(selected)