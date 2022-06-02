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

def pull_genes_count(input_df, p_smr_multi = None, p_HEIDI = None, count = False, *genes):
    # filter df for user defined p_values if provided
    if p_smr_multi != None and p_HEIDI != None:
        filter_df = input_df.query(f'p_SMR_multi < {p_smr_multi} & p_SMR_multi != -9999.000000 & p_HEIDI > {p_HEIDI} & p_HEIDI != -9999.000000').sort_values('p_HEIDI', ascending = False)
    else:
        filter_df = input_df.query('p_SMR_multi != -9999.000000 & p_HEIDI != -9999.000000').sort_values('p_HEIDI', ascending = False)
    # narrow down to main columns of interest
    filter_df = filter_df[['Omic', 'Disease', 'Gene_rename', 'topRSID', 'b_GWAS','p_GWAS', 'b_SMR', 'p_SMR_multi', 'p_HEIDI']]
    
    # filter for genes
    final_df = pd.DataFrame()
    for gene in genes:
        df = filter_df.query(f"Gene_rename == '{gene}'")
        final_df = pd.concat([final_df,df])
    
    final_df.sort_values('p_HEIDI', ascending = False, inplace = True)
    final_df_fil = final_df.set_index(['Gene_rename', 'Disease', 'Omic']).sort_index()
    
    if count == True:
        final_df.sort_values('Omic', ascending = False, inplace = True)
        final_df_count = final_df.groupby(['Gene_rename', 'Omic']).count()
        return final_df_fil, final_df_count
    else:
        return final_df_fil

def app():
    # CSS to inject contained in a string
    #hide_dataframe_row_index = """<style>.row_heading.level0 {display:none}.blank {display:none}</style>"""

    # Inject CSS with Markdown
    #st.markdown(hide_dataframe_row_index, unsafe_allow_html=True)

    # session state variables for  gene dataframe
    if 'genes_list' not in st.session_state: # genes in dataase
        st.session_state['genes_list'] = None
    if 'gene_results_df' not in st.session_state:
        st.session_state['gene_results_df'] = None

    # session state variables for parameter buttons
    
    if 'omic_count' not in st.session_state:
        st.session_state['omic_count'] = False
    if 'dx_count' not in st.session_state:
        st.session_state['dx_count'] = False
    
    # session state variables forcount dfs
    if 'omic_df' not in st.session_state:
        st.session_state['omic_df'] = False
    if 'dx_df' not in st.session_state:
        st.session_state['dx_df'] = False
    
    
    # parameter session state
    if 'pval' not in st.session_state:
        st.session_state['pval'] = 0.05
    if 'heidi' not in st.session_state:
        st.session_state['heidi'] = 0.01
    if 'peqtl' not in st.session_state:
        st.session_state['peqtl'] = 0.05
    if 'genes' not in st.session_state:
        st.session_state['genes'] = ''
    if 'z_status' not in st.session_state:
        st.session_state['z_status'] = 'not_run'
    


    # session state variables for z-score copmutation
    
    # if 'z_df' not in st.session_state:
    #     st.session_state['mtc_df'] = None
    

    main_df = st.session_state['main_data']
    main_gene_list = list(main_df['annotated_gene'].unique())
    st.session_state['genes_list']  = main_gene_list

    st.title('Gene Search')

    st.write("Use this interactive tool to search for your gene(s) of interest against our database.")
    
    with st.form("Gene_Search"):
        st.subheader('Parameter Selection')

        # user input - genes list
        genes_list = st.text_area('Input a gene or list of genes (seperated by comma):')
        # split list
        genes_list = genes_list.split(',')
        new_genes = []
        for gene in genes_list:
            new_genes.append(gene.strip())

        st.session_state['genes'] = new_genes # add to session state

        # user -  p_val and heidi
        pvalue = st.number_input('P-value threshold', value = 0.05)
        st.session_state['pval'] = pvalue

        heidival = st.number_input('HEIDI p-value threshold', value = 0.01)
        st.session_state['heidi'] = heidival

        peqtlval = st.number_input('p-eQTL threshold', value = 0.05)
        st.session_state['peqtl'] = peqtlval


        # return count df
        count_omic = st.radio('Would you like to return a dataframe that shows how many times a gene shows up in each omic?', ('Yes', 'No'))
        st.session_state['omic_count'] = count_omic
        count_dx = st.radio('Would you like to return a dataframe that shows how many times a gene shows up in each disease?', ('Yes', 'No'))
        st.session_state['dx_count'] = count_dx


        submitted = st.form_submit_button("Search!")
        
        
    if submitted:
        st.session_state['z_status'] = 'run'

    if st.session_state['z_status'] == 'run':
        # make sure gene list is correct
        not_genes = []
        for gene in st.session_state['genes']:
            if gene not in st.session_state['genes_list']:
                not_genes.append(gene)
        
        if len(not_genes) >= 2:
            st.write(f'{", ".join(not_genes)} are not recognized or are not in our dataset. Please check your list and try again.')
        elif len(not_genes) == 1:
            st.write(f'{not_genes[0]} is not recognized or is not in our dataset. Please check your gene and try again.')
        else:
            gene_results_df = main_df[main_df['annotated_gene'].isin(st.session_state['genes'])]
            p_val = st.session_state['pval']
            heidi_val = st.session_state['heidi']
            gene_results_df = gene_results_df.query(f'p_SMR_multi < {p_val} & p_HEIDI > {heidi_val} & p_eQTL > {peqtl}')
            st.session_state['gene_results_df'] = gene_results_df
            st.subheader('Results dataframe')
            st.dataframe(st.session_state['gene_results_df'])

            st.download_button(label="Download results as CSV", data=convert_df(st.session_state['gene_results_df']), mime='text/csv')

        if st.session_state['omic_count'] == 'Yes' and st.session_state['dx_count'] == 'Yes':
            # omic count
            omic_count_df = pd.crosstab(st.session_state['gene_results_df'].Omic, st.session_state['gene_results_df'].annotated_gene)
            st.session_state['omic_df'] = omic_count_df
            st.subheader('Omic count')
            st.dataframe(st.session_state['omic_df'])
            
            st.download_button(label="Download omic count results as CSV", data=convert_df(st.session_state['omic_df']), mime='text/csv')

            # disease count
            dx_count_df = pd.crosstab(st.session_state['gene_results_df'].Disease, st.session_state['gene_results_df'].annotated_gene)
            st.session_state['dx_df'] = dx_count_df
            st.subheader('Disease count')
            st.dataframe(st.session_state['dx_df'])

            st.download_button(label="Download disease count results as CSV", data=convert_df(st.session_state['dx_df']), mime='text/csv')

        elif st.session_state['omic_count'] == 'Yes' and st.session_state['dx_count'] == 'No':
            # omic count
            omic_count_df = pd.crosstab(st.session_state['gene_results_df'].Omic, st.session_state['gene_results_df'].annotated_gene)
            st.session_state['omic_df'] = omic_count_df
            st.subheader('Omic count')
            st.dataframe(st.session_state['omic_df'])

            st.download_button(label="Download omic count results as CSV", data=convert_df(st.session_state['omic_df']), mime='text/csv')

        elif st.session_state['omic_count'] == 'No' and st.session_state['dx_count'] == 'Yes':
            # disease count
            dx_count_df = pd.crosstab(st.session_state['gene_results_df'].Disease, st.session_state['gene_results_df'].annotated_gene)
            st.session_state['dx_df'] = dx_count_df
            st.subheader('Disease count')
            st.dataframe(st.session_state['dx_df'])

            st.download_button(label="Download disease count results as CSV", data=convert_df(st.session_state['dx_df']), mime='text/csv')
