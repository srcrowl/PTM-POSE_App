import pandas as pd
from ptm_pose import pose_config, project, flanking_sequences

import streamlit as st

#st.session_state['ptm_coordinates'] = pose_config.download_ptm_coordinates()

@st.cache_data
def convert_df(df):
    return df.to_csv().encode('utf-8')

@st.cache_data(show_spinner = False)
def project_cache(splicing_data, chromosome_col, strand_col, region_start_col, region_end_col, event_id_col, dPSI_col, sig_col, gene_col, coordinate_type):
    """
    Projection function with caching
    """
    #download ptm coordinates if not already present
    if 'ptm_coordinates' not in st.session_state:
        st.session_state['ptm_coordinates'] = pose_config.download_ptm_coordinates()

    #run projection
    splicing_data, spliced_ptms = project.project_ptms_onto_splice_events(splicing_data, ptm_coordinates = st.session_state['ptm_coordinates'], chromosome_col = chromosome_col, strand_col = strand_col, region_start_col = region_start_col, region_end_col = region_end_col, event_id_col = event_id_col, dPSI_col=dPSI_col, sig_col = sig_col, gene_col = gene_col, coordinate_type=coordinate_type, taskbar_label = 'Projecting PTMs', PROCESSES = 1)
    return splicing_data, spliced_ptms

@st.cache_data(show_spinner = False)
def get_flanking_changes(splicing_data, chromosome_col, strand_col, spliced_region_start_col, spliced_region_end_col, first_flank_start_col, first_flank_end_col, second_flank_start_col, second_flank_end_col, dPSI_col, sig_col, gene_col,  event_id_col, coordinate_type):
    """
    Flanking sequence function with caching
    """
    altered_flanks = flanking_sequences.get_flanking_changes_from_splice_data(splicing_data, ptm_coordinates = st.session_state['ptm_coordinates'], chromosome_col = chromosome_col, strand_col = strand_col, spliced_region_start_col = spliced_region_start_col, spliced_region_end_col = spliced_region_end_col, first_flank_start_col = first_flank_start_col, first_flank_end_col = first_flank_end_col, second_flank_start_col = second_flank_start_col, second_flank_end_col = second_flank_end_col, dPSI_col=dPSI_col, sig_col = sig_col, gene_col = gene_col,  event_id_col = event_id_col, coordinate_type=coordinate_type)
    return altered_flanks
    
st.title('PTM-POSE: Projecting PTMs onto splicing events')

st.markdown("Welcome to the PTM-POSE web application! This tool offers a user-friendly alternative to running PTM-POSE that requires no coding experience and will allow users to quickly annotate their data with PTMs. Simply upload your splicing dataset, indicate the columns that contain relevant information, and let PTM-POSE do the rest. This web app is best suited to smaller datasets, as larger datasets may take a long time to process and is more subject to web app failures. For larger datasets or to annotate the data with functional information, we recommend running PTM-POSE locally. For more information on PTM-POSE and how it works, see the below links.")

col1, col2 = st.columns(2)
col1.link_button('PTM-POSE Documentation and Source Code', 'https://github.com/NaegleLab/PTM-POSE')
col2.link_button('PTMs and Splicing Manuscript', 'https://www.biorxiv.org/content/10.1101/2024.01.10.575062v2')

st.header('Upload data')
splicing_file =  st.file_uploader('Upload splicing dataset (as .csv) which includes genomic coordinates for each splicing event:', type = ['csv'])

if splicing_file is not None:
    splicing_data = pd.read_csv(splicing_file)

    st.write('Splicing data uploaded successfully!')
    show_data = st.checkbox('Show dataset', key = 'show_data')
    if show_data:
        st.write(splicing_data.head())

    st.header('Indicate columns to use in dataset')
    st.write('Please indicate the columns that contain the following information (required):')
    cols = st.columns(4)
    chromosome_col = cols[0].selectbox('Chromosome column:', list(splicing_data.columns), key = 'chromosome_col', index = None, placeholder = 'Choose a column')
    strand_col = cols[1].selectbox('DNA strand:', list(splicing_data.columns), key = 'strand_col', index = False, placeholder = 'Choose a column')
    region_start_col = cols[2].selectbox('Spliced region start:', list(splicing_data.columns), key = 'region_start_col',  placeholder = 'Choose a column', index = None)
    region_end_col = cols[3].selectbox('Spliced region end:', list(splicing_data.columns), key = 'region_end_col', index = None,  placeholder = 'Choose a column')
    st.write('\n\n\n\n\n')

    if None not in [chromosome_col, strand_col, region_start_col, region_end_col]:
        coordinate_type = st.radio('Coordinate system to use:', ['hg38', 'hg19', 'hg18'], horizontal = True)

        st.write('Additional layers of information for each splicing event (optional):')
        cols = st.columns(2)
        gene_name_col = cols[0].selectbox('Gene name column:', [None]  + list(splicing_data.columns), key = 'gene_name_col', index = 0)
        event_id_col = cols[1].selectbox('Event/Region ID column:', [None] + list(splicing_data.columns), key = 'event_id_col')
        cols = st.columns(2)
        dPSI_col = cols[0].selectbox('Change in splicing (deltaPSI) column:', [None] + list(splicing_data.columns), key = 'dPSI_col')
        sig_col = cols[1].selectbox('Significance column (e.g. p-value or FDR):', [None] + list(splicing_data.columns), key = 'sig_col')
        #additional columns
        extra_cols = st.multiselect('Additional columns to include in output (can select multiple):', [col for col in splicing_data.columns if col not in [chromosome_col, strand_col, region_start_col, region_end_col, gene_name_col, event_id_col, dPSI_col, sig_col]])
    

        st.header('Run PTM-POSE')
        find_flanks = st.checkbox('I would also like to identify flanking sequences that may be altered due to splicing events')
        st.subheader('PTMs with differential inclusion')
        project_ptms = st.button('Project PTMs onto splicing events!')
        if project_ptms:
            #project ptms
            with st.spinner('Finding PTMs...'):
                st.session_state['splicing_data'], st.session_state['spliced_ptms'] = project_cache(splicing_data, chromosome_col = chromosome_col, strand_col = strand_col, region_start_col = region_start_col, region_end_col = region_end_col, event_id_col = event_id_col, dPSI_col=dPSI_col, sig_col = sig_col, gene_col = gene_name_col, coordinate_type=coordinate_type)

        if 'spliced_ptms' in st.session_state:
            #output spliced ptms
            st.success(f"Done! {st.session_state['spliced_ptms'].shape[0]} found PTMs.")
            show_spliced_ptms = st.checkbox('Show spliced PTM information', key = 'show_spliced_ptms')
            if show_spliced_ptms:
                st.write(st.session_state['spliced_ptms'])



        if find_flanks:
            st.subheader('PTMs with altered flanking sequences')
            st.write('Please indicate the columns that contain the following information (required):')
            cols = st.columns(4)
            first_flank_start_col = cols[0].selectbox('Upstream flanking region start:', list(splicing_data.columns), key = 'first_flank_start_col', index = None, placeholder = 'Choose a column')
            first_flank_end_col = cols[1].selectbox('Upstream flanking region end:', list(splicing_data.columns), key = 'first_flank_end_col', index = False, placeholder = 'Choose a column')
            second_flank_start_col = cols[2].selectbox('Downstream flanking region start:', list(splicing_data.columns), key = 'second_flank_region_start_col',  placeholder = 'Choose a column', index = None)
            second_flank_end_col = cols[3].selectbox('Downstream flanking region end:', list(splicing_data.columns), key = 'second_flank_region_end_col', index = None,  placeholder = 'Choose a column')

            #find flinking sequences
            find_flanking_sequences = st.button('Find altered flanking sequences due to splice events', key = 'find_flanking_sequences')
            if find_flanking_sequences:
                with st.spinner('Finding altered flanking sequences...'):
                    st.session_state['altered_flanks'] = get_flanking_changes(splicing_data, chromosome_col = chromosome_col, strand_col = strand_col, spliced_region_start_col = region_start_col, spliced_region_end_col = region_end_col, first_flank_start_col = first_flank_start_col, first_flank_end_col = first_flank_end_col, second_flank_start_col = second_flank_start_col, second_flank_end_col = second_flank_end_col, dPSI_col=dPSI_col, sig_col = sig_col, gene_col = gene_name_col,  event_id_col = event_id_col, coordinate_type=coordinate_type)
                    st.session_state['altered_flanks'] = st.session_state['altered_flanks'][~st.session_state['altered_flanks']['Matched']]

            if 'altered_flanks' in st.session_state:
                #output spliced ptms
                st.success(f"Done! {st.session_state['altered_flanks'].shape[0]} PTMs with altered flanking sequences.")
                show_spliced_ptms = st.checkbox('Show PTMs with altered flanks information', key = 'show_altered_flanks')
                if show_spliced_ptms:
                    st.write(st.session_state['altered_flanks'])



        st.subheader('Download outputs')
        if 'spliced_ptms' in st.session_state:
            #output annotated splicing data
            col1, col2 = st.columns(2)
            col1.download_button('Download differentially included information', data = convert_df(st.session_state['spliced_ptms']), file_name = 'spliced_ptms.csv', mime='text/csv')
            col2.download_button('Download annotated splicing data', data = convert_df(st.session_state['splicing_data']), file_name = 'annotated_splicing_data.csv', mime='text/csv')

        if 'altered_flanks' in st.session_state:
            #output annotated splicing data
            st.download_button('Download altered flank information', data = convert_df(st.session_state['altered_flanks']), file_name = 'altered_flanks.csv', mime='text/csv')

    
