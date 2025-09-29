import base64
import datetime

import streamlit as st
# from backend import ResearchBackend
import os
import pandas as pd

import streamlit as st
import py3Dmol

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis import rms
import seaborn as sns

from Bio.PDB import Superimposer
from st_aggrid import AgGrid, GridOptionsBuilder, JsCode, GridUpdateMode, ColumnsAutoSizeMode
# from st_table_select_cell import st_table_select_cell

from Bio.PDB import PDBParser

import pymol2

import plotly.graph_objects as go
from io import BytesIO

import math
import requests
import io
from difflib import SequenceMatcher


STOP_FLAG = False
FIG_WIDTH = 400
VISITOR_COUNT = 0

def _get_go_category_from_id(go_id: str) -> str:
        """Get GO category using QuickGO API."""
        if not go_id or not go_id.startswith('GO:'):
            return 'Unknown'
            
        url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{go_id}"
        
        try:
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            data = response.json()
            
            if 'results' in data and data['results']:
                aspect = data['results'][0].get('aspect', '')
                category_map = {
                    'biological_process': 'BP',
                    'molecular_function': 'MF',
                    'cellular_component': 'CC'
                }
                return category_map.get(aspect, 'Unknown')
                
        except:
            pass
        
        return 'Unknown'


def _extract_feature_evidence(feature):
    if 'evidences' in feature:
        return ', '.join([ev.get('evidenceCode', '') for ev in feature['evidences']])
    return 'Unknown'

def _extract_features(data):
    """Extract sequence features from UniProt data."""
    features = []
    
    if 'features' in data:
        for feature in data['features']:
            feature_data = {
                'type': feature.get('type', ''),
                'category': feature.get('category', ''),
                'description': feature.get('description', ''),
                'start': None,
                'end': None,
                'evidence': _extract_feature_evidence(feature)
            }
            
            # Extract position information
            if 'location' in feature:
                location = feature['location']
                if 'start' in location and 'end' in location:
                    feature_data['start'] = location['start'].get('value')
                    feature_data['end'] = location['end'].get('value')
            
            if feature_data['start'] and feature_data['end']:
                features.append(feature_data)
    
    return features

def _extract_go_terms(data):
    go_terms = []
    if 'uniProtKBCrossReferences' in data:
        for ref in data['uniProtKBCrossReferences']:
            if ref['database'] == 'GO':
                go_id = ref['id']
                properties = {prop['key']: prop['value'] for prop in ref.get('properties', [])}
                evidences = ref.get('evidences', [])[0] if len(ref.get('evidences', [])) > 0 else {}
                
                go_terms.append({
                    'go_id': go_id,
                    'term': properties.get('GoTerm', ''),
                    'evidence': properties.get('GoEvidenceType', ''),
                    'source': evidences.get('source', ''),
                    'category': _get_go_category_from_id(go_id)
                })
    
    return go_terms

def _extract_keywords(data):
        """Extract keywords from UniProt data."""
        keywords = []
        if 'keywords' in data:
            for keyword in data['keywords']:
                keywords.append(keyword.get('name', ''))
        return keywords

def _extract_comments(data):
    """Extract functional comments from UniProt data."""
    comments = []
    if 'comments' in data:
        for comment in data['comments']:
            comment_type = comment.get('commentType', '')
            if comment_type in ['FUNCTION', 'CATALYTIC ACTIVITY', 'SUBCELLULAR LOCATION', 'DOMAIN']:
                texts = comment.get('texts', [])
                for text in texts:
                    comments.append({
                        'type': comment_type,
                        'text': text.get('value', ''),
                        'evidence': text.get('evidences', [])
                    })
    return comments

def _extract_protein_name(data):
    try:
        return data['proteinDescription']['recommendedName']['fullName']['value']
    except:
        return 'Unknown'
    
def _extract_organism(data):
    try:
        return data['organism']['scientificName']
    except:
        return 'Unknown'

def get_go_terms_from_uniprot(uniprot_id: str):
        """Retrieve GO terms for a UniProt ID."""
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        data = response.json()

        protein_data = {
            'uniprot_id': uniprot_id,
            'sequence': data.get('sequence', {}).get('value', ''),
            'length': data.get('sequence', {}).get('length', 0),
            'protein_name': _extract_protein_name(data),
            'organism': _extract_organism(data),
            'features': _extract_features(data),
            'go_terms': _extract_go_terms(data),
            'keywords': _extract_keywords(data),
            'comments': _extract_comments(data)
        }
        
        return protein_data


def show_pdb_with_disorder(pdb_path, width=FIG_WIDTH):
    viewer = py3Dmol.view(width=width, height=width)

    if pdb_path.endswith(".pdb"):
        # Parse with Biopython
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", pdb_path)

        # Collect B-factor info
        disordered_residues = set()
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        if atom.bfactor > 90.0:  # Threshold, tweak as needed
                            disordered_residues.add((chain.id, residue.id[1]))

        pdb_data = load_pdb_file(pdb_path)
        viewer.addModel(pdb_data, "pdb")  # Add the PDB model
    else:
        # Load your 3D TIF file
        with open(pdb_path, 'r') as f:
            pdb_data = f.read()
        viewer.addModel(pdb_data, "cif")  # Add the PDB model
    
    viewer.setStyle({"cartoon": {"color": "spectrum", "colorValue": "bfactor"}})  # Color based on B-factor
    viewer.zoom(0.05)  # Auto-zoom to fit the structure to the viewer size
    viewer.setBackgroundColor("white")
    viewer.show()

    return viewer, pdb_data

def show_pdb_align(pdb1_path, pdb2_path, width=FIG_WIDTH):
    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(pdb1_path, "mol1")
        pymol.cmd.load(pdb2_path, "mol2")
        rmsd = pymol.cmd.align("mol1 and name CA", "mol2 and name CA")[0]
        pymol.cmd.save("tmp/mol1.pdb", "mol1")
        pymol.cmd.save("tmp/mol2.pdb", "mol2")

    viewer = py3Dmol.view(width=width, height=width)

    # Add first structure in red
    pdb1 = load_pdb_file("tmp/mol1.pdb")
    pdb2 = load_pdb_file("tmp/mol2.pdb")

    viewer.addModel(pdb1, "pdb")
    viewer.setStyle({'model': 0}, {"cartoon": {"color": "red"}})

    # Add second structure in blue
    viewer.addModel(pdb2, "pdb")
    viewer.setStyle({'model': 1}, {"cartoon": {"color": "blue"}})

    viewer.zoom(0.05)  # Auto-zoom

    os.remove("tmp/mol1.pdb")
    os.remove("tmp/mol2.pdb")

    return viewer

def get_alignment_score(pdb_data1, pdb_data2):
    # Load the two structures
    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(pdb_data1, "mol1")
        pymol.cmd.load(pdb_data2, "mol2")

        # # Align mol2 to mol1 and get RMSD
        # # align() returns a list: [RMSD, aligned atoms, target atoms, RMS deviation, rotation matrix, translation vector]
        rmsd_info = pymol.cmd.align("mol2", "mol1")

    return round(rmsd_info[0], None)

def generate_color_map():
    values = np.linspace(0, 100, 256)  # Range from 0 to 100 for B-factor
    cmap = plt.get_cmap('viridis')  # You can change this colormap to 'viridis', 'plasma', etc.
    rgba_values = cmap(values / 100)  # Normalize to [0, 1] range

    # Create a gradient bar plot
    fig, ax = plt.subplots(figsize=(80, 1))
    cbar = ax.imshow([values], cmap='viridis', aspect='auto')
    ax.set_axis_off()
    # cbar_ax = fig.colorbar(cbar, ax=ax, orientation='vertical', ticks=[], drawedges=False)
    # cbar_ax.set_ticks([])  # Remove ticks
    # cbar_ax.set_xticks([])  # Remove x-axis ticks
    # cbar_ax.set_yticks([])  # Remove y-axis ticks

    # Remove the x and y axis labels
    # cbar_ax.set_ticklabels([])  # Remove x-axis labels

    # Add text annotations
    ax.text(0.01, 1.5, "Low", fontsize=80, color="black", verticalalignment="center", horizontalalignment="center", transform=ax.transAxes)
    ax.text(0.99, 1.5, "High", fontsize=80, color="black", verticalalignment="center", horizontalalignment="center", transform=ax.transAxes)

    # Adjust the plot margins to ensure the text is visible
    # plt.subplots_adjust(left=0.1, right=0.9, top=0.8, bottom=0.2)
    st.pyplot(fig)

def load_pdb_file(filepath):
    with open(filepath, "r") as file:
        pdb_data = file.read()
    return pdb_data

def get_seq_path(target_seq, task, model_name):
    pdb_path = "/{}/{}/{}_model.pdb".format(task, model_name, target_seq.upper())  # Change this to your actual directory path
    # print(pdb_path)
    if not os.path.exists(pdb_path):
        pdb_path = "/{}/{}/{}_model.pdb".format(task, model_name, target_seq.lower())
        if not os.path.exists(pdb_path):
            pdb_path = "/{}/{}/{}_model.cif".format(task, model_name, target_seq.upper())
            if not os.path.exists(pdb_path):
                pdb_path = "/{}/{}/{}_model.cif".format(task, model_name, target_seq)
                if not os.path.exists(pdb_path):
                    if model_name == "AF2":
                        pdb_path = "/{}/{}/{}.pdb".format(task, model_name, target_seq.upper())
                        if not os.path.exists(pdb_path):
                            pdb_path = "/{}/{}/{}.pdb".format(task, model_name, target_seq.lower())
                            if not os.path.exists(pdb_path):
                                st.info("Seq {} not found in {}, please try again later!".format(target_seq.upper(), model_name))
                                return None
                    elif model_name == "AF3":
                        pdb_path = "/{}/{}/{}.pdb".format(task, model_name, target_seq.upper())
                        if not os.path.exists(pdb_path):
                            pdb_path = "/{}/{}/{}.pdb".format(task, model_name, target_seq.lower())
                            if not os.path.exists(pdb_path):
                                st.info("Seq {} not found in {}, please try again later!".format(target_seq.upper(), model_name))
                                return None
                    else:
                        return None

    return pdb_path

def get_id_path(target_id, task, model_name):
    pdb_path = "/{}/{}/{}_model.pdb".format(task, model_name, target_id.upper())  # Change this to your actual directory path

    if not os.path.exists(pdb_path):
        pdb_path = "/{}/{}/{}_model.pdb".format(task, model_name, target_id)
        if not os.path.exists(pdb_path):
            pdb_path = "/{}/{}/{}_model.cif".format(task, model_name, target_id.upper())
            if not os.path.exists(pdb_path):
                pdb_path = "/{}/{}/{}_model.cif".format(task, model_name, target_id)
                if not os.path.exists(pdb_path):
                    if model_name == "AF3":
                        pdb_path = "/{}/{}/{}.pdb".format(task, model_name, target_id.lower())
                        if not os.path.exists(pdb_path):
                            st.info("Seq {} not found in {}, please try again later!".format(target_id.upper(), model_name))
                            return None

    return pdb_path

def get_prediction_path(task, model_name):
    if task == "PPI":
        pdb_path = "/{}/{}/{}.csv".format(task, model_name, model_name)  # Change this to your actual directory path

        if not os.path.exists(pdb_path):
            st.info("Prediction results of {} task in model {} are not available now, please try again later!".format(task, model_name))
            return None
    
    elif task == "Drug":
        pdb_path = "/{}/Prediction_LISA/{}/5HT1A.csv".format(task, model_name)  # Change this to your actual directory path

        if not os.path.exists(pdb_path):
            st.info("Prediction results of {} task in model {} are not available now, please try again later!".format(task, model_name))
            return None

    return pdb_path

def get_server_path(target_seq, task):
    full_path = "/{}/Full/tmp/{}.cif".format(task, target_seq.upper().replace("-", "_"))  # Change this to your actual directory path

    if not os.path.exists(full_path):
        full_path = "/{}/Full/tmp/{}.cif".format(task, target_seq.replace("-", "_").lower())
        if not os.path.exists(full_path):
            # st.error("Seq {} not found".format(target_seq.upper()))
            full_path = None

    disorder_path = "/{}/Disorder/tmp/{}_disorder.cif".format(task, target_seq.upper().replace("-", "_"))  # Change this to your actual directory path

    if not os.path.exists(disorder_path):
        disorder_path = "/{}/Disorder/tmp/{}_disorder.cif".format(task, target_seq.replace("-", "_").lower())
        if not os.path.exists(disorder_path):
            # st.error("Seq {} disorder file not found".format(target_seq.upper()))
            disorder_path = None

    disorder_path2 = "/{}/Disorder/tmp/{}_disorder_2.cif".format(task, target_seq.upper().replace("-", "_"))  # Change this to your actual directory path

    if not os.path.exists(disorder_path2):
        disorder_path2 = "/{}/Disorder/tmp/{}_disorder_2.cif".format(task, target_seq.replace("-", "_").lower())
        if not os.path.exists(disorder_path2):
            # st.error("Seq {} not found".format(target_seq.upper()))
            disorder_path2 = None

    return full_path, disorder_path, disorder_path2

def make_aggrid(df, fields):
    df = df.copy()

    # casting to string wasn't needed for streamlit-aggrid<=0.3.4
    # df.style.highlight_max(axis=0)
    df = df.astype(str)

    gb = GridOptionsBuilder.from_dataframe(df)
    gb.configure_grid_options(domLayout='autoHeight')
    for field in fields:
        gb.configure_column(field, cellRenderer=BtnCellRenderer)
    grid_options = gb.build()

    response = AgGrid(
        df,
        gridOptions=grid_options,
        allow_unsafe_jscode=True,
        try_to_convert_back_to_original_types=False,  # otherwise we lose [clicked] strings in numerical columns
        columns_auto_size_mode=ColumnsAutoSizeMode.FIT_CONTENTS,
        height=math.ceil(40 * len(fields)),
        fit_columns_on_grid_load=True
    )

    return response

def generate_rmsd_color_map():
    red_colors = np.linspace(0, 1, 100)
    colors = [(1 - color, 128.0 / 256.0, 50.0 / 256) for color in red_colors]

    # Create the ListedColormap
    my_cmap = ListedColormap(colors)

    # Create a gradient bar plot
    fig, ax = plt.subplots(figsize=(100, 1))
    cbar = ax.imshow([np.linspace(0, 100)], cmap=my_cmap, aspect='auto')
    ax.set_axis_off()

    # Add text annotations
    ax.text(0.01, 1.5, "Low", fontsize=80, color="black", verticalalignment="center", horizontalalignment="center", transform=ax.transAxes)
    ax.text(0.99, 1.5, "High", fontsize=80, color="black", verticalalignment="center", horizontalalignment="center", transform=ax.transAxes)

    # Adjust the plot margins to ensure the text is visible
    # plt.subplots_adjust(left=0.1, right=0.9, top=0.8, bottom=0.2)
    st.pyplot(fig)

BtnCellRenderer = JsCode(
    """
    class BtnCellRenderer {
        init(params) {
            this.params = params;
            this.eGui = document.createElement('div');
            this.eGui.style.position = 'relative';  // To position elements within the cell
            const bgColor = this.getHeatmapColor(this.params.value);
            this.eGui.style.backgroundColor = bgColor;
            this.NoMoreClick = false;

            if (String(this.params.value).includes('[clicked]')) {
                this.params.value = this.params.value.replace('[clicked]','');
                this.params.originalValue = this.params.value;
                this.makeButton('📍');
            } else {
                this.params.originalValue = this.params.value;
                this.makeButton('🔍');
            }
        }

        makeButton(symbol) {
            this.destroy();
            this.eGui.innerHTML = `
                <span id='click-button' style="position: absolute; left: 0; top: 50%; transform: translateY(-50%);">
                    ${symbol}
                </span>
                <span style="margin-left: 20px;">${this.params.value}</span>`;
            this.eButton = this.eGui.querySelector('#click-button');
            this.btnClickedHandler = this.btnClickedHandler.bind(this);
            this.eButton.addEventListener('click', this.btnClickedHandler);
        }

        getGui() { return this.eGui; }

        refresh() { return true; }

        destroy() {
            if (this.eButton) { this.eGui.removeEventListener('click', this.btnClickedHandler); }
        }

        refreshTable(value) { this.params.setValue(value); }

        btnClickedHandler(event) {
            if(String(this.params.getValue()).includes('[clicked]')) {
                this.refreshTable(this.params.originalValue);
                this.makeButton('🔍');
            } else {
                this.destroy();
                this.refreshTable('[clicked]'+this.params.originalValue);
                this.makeButton('📍');
            }
        }

        getHeatmapColor(value) {
            // Normalize value from 0 (bad) to 100 (good)
            const min = 0;
            const max = 100;
            const val = Math.min(Math.max(value, min), max);  // clamp
            const percent = (val - min) / (max - min);

            const r = Math.round(256 * (1 - percent));
            const g = 128;
            const b = 50;

            return `rgb(${r},${g},${b})`;
        }
    };
    """
)

def get_prediction_result(task, prediction_path, target_seq):
    if task == "PPI":
        if not os.path.exists(prediction_path):
            st.info("Prediction results of {} task is not available now, please try again later!".format(task))
            return
        df = pd.read_csv(prediction_path, header=0)

        target_seq1 = target_seq.split("-")[0].upper()
        target_seq2 = target_seq.split("-")[1].upper()

        result = df.loc[(df['Seq_A'] == target_seq1) & (df['Seq_B'] == target_seq2), 'Label']
        if result is None or len(result) == 0:
            result = df.loc[(df['Seq_A'] == target_seq2) & (df['Seq_B'] == target_seq1), 'Label']

            if result is None or len(result) == 0:
                return None
    elif task == "Drug":
        if not os.path.exists(prediction_path):
            st.info("Prediction results of {} task is not available now, please try again later!".format(task))
            return
        df = pd.read_csv(prediction_path, header=0)
        result = df.head()
        
    return result

def show_page_1():
    st.session_state.result = []
    models_list = ["AF2", "AF3", "Boltz", "Chai", "OpenFold", "Proteinx", "UniFold"]

    global STOP_FLAG

    st.markdown("""
        <style>
        /* Outer div controlling the input box */
        div[data-baseweb="input"] {
            height: 300px; /* set height of the container */
        }

        /* The div inside that centers the input field */
        div[data-baseweb="input"] > div {
            height: 300px; /* match height */
            display: flex;
            align-items: flex-start;
            padding-top: 10px;
        }

        # /* The actual text input */
        # input[type="text"] {
        #     font-size: 24px !important;
        # }
        # </style>
    """, unsafe_allow_html=True)

    st.info("Enter protein sequence and select models to compare.", icon="📝")
    col1, col2 = st.columns(2, gap="small")

    with col1:
        st.markdown("⌨️**Enter protein sequence**", help="Please follow the protein-protein form as the input sequence. \nExample: A6NIX2-O60663")
        target_seq = st.text_input(label="Enter the protein sequence below", value="A6NIX2-O60663")

        if not "-" in target_seq:
            st.error("🚫The input sequence format contains mistakes. Please check the input sequence!")
            STOP_FLAG = True

    progress_bar = st.progress(0)
    progress_text = st.empty()

    if not STOP_FLAG:
        with col2:
            st.markdown("📋**Select models to compare**", help="Please choose at least two models!")

            options_list = []
            for index, model in enumerate(models_list):
                options_list.append(st.checkbox(model, value=True))

            select_models_index = []
            for index, option in enumerate(options_list):
                progress_bar.progress(index)
                progress_text.text("{}% - Initializing...".format(index))
                if int(option) == 1:
                    select_models_index.append(index)
            
            models_selected_list = []
            for model_index in select_models_index:
                models_selected_list.append(models_list[model_index])

            STOP_FLAG = False
    
    if not STOP_FLAG:
        df_path = pd.DataFrame(columns=["model_name", "pdb_path"])
        for index, model_name in enumerate(models_selected_list):
            progress_bar.progress(15 + index)
            progress_text.text("{}% - Loading models...".format(str(15 + index)))

            tmp_path = get_seq_path(target_seq, 'PPI', model_name)

            if tmp_path is None:
                st.error("File not found for {}, please deselect {}".format(model_name, model_name))
                STOP_FLAG = True
                st.stop()
            df_path.loc[df_path.shape[0]] = [model_name, tmp_path]

    if not STOP_FLAG:
        df_alignment_score = pd.DataFrame(columns=models_selected_list)
        
        p_i = 1
        for row, model_name_row in enumerate(models_selected_list):
            for col, model_name_col in enumerate(models_selected_list):
                progress_bar.progress(30 + int(p_i / 10))
                progress_text.text("{}% - Calculating RMSD...".format(30 + int(p_i / 10)))
                df_alignment_score.loc[row, model_name_col] = get_alignment_score(df_path.loc[row, "pdb_path"], df_path.loc[col, "pdb_path"])
                p_i += 1

    progress_bar.progress(50)
    progress_text.text("{}% - Drawing heatmap...".format(50))

    if not STOP_FLAG:
        rowId = 0
        colIndex = 0

        if len(select_models_index) > 1:
            st.markdown("**Aligment RMSD Heatmap**", help="The RMSD comes from PyMOL package.")
            st.info("The numerical value in the heatmap indicates the RMSD between the two models. Click icon 🔍 before the numerical value to indicate the model pairs to compare.")

            generate_rmsd_color_map()
            
            df_tmp = pd.DataFrame(columns=['Models'] + models_selected_list)
            p_i = 1
            for model_index, model_name in enumerate(models_selected_list):
                alignment_score_list = []
                for index, tmp_model_name in enumerate(models_selected_list):
                    progress_bar.progress(50 + p_i)
                    progress_text.text("{}% - Loading RMSD values...".format(str(50 + p_i)))
                    alignment_score_list.append(df_alignment_score.loc[model_index, tmp_model_name])
                    p_i += 1

                df_tmp.loc[df_tmp.shape[0]] = [model_name] + alignment_score_list

            response = make_aggrid(df_tmp, models_selected_list)

            progress_bar.progress(90)
            progress_text.text("{}% - Checking model pairs...".format(str(90)))

            df = response["data"]

            # st.write(st.session_state.result)
            for col in df.columns:
                for idx in df.index:
                    if isinstance(df.at[idx, col], str):
                        try:
                            rowId = models_list.index(models_selected_list[int(idx)])
                            colIndex = models_list.index(col)
                        except Exception as e:
                            rowId = 0
                            colIndex = 0
                        
                        if df.at[idx, col].startswith("[clicked]"):
                            # st.write(rowId, colIndex)
                            if (rowId, colIndex) not in st.session_state.result:
                                st.session_state.result.append((rowId, colIndex))
                        else:
                            if (rowId, colIndex) in st.session_state.result:
                                st.session_state.result.remove((rowId, colIndex))
                        # break
            progress_bar.progress(100)
            progress_text.text("{}% - Complete!".format(str(100)))

        # Visulization
        if len(st.session_state.result) > 0:
            rowId = int(st.session_state.result[-1][0])
            colIndex = int(st.session_state.result[-1][1])
        else:
            rowId = 0
            colIndex = 0

        protein_data0 = get_go_terms_from_uniprot(target_seq.split("-")[0].upper())
        protein_data1 = get_go_terms_from_uniprot(target_seq.split("-")[1].upper())

        def match_go_to_regions(regions, go_terms, threshold=0.4):
            region_to_go = []
            for region in regions:
                matched_go = []
                for go in go_terms:
                    score = SequenceMatcher(None, region['description'].lower(), go['term'].lower()).ratio()
                    if score > threshold:
                        matched_go.append(go)
                region_to_go.append({
                    "start": region["start"],
                    "end": region["end"],
                    "description": region["description"],
                    "go_terms": matched_go
                })
            return region_to_go

        go_position_mapping0 = match_go_to_regions(protein_data0['features'], protein_data0['go_terms'], threshold=0.4)
        go_position_mapping1 = match_go_to_regions(protein_data1['features'], protein_data1['go_terms'], threshold=0.4)

        if (rowId > -1) & (colIndex > -1) & (rowId < len(models_selected_list)) & (colIndex < len(models_selected_list)):
            st.info("B-factor Color Mapping", icon="🌈")
            generate_color_map()  # Generate and display the color map

            if colIndex == rowId:
                pdb_path_col1 = df_path.loc[rowId, 'pdb_path']
                model_name_col1 = df_path.loc[rowId, 'model_name']
                st.info("{}".format(model_name_col1), icon="🧬")
                viewer, pdb_content_col1 = show_pdb_with_disorder(pdb_path_col1)

                buffer = BytesIO()
                buffer.write(pdb_content_col1.encode('utf-8'))
                buffer.seek(0)

                if viewer is not None:
                    st.components.v1.html(viewer._make_html(), height=FIG_WIDTH, width=FIG_WIDTH)

                    if pdb_path_col1.endswith(".pdb"):
                        st.download_button(
                            label="📥 Download PDB file",
                            data=buffer,
                            file_name="{}_{}_structure.pdb".format(target_seq, model_name_col1),
                            mime='Biomedical/x-pdb'
                        )
                    else:
                        st.download_button(
                            label="📥 Download CIF file",
                            data=buffer,
                            file_name="{}_{}_structure.cif".format(target_seq, model_name_col1),
                            mime='Biomedical/x-cif'
                        )
                    # Show GO terms
                    df_goterms = pd.DataFrame(columns=["position", "description", "go_terms"])
                    for id, item in enumerate(go_position_mapping0):
                        position = str(item['start']) + "-" + str(item['end'])
                        description = item['description']
                        go_terms = ", ".join([go['term'] for go in item['go_terms']])
                        df_goterms.loc[id] = [position, description, go_terms]

                    st.info("GO terms for the target sequence **{}**".format(target_seq.split("-")[0].upper()))
                    st.dataframe(df_goterms)

                    df_goterms = pd.DataFrame(columns=["position", "description", "go_terms"])
                    for id, item in enumerate(go_position_mapping1):
                        position = str(item['start']) + "-" + str(item['end'])
                        description = item['description']
                        go_terms = ", ".join([go['term'] for go in item['go_terms']])
                        df_goterms.loc[id] = [position, description, go_terms]

                    st.info("GO terms for the target sequence **{}**".format(target_seq.split("-")[1].upper()))
                    st.dataframe(df_goterms)
                
                st.markdown("**PPI Prediction**")

                model_name_col1 = df_path.loc[rowId, 'model_name']

                prediction_path_col1 = get_prediction_path("PPI", model_name_col1)

                result_col1 = None

                if prediction_path_col1 is not None:
                    result_col1 = get_prediction_result("PPI", prediction_path_col1, target_seq)

                if result_col1 is not None and len(result_col1) != 0:
                    st.info("The PPI prediction for target sequence **{}** in **{}** is **{}**.".format(target_seq, model_name_col1, result_col1.values[0]))
            else:
                col1, col2 = st.columns(2, gap="small")

                with col1:
                    pdb_path_col1 = df_path.loc[rowId, 'pdb_path']
                    model_name_col1 = df_path.loc[rowId, 'model_name']
                    st.info("{} structure".format(model_name_col1), icon="🧬")
                    viewer, pdb_content_col1 = show_pdb_with_disorder(pdb_path_col1)

                    buffer = BytesIO()
                    buffer.write(pdb_content_col1.encode('utf-8'))
                    buffer.seek(0)

                    if viewer is not None:
                        st.components.v1.html(viewer._make_html(), height=FIG_WIDTH, width=FIG_WIDTH)

                        if pdb_path_col1.endswith(".pdb"):
                            st.download_button(
                                label="📥 Download PDB file",
                                data=buffer,
                                file_name="{}_{}_structure.pdb".format(target_seq, model_name_col1),
                                mime='Biomedical/x-pdb'
                            )
                        else:
                            st.download_button(
                                label="📥 Download CIF file",
                                data=buffer,
                                file_name="{}_{}_structure.cif".format(target_seq, model_name_col1),
                                mime='Biomedical/x-cif'
                            )
                
                with col2:
                    pdb_path_col2 = df_path.loc[colIndex, 'pdb_path']
                    model_name_col2 = df_path.loc[colIndex, 'model_name']
                    st.info("{} structure".format(model_name_col2), icon="🧬")
                    viewer, pdb_content_col2 = show_pdb_with_disorder(pdb_path_col2)

                    buffer = BytesIO()
                    buffer.write(pdb_content_col2.encode('utf-8'))
                    buffer.seek(0)

                    if viewer is not None:
                        st.components.v1.html(viewer._make_html(), height=FIG_WIDTH, width=FIG_WIDTH)

                        if pdb_path_col2.endswith(".pdb"):
                            st.download_button(
                                label="📥 Download PDB file",
                                data=buffer,
                                file_name="{}_{}_structure.pdb".format(target_seq, model_name_col2),
                                mime='Biomedical/x-pdb'
                            )
                        else:
                            st.download_button(
                                label="📥 Download CIF file",
                                data=buffer,
                                file_name="{}_{}_structure.cif".format(target_seq, model_name_col2),
                                mime='Biomedical/x-cif'
                            )
                
                
                st.info("Aligned structure, {} in red, {} in blue".format(model_name_col1, model_name_col2), icon="🧬")
                viewer = show_pdb_align(pdb_path_col1, pdb_path_col2, width=FIG_WIDTH * 2)
                # Show structure
                st.components.v1.html(viewer._make_html(), height=FIG_WIDTH * 2, width=FIG_WIDTH * 2)

                df_goterms = pd.DataFrame(columns=["position", "description", "go_terms"])
                for id, item in enumerate(go_position_mapping0):
                    position = str(item['start']) + "-" + str(item['end'])
                    description = item['description']
                    go_terms = ", ".join([go['term'] for go in item['go_terms']])
                    df_goterms.loc[id] = [position, description, go_terms]

                st.info("GO terms for the target sequence **{}**".format(target_seq.split("-")[0].upper()))
                st.dataframe(df_goterms)

                df_goterms = pd.DataFrame(columns=["position", "description", "go_terms"])
                for id, item in enumerate(go_position_mapping1):
                    position = str(item['start']) + "-" + str(item['end'])
                    description = item['description']
                    go_terms = ", ".join([go['term'] for go in item['go_terms']])
                    df_goterms.loc[id] = [position, description, go_terms]

                st.info("GO terms for the target sequence **{}**".format(target_seq.split("-")[1].upper()))
                st.dataframe(df_goterms)

                model_name_col1 = df_path.loc[rowId, 'model_name']
                model_name_col2 = df_path.loc[colIndex, 'model_name']

                # st.markdown('<p style="font-size:24px;</p>', unsafe_allow_html=True)
                st.markdown("**PPI Prediction**")
                prediction_path_col1 = get_prediction_path("PPI", model_name_col1)
                prediction_path_col2 = get_prediction_path("PPI", model_name_col2)

                result_col1 = None
                result_col2 = None

                if prediction_path_col1 is not None:
                    result_col1 = get_prediction_result("PPI", prediction_path_col1, target_seq)

                if prediction_path_col2 is not None:
                    result_col2 = get_prediction_result("PPI", prediction_path_col2, target_seq)
                
                if result_col1 is not None and len(result_col1) != 0:
                    st.info("The PPI prediction for target sequence **{}** in **{}** is **{}**.".format(target_seq, model_name_col1, result_col1.values[0]))
                
                if result_col2 is not None and len(result_col2) != 0:
                    st.info("The PPI prediction for target sequence **{}** in **{}** is **{}**.".format(target_seq, model_name_col2, result_col2.values[0]))
            
            st.markdown("---")
            st.markdown("In PPI task, if the results is positive, it means that the two proteins are more likely to interacte with each other, otherwise, they are less likely to interact.")
            st.markdown("### 📚 Help Resources")
            st.markdown("""
                        - [Documentation](https://anonymous.4open.science/r/DisProtBench/)
                        - [Tutorial](https://anonymous.4open.science/r/DisProtBench/)
                        - [FAQ](https://anonymous.4open.science/r/DisProtBench/)
                        """)

def show_page_2():
    st.session_state.result = []
    models_list = ["AF2", "AF3", "Boltz", "Chai", "DeepFold", "ESMFold", "OmegaFold", "OpenFold", "Proteinx", "RoseTTAFold", "UniFold"] # temporally remove ESMFold

    global STOP_FLAG

    st.info("Enter protein ID and select models to compare.", icon="📝")
    col1, col2 = st.columns(2, gap="small")

    st.markdown("""
        <style>
        /* Outer div controlling the input box */
        div[data-baseweb="input"] {
            height: 400px; /* set height of the container */
        }

        /* The div inside that centers the input field */
        div[data-baseweb="input"] > div {
            height: 400px; /* match height */
            display: flex;
            align-items: flex-start;
            padding-top: 10px;
        }

        # /* The actual text input */
        # input[type="text"] {
        #     font-size: 24px !important;
        # }
        </style>
    """, unsafe_allow_html=True)

    with col1:
        st.markdown("⌨️**Enter protein id**", help="Please follow the protein form as the input sequence. \nExample: O43614")
        target_seq = st.text_input(label="Enter the protein sequence below", value="O43614").lower()

    progress_bar = st.progress(0)
    progress_text = st.empty()

    if not STOP_FLAG:
        with col2:
            st.markdown("📋**Select models to compare**", help="Please choose at least two models!")

            options_list = []
            for index, model in enumerate(models_list):
                options_list.append(st.checkbox(model, value=True))

            select_models_index = []
            for index, option in enumerate(options_list):
                progress_bar.progress(index)
                progress_text.text("{}% - Initializing...".format(index))
                if int(option) == 1:
                    select_models_index.append(index)

            models_selected_list = []
            for model_index in select_models_index:
                models_selected_list.append(models_list[model_index])

            STOP_FLAG = False
    
    if not STOP_FLAG:
        df_path = pd.DataFrame(columns=["model_name", "pdb_path"])
        for index, model_name in enumerate(models_selected_list):
            progress_bar.progress(15 + index)
            progress_text.text("{}% - Loading models...".format(str(15 + index)))
            
            tmp_path = get_id_path(target_seq, 'Drug', model_name)

            if tmp_path is None:
                st.error("File not found for {}, please deselect {}".format(model_name, model_name))
                STOP_FLAG = True
                st.stop()
            df_path.loc[df_path.shape[0]] = [model_name, tmp_path]

    if not STOP_FLAG:
        df_alignment_score = pd.DataFrame(columns=models_selected_list)
        
        p_i = 1
        for row, model_name_row in enumerate(models_selected_list):
            for col, model_name_col in enumerate(models_selected_list):
                progress_bar.progress(30 + int(p_i / 10))
                progress_text.text("{}% - Calculating RMSD...".format(30 + int(p_i / 10)))
                df_alignment_score.loc[row, model_name_col] = get_alignment_score(df_path.loc[row, "pdb_path"], df_path.loc[col, "pdb_path"])
                p_i += 1
    
    # st.write(st.session_state.result)

    progress_bar.progress(50)
    progress_text.text("{}% - Drawing heatmap...".format(50))
    
    if not STOP_FLAG:
        rowId = 0
        colIndex = 0

        if len(select_models_index) > 1:
            st.markdown("**Aligment RMSD Heatmap**", help="The RMSD comes from PyMOL package.")
            st.info("The numerical value in the heatmap indicates the RMSD between the two models. Click icon 🔍 before the numerical valu to indicate the model pairs to compare.")
            generate_rmsd_color_map()

            # models_selected_list = []
            # for model_index in select_models_index:
            #     models_selected_list.append(models_list[model_index])
            
            df_tmp = pd.DataFrame(columns=['Models'] + models_selected_list)

            p_i = 1
            for model_index, model_name in enumerate(models_selected_list):
                alignment_score_list = []
                for index, tmp_model_name in enumerate(models_selected_list):
                    progress_bar.progress(50 + int(p_i / 10))
                    progress_text.text("{}% - Loading RMSD values...".format(str(50 + int(p_i / 10))))
                    alignment_score_list.append(df_alignment_score.loc[model_index, tmp_model_name])
                    p_i += 1

                df_tmp.loc[df_tmp.shape[0]] = [model_name] + alignment_score_list

            response = make_aggrid(df_tmp, models_selected_list)

            progress_bar.progress(90)
            progress_text.text("{}% - Checking model pairs...".format(str(90)))

            df = response["data"]

            for col in df.columns:
                for idx in df.index:
                    if isinstance(df.at[idx, col], str):
                        try:
                            rowId = models_list.index(models_selected_list[int(idx)])
                            colIndex = models_list.index(col)
                        except Exception as e:
                            rowId = -1
                            colIndex = -1
                        
                        if df.at[idx, col].startswith("[clicked]"):
                            # st.write(rowId, colIndex)
                            if (rowId, colIndex) not in st.session_state.result:
                                st.session_state.result.append((rowId, colIndex))
                        else:
                            if (rowId, colIndex) in st.session_state.result:
                                st.session_state.result.remove((rowId, colIndex))

            # st.write(st.session_state.result)
            
            progress_bar.progress(100)
            progress_text.text("{}% - Complete!".format(str(100)))

            # Visulization
            if len(st.session_state.result) > 0:
                rowId = int(st.session_state.result[-1][0])
                colIndex = int(st.session_state.result[-1][1])
            else:
                rowId = 0
                colIndex = 0
            
            if (rowId > -1) & (colIndex > -1) & (rowId < len(models_selected_list)) & (colIndex < len(models_selected_list)):
                st.info("B-factor Color Mapping", icon="🌈")
                generate_color_map()  # Generate and display the color map
                if colIndex == rowId:
                    pdb_path_col1 = df_path.loc[rowId, 'pdb_path']
                    model_name_col1 = df_path.loc[rowId, 'model_name']
                    st.info("{}".format(model_name_col1), icon="🧬")
                    viewer, pdb_content_col1 = show_pdb_with_disorder(pdb_path_col1)

                    buffer = BytesIO()
                    buffer.write(pdb_content_col1.encode('utf-8'))
                    buffer.seek(0)

                    if viewer is not None:
                        st.components.v1.html(viewer._make_html(), height=FIG_WIDTH, width=FIG_WIDTH)

                        if pdb_path_col1.endswith(".pdb"):
                            st.download_button(
                                label="📥 Download PDB file",
                                data=buffer,
                                file_name="{}_{}_structure.pdb".format(target_seq, model_name_col1),
                                mime='Biomedical/x-pdb'
                            )
                        else:
                            st.download_button(
                                label="📥 Download CIF file",
                                data=buffer,
                                file_name="{}_{}_structure.cif".format(target_seq, model_name_col1),
                                mime='Biomedical/x-cif'
                            )
                    
                    st.markdown("**Drug Prediction**")

                    model_name_col1 = df_path.loc[rowId, 'model_name']

                    prediction_path_col1 = get_prediction_path("Drug", model_name_col1)

                    result_col1 = None

                    if prediction_path_col1 is not None:
                        result_col1 = get_prediction_result("Drug", prediction_path_col1, target_seq)

                    # st.markdown('<p style="font-size:24px;</p>', unsafe_allow_html=True)

                    if prediction_path_col1 is not None:
                        result_col1 = get_prediction_result("Drug", prediction_path_col1, target_seq)

                    if result_col1 is not None and len(result_col1) != 0:
                        st.info("The Drug prediction for target protein **{}** in **{}**:\n".format(target_seq.upper(), model_name_col1))
                        st.write(result_col1)
                else:
                    col1, col2 = st.columns(2, gap="small")
                    with col1:
                        pdb_path_col1 = df_path.loc[rowId, 'pdb_path']
                        model_name_col1 = df_path.loc[rowId, 'model_name']
                        st.info("{} structure".format(model_name_col1), icon="🧬")
                        viewer, pdb_content_col1 = show_pdb_with_disorder(pdb_path_col1)

                        buffer = BytesIO()
                        buffer.write(pdb_content_col1.encode('utf-8'))
                        buffer.seek(0)

                        if viewer is not None:
                            st.components.v1.html(viewer._make_html(), height=FIG_WIDTH, width=FIG_WIDTH)

                            if pdb_path_col1.endswith(".pdb"):
                                # st.download_button(label="📥 Download PDB File", data=pdb_content_col1, file_name="{}_{}_structure.pdb".format(target_seq, model_name_col1), mime="Biomedical/x-pdb")
                                st.download_button(
                                    label="📥 Download PDB file",
                                    data=buffer,
                                    file_name="{}_{}_structure.pdb".format(target_seq, model_name_col1),
                                    mime='Biomedical/x-pdb'
                                )
                            else:
                                st.download_button(
                                    label="📥 Download CIF file",
                                    data=buffer,
                                    file_name="{}_{}_structure.cif".format(target_seq, model_name_col1),
                                    mime='Biomedical/x-cif'
                                )
                    
                    with col2:
                        pdb_path_col2 = df_path.loc[colIndex, 'pdb_path']
                        model_name_col2 = df_path.loc[colIndex, 'model_name']
                        st.info("{} structure".format(model_name_col2), icon="🧬")
                        viewer, pdb_content_col2 = show_pdb_with_disorder(pdb_path_col2)

                        buffer = BytesIO()
                        buffer.write(pdb_content_col2.encode('utf-8'))
                        buffer.seek(0)

                        if viewer is not None:
                            st.components.v1.html(viewer._make_html(), height=FIG_WIDTH, width=FIG_WIDTH)

                            if pdb_path_col2.endswith(".pdb"):
                                # st.download_button(label="📥 Download PDB File", data=pdb_content_col1, file_name="{}_{}_structure.pdb".format(target_seq, model_name_col1), mime="Biomedical/x-pdb")
                                st.download_button(
                                    label="📥 Download PDB file",
                                    data=buffer,
                                    file_name="{}_{}_structure.pdb".format(target_seq, model_name_col2),
                                    mime='Biomedical/x-pdb'
                                )
                            else:
                                st.download_button(
                                    label="📥 Download CIF file",
                                    data=buffer,
                                    file_name="{}_{}_structure.cif".format(target_seq, model_name_col2),
                                    mime='Biomedical/x-cif'
                                )
                    
                    st.info("Aligned structure, {} in red, {} in blue".format(model_name_col1, model_name_col2), icon="🧬")
                    viewer = show_pdb_align(pdb_path_col1, pdb_path_col2, width=FIG_WIDTH * 2)
                    # Show structure
                    st.components.v1.html(viewer._make_html(), height=FIG_WIDTH * 2, width=FIG_WIDTH * 2)

                    st.markdown("**Drug Prediction**")

                    model_name_col1 = df_path.loc[rowId, 'model_name']
                    model_name_col2 = df_path.loc[colIndex, 'model_name']

                    prediction_path_col1 = get_prediction_path("Drug", model_name_col1)
                    prediction_path_col2 = get_prediction_path("Drug", model_name_col2)

                    result_col1 = None
                    result_col2 = None

                    if prediction_path_col1 is not None:
                        result_col1 = get_prediction_result("Drug", prediction_path_col1, target_seq)
                    
                    if prediction_path_col2 is not None:
                        result_col2 = get_prediction_result("Drug", prediction_path_col2, target_seq)

                    # st.markdown('<p style="font-size:24px;</p>', unsafe_allow_html=True)

                    if result_col1 is not None or len(result_col1) != 0:
                        st.info("The Drug prediction for target protein **{}** in **{}**:\n".format(target_seq.upper(), model_name_col1))
                        st.write(result_col1)
                    
                    if result_col2 is not None or len(result_col2) != 0:
                        st.info("The PPI prediction for target sequence **{}** in **{}**:\n".format(target_seq.upper(), model_name_col2))
                        st.write(result_col2)

                st.markdown("---")
                st.markdown("In Drug discover task, pKi stands for Protein Kinase Inhibitor. PKIs are small molecules or peptides that specifically inhibit the activity of protein kinases, which are enzymes that play a crucial role in cell signaling and many cellular processes. Targeting kinases is a major area of drug development, particularly in oncology, due to their involvement in various diseases.")
                st.markdown("### 📚 Help Resources")
                st.markdown("""
                            - [Documentation](https://anonymous.4open.science/r/DisProtBench/)
                            - [Tutorial](https://anonymous.4open.science/r/DisProtBench/)
                            - [FAQ](https://anonymous.4open.science/r/DisProtBench/)
                            """)

def show_page_3():
    global STOP_FLAG

    st.markdown("""
        <style>
        /* Outer div controlling the input box */
        div[data-baseweb="input"] {
            height: 150px; /* set height of the container */
        }

        /* The div inside that centers the input field */
        div[data-baseweb="input"] > div {
            height: 150px; /* match height */
            display: flex;
            align-items: flex-start;
            padding-top: 10px;
        }

        # /* The actual text input */
        # input[type="text"] {
        #     font-size: 24px !important;
        # }
        # </style>
    """, unsafe_allow_html=True)

    st.info("Enter protein sequence and select models to compare.", icon="📝")

    st.markdown("⌨️**Enter protein sequence**", help="Please follow the protein-protein form as the input sequence. \nExample: O60880-P10275")
    target_seq = st.text_input(label="Enter the protein sequence below", value="O60880-P10275")

    if "-" not in target_seq:
        st.error("🚫The input sequence format contains mistakes. Please check the input sequence!")
        STOP_FLAG = True

    progress_bar = st.progress(0)
    progress_text = st.empty()

    if not STOP_FLAG:
        full_seq_path, disorder_path, disorder_path2 = get_server_path(target_seq, 'Server')
        
        progress_bar.progress(66)
        progress_text.text("{}% - Loading RMSD values...".format(66))

        # Visulization
        st.info("B-factor Color Mapping", icon="🌈")
        generate_color_map()  # Generate and display the color map

        col1, col2, col3 = st.columns(3, gap="small")

        with col1:
            if full_seq_path is not None:
                st.info("{} structure".format(target_seq.upper()), icon="🧬")
                viewer, pdb_content_col1 = show_pdb_with_disorder(full_seq_path, width=int(FIG_WIDTH))

                buffer = BytesIO()
                buffer.write(pdb_content_col1.encode('utf-8'))
                buffer.seek(0)

                if viewer is not None:
                    st.components.v1.html(viewer._make_html(), height=int(FIG_WIDTH), width=int(FIG_WIDTH))

                    if full_seq_path.endswith(".pdb"):
                        st.download_button(
                            label="📥 Download PDB file",
                            data=buffer,
                            file_name="{}_structure.pdb".format(target_seq),
                            mime='Biomedical/x-pdb'
                        )
                    else:
                        st.download_button(
                            label="📥 Download CIF file",
                            data=buffer,
                            file_name="{}_structure.cif".format(target_seq),
                            mime='Biomedical/x-cif'
                        )
            else:
                st.info("Full sequence results for {} are not available now, please try again later!".format(target_seq.upper()))
        
        with col2:
            if disorder_path is not None:
                st.info("{} based on disordered {}".format(target_seq.upper(), target_seq.split("-")[0].upper()), icon="🧬")
                viewer, pdb_content_col2 = show_pdb_with_disorder(disorder_path, width=int(FIG_WIDTH))

                buffer = BytesIO()
                buffer.write(pdb_content_col2.encode('utf-8'))
                buffer.seek(0)

                if viewer is not None:
                    st.components.v1.html(viewer._make_html(), height=int(FIG_WIDTH), width=int(FIG_WIDTH))

                    if disorder_path.endswith(".pdb"):
                        # st.download_button(label="📥 Download PDB File", data=pdb_content_col1, file_name="{}_{}_structure.pdb".format(target_seq, model_name_col1), mime="Biomedical/x-pdb")
                        st.download_button(
                            label="📥 Download PDB file",
                            data=buffer,
                            file_name="{}_disorder_{}_structure.pdb".format(target_seq, target_seq.split("-")[0]),
                            mime='Biomedical/x-pdb'
                        )
                    else:
                        st.download_button(
                            label="📥 Download CIF file",
                            data=buffer,
                            file_name="{}_disorder_{}_structure.cif".format(target_seq, target_seq.split("-")[0]),
                            mime='Biomedical/x-cif'
                        )
            else:
                st.info("Disorder results for {} are not available now, please try again later!".format(target_seq.upper()))
        
        with col3:
            if disorder_path2 is not None:
                st.info("{} based on disordered {}".format(target_seq.upper(), target_seq.split("-")[1].upper()), icon="🧬")
                viewer, pdb_content_col3 = show_pdb_with_disorder(disorder_path2, width=int(FIG_WIDTH))

                buffer = BytesIO()
                buffer.write(pdb_content_col3.encode('utf-8'))
                buffer.seek(0)

                if viewer is not None:
                    st.components.v1.html(viewer._make_html(), height=int(FIG_WIDTH), width=int(FIG_WIDTH))

                    if disorder_path2.endswith(".pdb"):
                        # st.download_button(label="📥 Download PDB File", data=pdb_content_col1, file_name="{}_{}_structure.pdb".format(target_seq, model_name_col1), mime="Biomedical/x-pdb")
                        st.download_button(
                            label="📥 Download PDB file",
                            data=buffer,
                            file_name="{}_disorder_{}_structure.pdb".format(target_seq, target_seq.split("-")[1]),
                            mime='Biomedical/x-pdb'
                        )
                    else:
                        st.download_button(
                            label="📥 Download CIF file",
                            data=buffer,
                            file_name="{}_disorder_{}_structure.cif".format(target_seq, target_seq.split("-")[1]),
                            mime='Biomedical/x-cif'
                        )
            else:
                st.info("Disorder results for {} are not available now, please try again later!".format(target_seq.upper()))
        
        progress_bar.progress(100)
        progress_text.text("{}% - Complete!".format(str(100)))

def count_sessions():
    """Initialize or increment the session counter"""
    if 'visitor_count' not in st.session_state:
        st.session_state.visitor_count = 1
    else:
        st.session_state.visitor_count += 1
    
    # Store the last update time
    st.session_state.last_update = datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    
    return st.session_state.visitor_count, st.session_state.last_update 

def main():
    # Configure page
    st.set_page_config(page_title="ProteinBench", layout="wide")

    visitor_count, last_update = count_sessions()

    # Session state initialization
    if 'page' not in st.session_state:
        st.session_state.page = 'Page 1'

    if 'result' not in st.session_state:
        st.session_state.result = []

    st.markdown("""
        <style>
        .visitor-counter {
            background-color: #f0f2f6;
            padding: 1rem;
            border-radius: 0.5rem;
            margin-bottom: 1rem;
            text-align: center;
        }
        .visitor-count {
            font-size: 1.5rem;
            font-weight: bold;
            color: #1f77b4;
        }
        </style>
    """, unsafe_allow_html=True)

    # Navigation
    page_options = {
        "Page 1": "PPI",
        "Page 2": "Drug Discovery",
        "Page 3": "3D Structure",
    }
    
    # Add custom CSS for navigation styling
    st.markdown("""
        <style>
        .nav-container {
            display: flex;
            justify-content: center;
            margin-bottom: 2rem;
            padding: 1rem;
            background-color: #f0f2f6;
            border-radius: 0.5rem;
        }
        .nav-container select {
            font-size: 1.2rem;
            padding: 0.5rem 1rem;
            border-radius: 0.3rem;
            border: 1px solid #ccc;
        }
        .visitor-counter {
            position: fixed;
            bottom: 0;
            left: 0;
            right: 0;
            background-color: #f0f2f6;
            padding: 0.5rem;
            text-align: center;
            border-top: 1px solid #e0e0e0;
            z-index: 1000;
        }
        .visitor-count {
            font-size: 1.2rem;
            font-weight: bold;
            color: #1f77b4;
        }
        /* Add padding to the bottom of the main content to prevent overlap with fixed counter */
        .main .block-container {
            padding-bottom: 3rem;
        }
        </style>
    """, unsafe_allow_html=True)

    # st.markdown("""
    #     <div class="visitor-counter">
    #         <div>🕒 Last updated: {}  👥 Total Visitors: <span class="visitor-count">{}</span></div>
    #     </div>
    # """.format(last_update, visitor_count), unsafe_allow_html=True)
    
    # Display the main title
    st.title("DisProtBench")
    
    # Create a container for navigation
    selected_page = st.selectbox("**Navigation**", list(page_options.keys()), format_func=lambda x: page_options[x])

    # Show selected page
    if selected_page == "Page 1":
        show_page_1()
    elif selected_page =='Page 2':
        show_page_2()
    elif selected_page =='Page 3':
        show_page_3()


if __name__ == "__main__":
    main()