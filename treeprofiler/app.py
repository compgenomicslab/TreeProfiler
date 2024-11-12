from bottle import route, run, request, redirect, template, static_file, response 
import threading
from tempfile import NamedTemporaryFile
from collections import defaultdict
import os
import json
import time

from ete4 import Tree
from treeprofiler.tree_annotate import run_tree_annotate, parse_csv  # or other functions you need
from treeprofiler import tree_plot
from treeprofiler.src import utils
from treeprofiler import layouts


# In-memory storage for chunks and complete files
trees = {}
uploaded_chunks = {}
paired_color = [
    '#9a312f', '#9b57d0', '#f8ce9a', '#f16017', '#28fef9', '#53707a',
    '#213b07', '#b5e5ac', '#9640b2', '#a9bd10', '#69e42b', '#b44d67',
    '#b110c1', '#0b08a3', '#d07671', '#29e23b', '#3f2bf4', '#9b2a08',
    '#b42b94', '#77566a', '#2dfee7', '#046904', '#e2835d', '#53db2b',
    '#0b97e9', '#e0f6e9', '#ba46d1', '#4aba53', '#d4d6db', '#7a5d7c',
    '#4b100e', '#9e6373', '#5f4945', '#7e057a', '#f8e372', '#209f87',
    '#383f59', '#9d59e9', '#40c9fb', '#4cfc8b', '#d94769', '#20feba',
    '#c53238', '#068b02', '#6b4c93', '#f1968e', '#86d720', '#076fa6',
    '#0dbcfe', '#4d74b2', '#7b3dd2', '#286d26', '#a0faca', '#97505d',
    '#159e7a', '#fc05df', '#5df454', '#9160e1', '#c2eb5e', '#304fce',
    '#033379', '#54770f', '#271211', '#ab8479', '#37d9a0', '#f12205',
    '#cdd7e2', '#578f56', '#5ad9be', '#8596e9', '#c999ee', '#5f6b8a',
    '#f5c3a1', '#8e0603', '#cc21cf', '#65e7d0', '#97b3b6', '#d6220c',
    '#29c1e1', '#a30139', '#c9a619', '#a19410', '#da874f', '#64246d',
    '#66f35d', '#b8366c', '#116c95', '#bd851a', '#27f7cb', '#512ca4',
    '#60e72e', '#d1941c', '#1045a8', '#c1b03a', '#0c62a5', '#7ac9b2',
    '#6bb9bd', '#cb30eb', '#26bad0', '#d9e557'
]


job_status = {}  # Dictionary to store job statuses
@route('/upload', method='POST')
def do_upload():
    treename = request.forms.get('treename')
    if not treename:
        response.status = 400
        return {"error": "Tree name is required"}

    if treename in job_status:
        response.status = 400
        return {"error": "A job with this tree name is already in progress"}

    job_status[treename] = "running"

    # Collect all form data
    job_args = {
        "treename": treename,
        "tree_data": request.forms.get('tree'),
        "treeparser": request.forms.get('treeparser'),
        "metadata": request.forms.get('metadata'),
        "separator": request.forms.get('separator'),
        "text_prop": request.forms.getlist('text_prop[]'),
        "num_prop": request.forms.getlist('num_prop[]'),
        "bool_prop": request.forms.getlist('bool_prop[]'),
        "multiple_text_prop": request.forms.getlist('multiple_text_prop[]'),
        "taxon_column": request.forms.get('taxonomicIdColumn'),
        "taxadb": request.forms.get('taxaDb'),
        "species_delimiter": request.forms.get('speciesFieldDelimiter'),
        "species_index": int(request.forms.get('speciesFieldIndex') or 0),
        "version": request.forms.get('version'),
        "ignore_unclassified": request.forms.get('ignoreUnclassified') == 'True',
        
        "acr_columns": request.forms.getlist('acrColumn[]'),
        "prediction_method": request.forms.get('acrMethod'),
        "model": request.forms.get('acrModel'),
        "phylosignal": request.forms.get('phylosignal'),
        "ent_type": request.forms.get('entType'),
        "iteration": request.forms.get('iteration'),
        "lambda0": request.forms.get('lambda0'),
        "se": request.forms.get('se'),
        "thin": request.forms.get('thin'),
        "burn": request.forms.get('burn'),

        "ls_columns": request.forms.getlist('lsColumn[]'),
        "prec_cutoff": request.forms.get('percision'),
        "sens_cutoff": request.forms.get('sensitivity'),

        "alignment": request.forms.get('alignment'),
        "consensus_cutoff": request.forms.get('consensusCutoff'),

        "pfam": request.forms.get('pfam'),
        "summary_methods": request.forms.get('summary_methods')
    }

    # Start the upload job in a new thread
    threading.Thread(target=process_upload_job, args=(job_args,)).start()
    response.content_type = 'application/json'
    return json.dumps({"job_id": treename})

def process_upload_job(job_args):
    """Function to handle the actual processing of the uploaded data."""
    treename = job_args['treename']
    # Initialize variables and parse arguments
    separator = "\t" if job_args.get("separator") == "<tab>" else job_args.get("separator", ",")
    tree_data = job_args.get("tree_data")
    treeparser = job_args.get("treeparser")
    summary_methods = job_args.get("summary_methods")
    
    # Parse summary methods JSON, if provided
    column2method = json.loads(summary_methods) if summary_methods else {}
    column2method['alignment'] = 'none'
    
    # Retrieve paths for uploaded files
    tree_file_path = uploaded_chunks.get(treename, {}).get("tree")
    metadata_file_path = uploaded_chunks.get(treename, {}).get("metadata")
    alignment_file_path = uploaded_chunks.get(treename, {}).get("alignment")
    pfam_file_path = uploaded_chunks.get(treename, {}).get("pfam")
    
    # Load the tree data
    if tree_data:
        tree = utils.ete4_parse(tree_data, internal_parser=treeparser)
    elif tree_file_path and os.path.exists(tree_file_path):
        with open(tree_file_path, 'r') as f:
            tree = utils.ete4_parse(f.read(), internal_parser=treeparser)
        os.remove(tree_file_path)

    # Prepare metadata processing if available
    metadata_bytes = job_args.get("metadata").encode('utf-8') if job_args.get("metadata") else None
    if metadata_file_path and os.path.exists(metadata_file_path):
        with open(metadata_file_path, 'rb') as f:
            metadata_bytes = f.read()
        os.remove(metadata_file_path)

    # Process metadata
    metadata_options = {}
    if metadata_bytes:
        with NamedTemporaryFile(suffix='.tsv') as f_annotation:
            f_annotation.write(metadata_bytes)
            f_annotation.flush()
            metadata_dict, node_props, columns, prop2type = parse_csv([f_annotation.name], delimiter=separator)

        metadata_options = {
            "metadata_dict": metadata_dict,
            "node_props": node_props,
            "columns": columns,
            "prop2type": prop2type,
            "text_prop": job_args.get("text_prop"),
            "num_prop": job_args.get("num_prop"),
            "bool_prop": job_args.get("bool_prop"),
            "multiple_text_prop": job_args.get("multiple_text_prop")
        }

    # Taxonomic annotation options
    taxonomic_options = {}
    if job_args.get("taxon_column"):
        taxonomic_options = {
            "taxon_column": job_args.get("taxon_column"),
            "taxadb": job_args.get("taxadb"),
            "gtdb_version": job_args.get("version"),
            "taxon_delimiter": job_args.get("species_delimiter"),
            "taxa_field": job_args.get("species_index"),
            "ignore_unclassified": job_args.get("ignore_unclassified")
        }

    # for analytic methods
    analytic_options = {}
    acr_discrete_columns = []
    acr_continuous_columns = []
    if job_args.get("acr_columns"):
        for acr_prop in job_args.get("acr_columns"):
            if prop2type.get(acr_prop) == str:
                acr_discrete_columns.append(acr_prop)
            elif prop2type.get(acr_prop) == float:
                acr_continuous_columns.append(acr_prop)
        analytic_options = {
            "acr_discrete_columns": acr_discrete_columns,
            "acr_continuous_columns": acr_continuous_columns,
            "prediction_method": job_args.get("prediction_method"),
            "model": job_args.get("model")
        }
        if job_args.get("phylosignal") == 'delta':
            analytic_options['delta_stats'] = True
            analytic_options['ent_type'] = job_args.get("ent_type")
            analytic_options['iteration'] = int(job_args.get("iteration"))
            analytic_options['lambda0'] = float(job_args.get("lambda0"))
            analytic_options['se'] = float(job_args.get("se"))
            analytic_options['thin'] = int(job_args.get("thin"))
            analytic_options['burn'] = int(job_args.get("burn"))

    ls_columns = []
    if job_args.get("ls_columns"):
        ls_columns = job_args.get("ls_columns")
        analytic_options['ls_columns'] = ls_columns
        analytic_options['prec_cutoff'] = float(job_args.get("prec_cutoff"))
        analytic_options['sens_cutoff'] = float(job_args.get("sens_cutoff"))

    # Alignment options
    alignment_options = {}
    if alignment_file_path:
        alignment_options = {
            "alignment": alignment_file_path,
            "consensus_cutoff": float(job_args.get("consensus_cutoff"))
        }

    # Emapper options
    emapper_options = {
        "emapper_mode": False,
        "emapper_pfam": pfam_file_path if pfam_file_path else None
    }

    # Run annotation
    annotated_tree, prop2type = run_tree_annotate(
        tree,
        **metadata_options,
        **emapper_options,
        **taxonomic_options,
        **analytic_options,
        **alignment_options,
        column2method=column2method
    )

    # Post-processing of annotated tree properties
    list_keys = [key for key, value in prop2type.items() if value == list]
    for node in annotated_tree.leaves():
        for key in list_keys:
            if node.props.get(key):
                node.add_prop(key, '||'.join(node.props[key]))

    avail_props = [key for key in prop2type.keys() if key not in ['name', 'dist', 'support']]
    annotated_newick = annotated_tree.write(props=avail_props, format_root_node=True)

    # Store the processed data
    trees[treename] = {
        'tree': tree_data,
        'treeparser': treeparser,
        'metadata': job_args.get("metadata"),
        'node_props': metadata_options.get('node_props', []) + [
            'rank', 'sci_name', 'taxid', 'evoltype', 'dup_sp', 'dup_percent'
        ],
        'annotated_tree': annotated_newick,
        'prop2type': prop2type,
        'layouts': []
    }

    # Cleanup temporary alignment file after usage
    if alignment_file_path and os.path.exists(alignment_file_path):
        os.remove(alignment_file_path)

    # Mark job as complete
    job_status[treename] = "complete"
        
    # except Exception as e:
    #     job_status[treename] = "failed"
    #     print(f"Error processing job {treename}: {e}")
    # finally:
    #     # Remove stored chunks for the completed job
    #     if treename in uploaded_chunks:
    #         del uploaded_chunks[treename]


# Route to serve static files like CSS and JS
@route('/static/<filepath:path>')
def server_static(filepath):
    return static_file(filepath, root='./static')

@route('/')
def upload_tree():
    return template('upload_tree')

@route('/upload_chunk', method='POST')
def upload_chunk():
    # Determine the file type based on the request
    chunk = request.files.get('treeFile') or request.files.get('metadataFile') or request.files.get('alignmentFile') or request.files.get('pfamFile')
    chunk_index = int(request.forms.get("chunkIndex"))
    total_chunks = int(request.forms.get("totalChunks"))
    treename = request.forms.get("treename")

    # Identify the file type
    if 'treeFile' in request.files:
        file_type = "tree"
    elif 'metadataFile' in request.files:
        file_type = "metadata"
    elif 'alignmentFile' in request.files:
        file_type = "alignment"
    elif 'pfamFile' in request.files:
        file_type = "pfam"
    else:
        return "Unknown file type", 400

    # Initialize chunk storage for each file type if it doesn't exist
    if treename not in uploaded_chunks:
        uploaded_chunks[treename] = {
            "tree": [], 
            "metadata": [], 
            "alignment": [], 
            "pfam": []
        }
    
    # Append the chunk to the corresponding file list
    uploaded_chunks[treename][file_type].append((chunk_index, chunk.file.read()))

    # Assemble the file when all chunks are received
    if len(uploaded_chunks[treename][file_type]) == total_chunks:
        uploaded_chunks[treename][file_type].sort()  # Ensure chunks are in order
        with NamedTemporaryFile(delete=False) as assembled_file:
            for _, part in uploaded_chunks[treename][file_type]:
                assembled_file.write(part)
            assembled_file_path = assembled_file.name
        uploaded_chunks[treename][file_type] = assembled_file_path  # Save assembled file path

    return "Chunk received"

@route('/job_status/<job_id>')
def job_status_check(job_id):
    status = job_status.get(job_id, "not_found")
    response.content_type = 'application/json'
    return {"status": status}


@route('/show_tree/<treename>')
def show_tree(treename):
    # Fetch and display the processed tree, if available
    tree_info = trees.get(treename)
    if tree_info:
        return f"Tree '{treename}' processed successfully! Displaying results."
    return "Processing still in progress or failed."

@route('/tree/<treename>/result')
def show_tree(treename):
    # Fetch the tree data by its name
    tree_info = trees.get(treename)
    if tree_info:
        return template('tree_details', treename=treename, tree_info=tree_info)
    else:
        return f"Tree '{treename}' not found."

@route('/tree/<treename>')
def tree_status(treename):
    """Serves the job_running.html template to show job status."""
    return template('job_running', treename=treename, job_id=treename)

@route('/check_job_status')
def check_job_status():
    """API endpoint to check the status of a given job."""
    job_id = request.query.get('job_id')
    status = job_status.get(job_id, "not_found")
    return status

@route('/explore_tree/<treename>', method=['GET', 'POST'])
def explore_tree(treename):
    tree_info = trees.get(treename)
    if not tree_info:
        return f"Tree '{treename}' not found."
    
    # Retrieve and initialize layouts, properties, and tree
    current_layouts = tree_info.get('layouts', [])
    current_props = list(tree_info['prop2type'].keys())
    t = Tree(tree_info['annotated_tree'])

    # Default values for layout settings (can be passed from form in the future)
    level = 0
    column_width = 70
    padding_x, padding_y = 1, 0
    color_config = {}
    internal_num_rep = 'avg'

    # Process POST request
    if request.method == 'POST':
        layers_data = request.forms.get('layers')
        if layers_data:
            layers = json.loads(layers_data)
            for layer in layers:
                # Process each layer individually
                current_layouts, current_props, level = process_layer(
                    t, layer, tree_info, current_layouts, current_props, level,
                    column_width, padding_x, padding_y, color_config, internal_num_rep, paired_color
                )
    
    # Start the ete exploration thread
    start_explore_thread(t, treename, current_layouts, current_props)
    
    # Render template
    return template('explore_tree', treename=treename, tree_info=tree_info, selected_props=current_props)

def process_layer(t, layer, tree_info, current_layouts, current_props, level, column_width, padding_x, padding_y, color_config, internal_num_rep, paired_color):
    """
    Processes a single layer, applying user selections and updating layouts, properties, and level.
    """
    selected_props = layer.get('props', [])
    selected_layout = layer.get('layout', '')
    query_type = layer.get('queryType', '')
    query_box = layer.get('query', '')

    current_props.extend(selected_props)

    # Process layout selection for the current layer
    level, current_layouts = apply_layouts(
        t, selected_layout, selected_props, tree_info, current_layouts, 
        level, column_width, padding_x, padding_y, color_config, internal_num_rep
    )

    # Process query actions for the current layer
    current_layouts = apply_queries(query_type, query_box, current_layouts, paired_color, tree_info, level)

    return current_layouts, current_props, level

def apply_layouts(t, selected_layout, selected_props, tree_info, current_layouts, level, column_width, padding_x, padding_y, color_config, internal_num_rep):
    """
    Applies the selected layout to the current tree and updates the layouts list.
    """
    # Additional layout types (Binary, Numerical, Analytic, Taxonomy, etc.)
    # can be added here following the same pattern
    # Categorical layouts
    if selected_layout == 'rectangle-layout':
        rectangle_layouts, level, _ = tree_plot.get_rectangle_layouts(t, selected_props, level, tree_info['prop2type'])
        current_layouts.extend(rectangle_layouts)
    elif selected_layout == 'label-layout':
        label_layouts, level, _ = tree_plot.get_label_layouts(t, selected_props, level, tree_info['prop2type'])
        current_layouts.extend(label_layouts)
    elif selected_layout == 'colorbranch-layout':
        colorbranch_layouts, level, _ = tree_plot.get_colorbranch_layouts(t, selected_props, level, tree_info['prop2type'])
        current_layouts.extend(colorbranch_layouts)
    elif selected_layout == 'category-bubble-layout':
        bubble_layouts, level, _ = tree_plot.get_categorical_bubble_layouts(t, selected_props, level, tree_info['prop2type'])
        current_layouts.extend(bubble_layouts)
    elif selected_layout == 'piechart-layout':
        piechart_layouts = tree_plot.get_piechart_layouts(t, selected_props, level, tree_info['prop2type'])
        current_layouts.extend(piechart_layouts)
    elif selected_layout == 'background-layout':
        background_layouts, level, _ = tree_plot.get_background_layouts(t, selected_props, level, tree_info['prop2type'])
        current_layouts.extend(background_layouts)
    elif selected_layout == 'profiling-layout':
        for profiling_prop in selected_props:
            matrix, value2color, all_profiling_values = tree_plot.multiple2matrix(t, profiling_prop, prop2type=tree_info['prop2type'], color_config=color_config)
            matrix_layout = tree_plot.profile_layouts.LayoutPropsMatrixBinary(name=f"Profiling_{profiling_prop}",
            matrix=matrix, matrix_props=all_profiling_values, value_range=[0,1],
            value_color=value2color, column=level, poswidth=column_width)
            current_layouts.append(matrix_layout)
            level += 1
    elif selected_layout == 'categorical-matrix-layout':
        # drawing as array in matrix
        matrix, value2color = tree_plot.categorical2matrix(t, selected_props, color_config=color_config)
        matrix_layout = tree_plot.profile_layouts.LayoutPropsMatrixOld(name=f"Categorical_matrix_{selected_props}",
        matrix=matrix, matrix_type='categorical', matrix_props=selected_props,
        value_color=value2color, column=level, poswidth=column_width)
        current_layouts.append(matrix_layout)
        level += 1

    # Binary layouts
    if selected_layout == 'binary-layout':
        binary_layouts, level, _ = tree_plot.get_binary_layouts(t, selected_props, 
            level, tree_info['prop2type'], column_width=column_width, 
            reverse=False, padding_x=padding_x, padding_y=padding_y,
            color_config=color_config, same_color=False, aggregate=False)
        current_layouts.extend(binary_layouts)
    elif selected_layout == 'binary-aggregate-layout':
        binary_layouts, level, _ = tree_plot.get_binary_layouts(t, selected_props, 
            level, tree_info['prop2type'], column_width=column_width, 
            reverse=False, padding_x=padding_x, padding_y=padding_y,
            color_config=color_config, same_color=False, aggregate=True)
        current_layouts.extend(binary_layouts)
    elif selected_layout == 'binary-unicolor-layout':
        binary_layouts, level, _ = tree_plot.get_binary_layouts(t, selected_props, 
            level, tree_info['prop2type'], column_width=column_width, 
            reverse=False, padding_x=padding_x, padding_y=padding_y,
            color_config=color_config, same_color=True, aggregate=False)
        current_layouts.extend(binary_layouts)
    elif selected_layout == 'binary-unicolor-aggregate-layout':
        binary_layouts, level, _ = tree_plot.get_binary_layouts(t, selected_props, 
            level, tree_info['prop2type'], column_width=column_width, 
            reverse=False, padding_x=padding_x, padding_y=padding_y,
            color_config=color_config, same_color=True, aggregate=True)
        current_layouts.extend(binary_layouts)

    # Numerical layouts
    if selected_layout == 'branchscore-layout':
        branchscore_layouts = tree_plot.get_branchscore_layouts(t, selected_props, level, 
            tree_info['prop2type'], internal_rep=internal_num_rep)
        current_layouts.extend(branchscore_layouts)
    elif selected_layout == 'barplot-layout':
        barplot_layouts, level, _ = tree_plot.get_barplot_layouts(t, selected_props, level, 
            tree_info['prop2type'], internal_rep=internal_num_rep)
        current_layouts.extend(barplot_layouts)
    elif selected_layout == 'numerical-bubble-layout':
        bubble_layouts_layouts, level, _ = tree_plot.get_numerical_bubble_layouts(t, selected_props, level, 
            tree_info['prop2type'], internal_rep=internal_num_rep, color_config=color_config)
        current_layouts.extend(bubble_layouts_layouts)
    elif selected_layout == 'heatmap-layout':
        heatmap_layouts, level = tree_plot.get_heatmap_layouts(t, selected_props, level,
            column_width=column_width, padding_x=padding_x, padding_y=padding_y, 
            internal_rep=internal_num_rep, color_config=color_config, norm_method='min-max',
            global_scaling=True)
        current_layouts.extend(heatmap_layouts)
    elif selected_layout == 'heatmap-mean-layout':
        heatmap_layouts, level = tree_plot.get_heatmap_layouts(t, selected_props, level,
            column_width=column_width, padding_x=padding_x, padding_y=padding_y, 
            internal_rep=internal_num_rep, color_config=color_config, norm_method='mean',
            global_scaling=True)
        current_layouts.extend(heatmap_layouts)
    elif selected_layout == 'heatmap-zscore-layout':
        heatmap_layouts, level = tree_plot.get_heatmap_layouts(t, selected_props, level,
            column_width=column_width, padding_x=padding_x, padding_y=padding_y, 
            internal_rep=internal_num_rep, color_config=color_config, norm_method='zscore',
            global_scaling=True)
        current_layouts.extend(heatmap_layouts)
    elif selected_layout == 'numerical-matrix-layout':
        matrix, minval, maxval, value2color, results_list, list_props, single_props = tree_plot.numerical2matrix(t, 
        selected_props, count_negative=True, internal_num_rep=internal_num_rep, 
        color_config=color_config, norm_method='min-max')
        if list_props:
            index_map = {value: idx for idx, value in enumerate(selected_props)}
            sorted_list_props = sorted(list_props, key=lambda x: index_map[x])
            for list_prop in sorted_list_props:
                matrix, minval, maxval, value2color = results_list[list_prop]
                matrix_layout = tree_plot.profile_layouts.LayoutPropsMatrixOld(name=f"Numerical_matrix_{list_prop}", 
                    matrix=matrix, matrix_type='numerical', matrix_props=[list_prop], is_list=True, 
                    value_color=value2color, value_range=[minval, maxval], column=level,
                    poswidth=column_width)
                level += 1
                current_layouts.append(matrix_layout)
        if single_props:
            index_map = {value: idx for idx, value in enumerate(selected_props)}
            sorted_single_props = sorted(single_props, key=lambda x: index_map[x])
            matrix_layout = tree_plot.profile_layouts.LayoutPropsMatrixOld(name=f"Numerical_matrix_{sorted_single_props}", 
                matrix=matrix, matrix_type='numerical', matrix_props=sorted_single_props, is_list=False, 
                value_color=value2color, value_range=[minval, maxval], column=level,
                poswidth=column_width)
            level += 1
            current_layouts.append(matrix_layout)

    # Analytic layouts
    if selected_layout == 'acr-discrete-layout':
        acr_layouts, level, _ = tree_plot.get_acr_discrete_layouts(t, selected_props, level, prop2type=tree_info['prop2type'],
        column_width=column_width, padding_x=padding_x, padding_y=padding_y, color_config=color_config) 
        #delta statistic 
        for suffix in ['delta', 'pval']:
            current_props.extend([utils.add_suffix(prop, suffix) for prop in selected_props])
        current_layouts.extend(acr_layouts)

    if selected_layout == 'acr-continuous-layout':
        acr_layouts = tree_plot.get_acr_continuous_layouts(t, selected_props, level, prop2type=tree_info['prop2type'],
        padding_x=padding_x, padding_y=padding_y)
        current_layouts.extend(acr_layouts)

    if selected_layout == 'ls-layout':
        ls_layouts, ls_props = tree_plot.get_ls_layouts(t, selected_props, level, prop2type=tree_info['prop2type'],
        padding_x=padding_x, padding_y=padding_y, color_config=color_config)
        current_props.extend(ls_props)
        current_layouts.extend(ls_layouts)

    # Taxonomy layouts
    if selected_layout == 'taxoncollapse-layout':
        taxon_color_dict = {}
        taxa_layouts = []
        # generate a rank2values dict for pre taxonomic annotated tree
        rank2values = defaultdict(list)
        for n in t.traverse():
            if n.props.get('lca'):
                lca_dict = utils.string_to_dict(n.props.get('lca'))
                for rank, sci_name in lca_dict.items():
                    rank2values[rank].append(sci_name)

            current_rank = n.props.get('rank')
            if current_rank and current_rank != 'Unknown':
                rank2values[current_rank].append(n.props.get('sci_name',''))

        # assign color for each value of each rank
        for rank, value in sorted(rank2values.items()):
            value = list(set(value))
            color_dict = utils.assign_color_to_values(value, paired_color)
            taxa_layout = layouts.taxon_layouts.TaxaCollapse(name = "TaxaCollapse_"+rank, rank=rank, 
            rect_width=column_width, color_dict=color_dict, column=level)
            taxa_layouts.append(taxa_layout)
            
        current_layouts = current_layouts + taxa_layouts
        level += 1

    if selected_layout == 'taxonclade-layout' or selected_layout == 'taxonrectangle-layout':
        taxon_color_dict = {}
        taxa_layouts = []
        rank2values = defaultdict(list)
        for n in t.traverse():
            if n.props.get('rank') and n.props.get('rank') != 'Unknown':
                rank = n.props.get('rank')
                rank2values[rank].append(n.props.get('sci_name',''))

        # assign color for each value of each rank
        for rank, value in sorted(rank2values.items()):
            value = list(set(value))
            color_dict = utils.assign_color_to_values(value, paired_color)
            if selected_layout == 'taxonclade-layout':
                taxa_layout = layouts.taxon_layouts.TaxaClade(name='TaxaClade_'+rank, level=level, rank = rank, color_dict=color_dict)
                taxa_layouts.append(taxa_layout)

            if selected_layout == 'taxonrectangle-layout':
                taxa_layout = layouts.taxon_layouts.TaxaRectangular(name = "TaxaRect_"+rank, rank=rank, rect_width=column_width, color_dict=color_dict, column=level)
                taxa_layouts.append(taxa_layout)
                #level += 1
            taxon_color_dict[rank] = color_dict

        #taxa_layouts.append(taxon_layouts.TaxaRectangular(name = "Last Common Ancester", color_dict=taxon_color_dict, column=level))
        taxa_layouts.append(layouts.taxon_layouts.LayoutSciName(name = 'Taxa Scientific name', color_dict=taxon_color_dict))
        taxa_layouts.append(layouts.taxon_layouts.LayoutEvolEvents(name='Taxa Evolutionary events', prop="evoltype",
            speciation_color="blue", 
            duplication_color="red", node_size = 3,
            legend=True))
        current_layouts.extend(taxa_layouts)
    
    # Alignment layout
    if selected_layout == 'alignment-layout':
        lengh = len(max(utils.tree_prop_array(t, 'alignment'),key=len))
        window = []
        aln_layout = layouts.seq_layouts.LayoutAlignment(name='Alignment_layout', 
                alignment_prop='alignment', column=level, scale_range=lengh, window=window,
                summarize_inner_nodes=True)
        current_layouts.append(aln_layout)
    
    # Domain layout
    if selected_layout == 'domain-layout':
        domain_layout = layouts.seq_layouts.LayoutDomain(name="Domain_layout", prop='dom_arq')
        current_layouts.append(domain_layout)
    
    return level, current_layouts

def apply_queries(query_type, query_box, current_layouts, paired_color, tree_info, level):
    """
    Applies conditional queries to the tree based on the provided query type and query box.
    """
    if query_type and query_box:
        query_strings = convert_query_string(query_box)  # convert query to list format

        if query_type == 'highlight':
            current_layouts = apply_highlight_queries(query_strings, current_layouts, paired_color, tree_info, level)
        
        elif query_type == 'collapse':
            current_layouts = apply_collapse_queries(query_strings, current_layouts, paired_color, tree_info, level)

    return current_layouts


def apply_highlight_queries(query_strings, current_layouts, paired_color, tree_info, level):
    """
    Processes highlight queries and appends them to current layouts.
    """
    for idx, condition in enumerate(query_strings):
        syntax_sep = ','
        condition_list = condition.split(syntax_sep)
        color2conditions = {paired_color[idx]: condition_list}
        
        s_layout = layouts.conditional_layouts.LayoutHighlight(
            name=f'Highlight_{condition}', color2conditions=color2conditions, 
            column=level, prop2type=tree_info['prop2type']
        )
        current_layouts.append(s_layout)
    
    return current_layouts


def apply_collapse_queries(query_strings, current_layouts, paired_color, tree_info, level):
    """
    Processes collapse queries and appends them to current layouts.
    """
    for idx, condition in enumerate(query_strings):
        syntax_sep = ','
        condition_list = condition.split(syntax_sep)
        color2conditions = {paired_color[idx]: condition_list}
        
        c_layout = layouts.conditional_layouts.LayoutCollapse(
            name=f'Collapse_{condition}', color2conditions=color2conditions, 
            column=level, prop2type=tree_info['prop2type']
        )
        current_layouts.append(c_layout)
    
    return current_layouts

def start_explore_thread(t, treename, current_layouts, current_props):
    """
    Starts the ete exploration in a separate thread.
    """

    def explore():
        t.explore(name=treename, layouts=current_layouts, port=5050, open_browser=False, include_props=current_props)
    
    explorer_thread = threading.Thread(target=explore)
    explorer_thread.start()

def convert_query_string(query_string):
    # Split the query by semicolon to get each statement
    queries = [q.strip() for q in query_string.split(';') if q.strip()]
    result = []
    for query in queries:
        if " OR " in query:
            # Split the query by 'OR' and convert each condition
            conditions = [cond.strip() for cond in query.split(" OR ")]
            result.extend(conditions)  # Add each condition as a separate entry for 'OR'
        elif " AND " in query:
            # Split the query by 'AND' and convert each condition
            conditions = [cond.strip() for cond in query.split(" AND ")]
            result.append(','.join(conditions))  # Join conditions with a comma for 'AND'
        else:
            result.append(query)
    return result

run(host='localhost', port=8080)

