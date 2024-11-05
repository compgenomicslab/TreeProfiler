from bottle import route, run, request, redirect, template, static_file
import threading
from tempfile import NamedTemporaryFile
from collections import defaultdict
import os

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

@route('/upload', method='POST')
def do_upload():
    treename = request.forms.get('treename')
    tree_data = request.forms.get('tree')
    treeparser = request.forms.get('treeparser')
    metadata = request.forms.get('metadata')
    separator = request.forms.get('separator')
    text_prop = request.forms.get('text_prop',[])
    num_prop = request.forms.get('num_prop',[])
    bool_prop = request.forms.get('bool_prop',[])
    multiple_text_prop = request.forms.get('multiple_text_prop',[])
    alignment = request.forms.get('alignment')
    pfam = request.forms.get('pfam')
    print(text_prop, type(text_prop))
    column2method = {
        'alignment': 'none',
    }
    default_props = [
        'name', 'dist', 'support', 'rank', 'sci_name', 'taxid', 'lineage', 'named_lineage',
        'evoltype', 'dup_sp', 'dup_percent', 'lca'
    ]

    # Get file paths for uploaded files
    tree_file_path = uploaded_chunks.get(treename, {}).get("tree")
    metadata_file_path = uploaded_chunks.get(treename, {}).get("metadata")
    alignment_file_path = uploaded_chunks.get(treename, {}).get("alignment")
    pfam_file_path = uploaded_chunks.get(treename, {}).get("pfam")

    # Load tree from text input or file
    if tree_data:
        tree = utils.ete4_parse(tree_data, internal_parser=treeparser)
    elif tree_file_path and os.path.exists(tree_file_path):
        with open(tree_file_path, 'r') as f:
            tree = utils.ete4_parse(f.read(), internal_parser=treeparser)
        os.remove(tree_file_path)

    # Load metadata from text input or file (if available)
    metadata_bytes = None
    if metadata:
        metadata_bytes = metadata.encode('utf-8')
    elif metadata_file_path and os.path.exists(metadata_file_path):
        with open(metadata_file_path, 'rb') as f:
            metadata_bytes = f.read()
        os.remove(metadata_file_path)

    # Parse metadata if available
    metadata_options = {}
    if separator == "<tab>":
        separator = "\t"

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
            "text_prop": text_prop,
            "num_prop": num_prop,
            "bool_prop": bool_prop
        }
    
    # Group emapper-related arguments
    if pfam_file_path:
        emapper_options = {
            "emapper_mode": False,
            "emapper_pfam": pfam_file_path
        }
    else:
        emapper_options = {
            "emapper_mode": False,
            "emapper_pfam": None
        }
    # Annotate tree if metadata is provided
    annotated_tree, prop2type = run_tree_annotate(tree, 
    **metadata_options,
    **emapper_options, 
    alignment=alignment_file_path,
    column2method=column2method
    )
    annotated_newick = annotated_tree.write(props=None, format_root_node=True)
    
    # Cleanup temporary alignment file after usage
    if alignment_file_path and os.path.exists(alignment_file_path):
        os.remove(alignment_file_path)

    # Store the tree data
    node_props = metadata_options.get('node_props', [])
    node_props.extend(default_props)
    trees[treename] = {
        'tree': tree_data,
        'treeparser': treeparser,
        'metadata': metadata,
        'node_props': node_props,
        'annotated_tree': annotated_newick,
        'prop2type': prop2type,
        'layouts': [],
    }

    # Cleanup chunks after use
    if treename in uploaded_chunks:
        del uploaded_chunks[treename]

    # Redirect to tree page
    return redirect(f'/tree/{treename}')

@route('/tree/<treename>')
def show_tree(treename):
    # Fetch the tree data by its name
    tree_info = trees.get(treename)
    if tree_info:
        return template('tree_details', treename=treename, tree_info=tree_info)
    else:
        return f"Tree '{treename}' not found."

@route('/explore_tree/<treename>', method=['GET', 'POST'])
def explore_tree(treename):
    tree_info = trees.get(treename)
    if tree_info:
        # Load previously stored layouts, props, and metadata
        current_layouts = tree_info.get('layouts', [])
        current_layouts = []
        current_props = tree_info.get('node_props', list(tree_info['prop2type'].keys()))
        t = Tree(tree_info['annotated_tree'])
        
        if request.method == 'POST':
            selected_props = request.forms.getall('props') or current_props
            current_props.extend(selected_props)
            selected_layout = request.forms.get('layout')
            # Update layouts based on user selection
            column_width = 70
            level = 0
            padding_x = 1
            padding_y = 0
            color_config = {}

            # categorical
            if selected_layout == 'rectangle-layout':
                current_layouts, level, _ = tree_plot.get_rectangle_layouts(t, selected_props, level, tree_info['prop2type'])
            if selected_layout == 'label-layout':
                current_layouts, level, _ = tree_plot.get_label_layouts(t, selected_props, level, tree_info['prop2type'])
            if selected_layout == 'colorbranch-layout':
                current_layouts, level, _ = tree_plot.get_colorbranch_layouts(t, selected_props, level, tree_info['prop2type'])
            if selected_layout == 'categorical-bubble-layout':
                current_layouts, level, _ = tree_plot.get_categorical_bubble_layouts(t, selected_props, level, tree_info['prop2type'])
            if selected_layout == 'piechart-layout':
                current_layouts = tree_plot.get_piechart_layouts(t, selected_props, level, tree_info['prop2type'])
            if selected_layout == 'background-layout':
                current_layouts, level, _ = tree_plot.get_background_layouts(t, selected_props, level, tree_info['prop2type'])
            if selected_layout == 'profiling-layout':
                for profiling_prop in selected_props:
                    matrix, value2color, all_profiling_values = tree_plot.multiple2matrix(t, profiling_prop, prop2type=tree_info['prop2type'], color_config=color_config)
                    matrix_layout = tree_plot.profile_layouts.LayoutPropsMatrixBinary(name=f"Profiling_{profiling_prop}",
                    matrix=matrix, matrix_props=all_profiling_values, value_range=[0,1],
                    value_color=value2color, column=level, poswidth=column_width)
                    current_layouts.append(matrix_layout)
                    level += 1
            if selected_layout == 'categorical-matrix-layout':
                # drawing as array in matrix
                matrix, value2color = tree_plot.categorical2matrix(t, selected_props, color_config=color_config)
                matrix_layout = tree_plot.profile_layouts.LayoutPropsMatrixOld(name=f"Categorical_matrix_{selected_props}",
                matrix=matrix, matrix_type='categorical', matrix_props=selected_props,
                value_color=value2color, column=level, poswidth=column_width)
                current_layouts.append(matrix_layout)
                level += 1

            # binary
            if selected_layout == 'binary-layout':
                current_layouts, level, _ = tree_plot.get_binary_layouts(t, selected_props, 
                                            level, tree_info['prop2type'], column_width=column_width, 
                                            reverse=False, padding_x=padding_x, padding_y=padding_y,
                                            color_config=color_config, same_color=False, aggregate=False)
            if selected_layout == 'binary-aggregate-layout':
                current_layouts, level, _ = tree_plot.get_binary_layouts(t, selected_props, 
                                            level, tree_info['prop2type'], column_width=column_width, 
                                            reverse=False, padding_x=padding_x, padding_y=padding_y,
                                            color_config=color_config, same_color=False, aggregate=True)
            if selected_layout == 'binary-unicolor-layout':
                current_layouts, level, _ = tree_plot.get_binary_layouts(t, selected_props, 
                                            level, tree_info['prop2type'], column_width=column_width, 
                                            reverse=False, padding_x=padding_x, padding_y=padding_y,
                                            color_config=color_config, same_color=True, aggregate=False)
            if selected_layout == 'binary-unicolor-aggregate-layout':
                current_layouts, level, _ = tree_plot.get_binary_layouts(t, selected_props, 
                                            level, tree_info['prop2type'], column_width=column_width, 
                                            reverse=False, padding_x=padding_x, padding_y=padding_y,
                                            color_config=color_config, same_color=True, aggregate=True)
            
            # numerical
            internal_num_rep = 'avg'
            if selected_layout == 'branchscore-layout':
                current_layouts = tree_plot.get_branchscore_layouts(t, selected_props, level, 
                                    tree_info['prop2type'], internal_rep=internal_num_rep)
            if selected_layout == 'barplot-layout':
                current_layouts, level, _ = tree_plot.get_barplot_layouts(t, selected_props, level, 
                                    tree_info['prop2type'], internal_rep=internal_num_rep)
            if selected_layout == 'numerical-bubble-layout':
                current_layouts, level, _ = tree_plot.get_numerical_bubble_layouts(t, selected_props, level, 
                                    tree_info['prop2type'], internal_rep=internal_num_rep, color_config=color_config)
            if selected_layout == 'heatmap-layout':
                current_layouts, level = tree_plot.get_heatmap_layouts(t, selected_props, level,
                                    column_width=column_width, padding_x=padding_x, padding_y=padding_y, 
                                    internal_rep=internal_num_rep, color_config=color_config, norm_method='min-max',
                                    global_scaling=True)
            if selected_layout == 'heatmap-mean-layout':
                current_layouts, level = tree_plot.get_heatmap_layouts(t, selected_props, level,
                                    column_width=column_width, padding_x=padding_x, padding_y=padding_y, 
                                    internal_rep=internal_num_rep, color_config=color_config, norm_method='mean',
                                    global_scaling=True)
            if selected_layout == 'heatmap-zscore-layout':
                current_layouts, level = tree_plot.get_heatmap_layouts(t, selected_props, level,
                                    column_width=column_width, padding_x=padding_x, padding_y=padding_y, 
                                    internal_rep=internal_num_rep, color_config=color_config, norm_method='zscore',
                                    global_scaling=True)
            if selected_layout == 'numerical-matrix-layout':
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

            # analytic
            if selected_layout == 'acr-discrete-layout':
                current_layouts, level, _ = tree_plot.get_acr_discrete_layouts(t, selected_props, level, prop2type=tree_info['prop2type'],
                column_width=column_width, padding_x=padding_x, padding_y=padding_y, color_config=color_config) 
                #delta statistic 
                for suffix in ['delta', 'pval']:
                    current_props.extend([utils.add_suffix(prop, suffix) for prop in selected_props])
            
            if selected_layout == 'acr-continuous-layout':
                current_layouts, level, _ = tree_plot.get_acr_continuous_layouts(t, selected_props, level, prop2type=tree_info['prop2type'],
                column_width=column_width, padding_x=padding_x, padding_y=padding_y, color_config=color_config, continuous=True)
            
            if selected_layout == 'ls-layout':
                current_layouts, ls_props = tree_plot.get_ls_layouts(t, selected_props, level, prop2type=tree_info['prop2type'],
                column_width=column_width, padding_x=padding_x, padding_y=padding_y, color_config=color_config)
                current_props.extend(ls_props)
            
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

            if selected_layout == 'taxonclade-layout' or 'taxonrectangle-layout':
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
                current_layouts = current_layouts + taxa_layouts

            if selected_layout == 'alignment-layout':
                lengh = len(max(utils.tree_prop_array(t, 'alignment'),key=len))
                aln_layout = layouts.seq_layouts.LayoutAlignment(name='Alignment_layout', 
                        alignment_prop='alignment', column=level, scale_range=lengh,
                        summarize_inner_nodes=True)
                current_layouts.append(aln_layout)

            if selected_layout == 'domain-layout':
                domain_layout = layouts.seq_layouts.LayoutDomain(name="Domain_layout", prop='dom_arq')
                current_layouts.append(domain_layout)

            # Store updated props and layouts back to the tree_info
            tree_info['layouts'] = current_layouts
            tree_info['props'] = selected_props

            # Run the ETE tree explorer in a separate thread with updated props and layout
            def start_explore():
                t.explore(name=treename, layouts=current_layouts, port=5050, open_browser=False, include_props=current_props)

            explorer_thread = threading.Thread(target=start_explore)
            explorer_thread.start()
        else:
            def start_explore():
                t.explore(name=treename, layouts=current_layouts, port=5050, open_browser=False, include_props=current_props)
            explorer_thread = threading.Thread(target=start_explore)
            explorer_thread.start()
            
        return template('explore_tree', treename=treename, tree_info=tree_info, selected_props=current_props)

    return f"Tree '{treename}' not found."


run(host='localhost', port=8080)

