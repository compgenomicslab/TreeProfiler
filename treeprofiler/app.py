from bottle import route, run, request, redirect, template, static_file
import threading
from tempfile import NamedTemporaryFile

from ete4 import Tree
from treeprofiler.tree_annotate import run_tree_annotate, parse_csv  # or other functions you need
from treeprofiler import tree_plot

from treeprofiler.src import utils
# In-memory storage for simplicity (you can replace this with a database or file storage)
trees = {}

# Route to serve static files like CSS and JS
@route('/static/<filepath:path>')
def server_static(filepath):
    return static_file(filepath, root='./static')

@route('/')
def upload_tree():
    return template('upload_tree')

@route('/upload', method='POST')
def do_upload():    
    treename = request.forms.get('treename')
    tree_data = request.forms.get('tree')
    metadata = request.forms.get('metadata')

    # here start the tree annotation process
    # load the tree
    tree = utils.ete4_parse(tree_data)
    # parse the metadata
    metadata_bytes = metadata.encode('utf-8')  # Convert string to bytes
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(metadata_bytes)
        f_annotation.flush()
        metadata_dict, node_props, columns, prop2type = parse_csv([f_annotation.name],delimiter=',')

    metadata_options = {
        "metadata_dict": metadata_dict,
        "node_props": node_props,
        "columns": columns,
        "prop2type": prop2type,
    }
    # run the tree_annotate function
    annotated_tree, prop2type = run_tree_annotate(
        tree,
        **metadata_options,
    )

    annotated_newick = annotated_tree.write(props=None, format_root_node=True)
    
    # layouts information
    # level = 0
    # layouts = []
    # column_width = 70
    # padding_x = 1
    # padding_y = 0
    # color_config = {}
    # precomputed_props = columns
    # rectangle_layouts, level, color_dict = tree_plot.get_rectangle_layouts(tree, node_props, 
    #         level, prop2type=prop2type, column_width=column_width, 
    #         padding_x=padding_x, padding_y=padding_y, color_config=color_config)
    # layouts.extend(rectangle_layouts)

    # Store the tree data (in memory or file/db)
    node_props.extend(['name', 'dist'])
    trees[treename] = {
        'tree': tree_data,
        'metadata': metadata,
        'node_props': node_props,
        'annotated_tree': annotated_newick,
        'prop2type': prop2type,
        'layouts': [],
    }

    # Redirect to the /tree/<treename> route after uploading
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
        current_props = tree_info.get('props', list(tree_info['prop2type'].keys()))
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
            
            # Store updated props and layouts back to the tree_info
            tree_info['layouts'] = current_layouts
            tree_info['props'] = selected_props

            # Run the ETE tree explorer in a separate thread with updated props and layout
            def start_explore():
                t.explore(name=treename, layouts=current_layouts, port=5050, open_browser=False, include_props=current_props)

            explorer_thread = threading.Thread(target=start_explore)
            explorer_thread.start()
        else:
            # Run the ETE tree explorer in a separate thread with updated props and layout
            def start_explore():
                t.explore(name=treename, layouts=current_layouts, port=5050, open_browser=False, include_props=current_props)
            explorer_thread = threading.Thread(target=start_explore)
            explorer_thread.start()
            
        return template('explore_tree', treename=treename, tree_info=tree_info, selected_props=current_props)

    return f"Tree '{treename}' not found."



run(host='localhost', port=8080)

