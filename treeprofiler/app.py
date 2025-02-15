from bottle import route, run, request, redirect, template, static_file, response, ServerAdapter
from bottle import Bottle
import requests
import threading
from tempfile import NamedTemporaryFile
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
import json
import time
from wsgiref.simple_server import WSGIServer, WSGIRequestHandler


from ete4 import Tree
from treeprofiler.tree_annotate import run_tree_annotate, run_array_annotate, parse_emapper_annotations, parse_csv, parse_tsv_to_array, name_nodes  # or other functions you need
from treeprofiler import tree_plot
from treeprofiler.src import utils
from treeprofiler import layouts

from bottle import TEMPLATE_PATH

# Set the template directory
TEMPLATE_PATH.append(os.path.join(os.path.dirname(__file__), 'views'))

# In-memory storage for chunks and complete files
trees = {}
uploaded_chunks = {}
default_paired_color = [
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
continuous_colormaps = [
    # Sequential
    'default', 'magma', 'inferno', 'plasma', 'viridis', 'cividis', 'twilight', 'twilight_shifted', 'turbo',
    'Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'OrRd', 'Oranges', 
    'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd',
    
    # Diverging
    'BrBG', 'PRGn', 'PiYG', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 
    'bwr', 'seismic',
    
    # Cyclic
    'twilight', 'twilight_shifted', 'hsv',
    
    # Miscellaneous
    'cubehelix', 'rainbow', 'nipy_spectral', 'gist_earth', 'terrain', 'ocean'
]
categorical_layout_list = [
    'rectangle-layout',
    'label-layout',
    'colorbranch-layout',
    'categorical-bubble-layout',
    'piechart-layout',
    'background-layout',
    'categorical-matrix-layout',
    'profiling-layout',
    'nodesymbol-layout',
]

categorical_prefix = [
    "rectangle",
    "label",
    "colorbranch",
    "categorical-bubble",
    "piechart",
    "background",
    "categorical-matrix",
    "profiling",
    'textbranch',
    'nodesymbol'
]

numerical_prefix = [
    "heatmap",
    "barplot",
    "branchscore",
    "numerical-bubble",
    "numerical-matrix"
]

binary_prefix = [
    "binary",
    "binary-matrix"
]

# Global control variable for restarting
job_status = {}  # Dictionary to store job statuses
# Global variables
app = Bottle()
server_thread = None
server_instance = None
stop_event = threading.Event()

class CustomWSGIServer(WSGIServer):
    """Custom WSGI server to enable graceful shutdown."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.running = True

    def serve_forever(self):
        """Serve requests until `running` is set to False."""
        while self.running:
            self.handle_request()

    def shutdown(self):
        """Stop the server."""
        self.running = False


class CustomServerAdapter(ServerAdapter):
    """Custom server adapter for Bottle to use CustomWSGIServer."""
    def run(self, handler):
        self.options['handler_class'] = WSGIRequestHandler
        self.server = CustomWSGIServer((self.host, self.port), self.options.get('handler_class', WSGIRequestHandler))
        self.server.set_app(handler)  # Set the Bottle app
        self.server.serve_forever()

    def stop(self):
        if hasattr(self, 'server'):
            self.server.shutdown()


def run_server():
    """Run the Bottle app."""
    global server_instance
    server_instance = CustomServerAdapter(host='localhost', port=8081)
    app.run(server=server_instance)


def start_server():
    """Start the Bottle app in a separate thread."""
    global server_thread, stop_event
    stop_event.clear()
    server_thread = threading.Thread(target=run_server)
    server_thread.start()


@app.route('/stop_explore')
def stop_explore():
    """Stop and restart the server."""
    global server_thread, stop_event, job_status
    job_status = {} # clean up the job status
    print("Stopping server for restart...")

    # Signal the server thread to stop
    stop_event.set()
    return redirect('/')

@app.route('/upload', method='POST')
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
        "is_annotated_tree": request.forms.get('isAnnotatedTree') == 'true',
        
        "metadata": request.forms.get('metadata'),
        "separator": request.forms.get('separator'),

        "matrix": request.forms.get('matrix'),
        "matrix_name": request.forms.get('matrixFiles'),
        "matrix_separator": request.forms.get('matrixSeparator'),

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

        "emapper": request.forms.get('emapper'),
        "pfam": request.forms.get('pfam'),
        "summary_methods": request.forms.get('summary_methods')
    }

    # Start the upload job in a new thread
    threading.Thread(target=process_upload_job, args=(job_args,)).start()
    response.content_type = 'application/json'
    return json.dumps({"job_id": treename})

def process_upload_job(job_args):
    """Function to handle the actual processing of the uploaded data."""
    eteformat_flag = False
    
    treename = job_args['treename']
    # Initialize variables and parse arguments
    separator = "\t" if job_args.get("separator") == "<tab>" else job_args.get("separator", ",")
    tree_data = job_args.get("tree_data")
    treeparser = job_args.get("treeparser")
    is_annotated_tree = job_args.get("is_annotated_tree")
    
    summary_methods = job_args.get("summary_methods")
    
    # Parse summary methods JSON, if provided
    column2method = json.loads(summary_methods) if summary_methods else {}
    column2method['alignment'] = 'none'
    
    # Retrieve paths for uploaded files
    tree_file_path = uploaded_chunks.get(treename, {}).get("tree")
    metadata_file_paths = uploaded_chunks.get(treename, {}).get("metadata")
    metadata_file_list = list(metadata_file_paths.values())
    
    matrix_file_paths = uploaded_chunks.get(treename, {}).get("matrix")

    matrix_file_list = list(matrix_file_paths.values())
    matrix_separator = "\t" if job_args.get("matrix_separator") == "<tab>" else job_args.get("matrix_separator", ",")

    alignment_file_path = uploaded_chunks.get(treename, {}).get("alignment")
    pfam_file_path = uploaded_chunks.get(treename, {}).get("pfam")
    emapper_file_path = uploaded_chunks.get(treename, {}).get("emapper")
    
    # Load the tree data
    if tree_data:
        tree = utils.ete4_parse(tree_data, internal_parser=treeparser)
    elif tree_file_path and os.path.exists(tree_file_path):
        
        # # load tree with newick
        # with open(tree_file_path, 'r') as f:
        #     tree = utils.ete4_parse(f.read(), internal_parser=treeparser)
        
        # Remove the temporary file after usage
        
        tree, eteformat_flag = utils.validate_tree(tree_file_path, 'auto', treeparser)
        os.remove(tree_file_path)

    if is_annotated_tree:
        prop2type = {}
        node_props = []
        for path, node in tree.iter_prepostorder():
            prop2type.update(utils.get_prop2type(node))

        avail_props = list(prop2type.keys())
        node_props.extend(avail_props)
        del avail_props[avail_props.index('name')]
        del avail_props[avail_props.index('dist')]
        if 'support' in avail_props:
            del avail_props[avail_props.index('support')]

        # convert list data to string
        list_keys = [key for key, value in prop2type.items() if value == list]
        
        # Replace all commas in the tree with '||'
        list_sep = '||'
        for node in tree.leaves():
            for key in list_keys:
                if node.props.get(key):
                    
                    list2str = list_sep.join(map(str, node.props.get(key)))
                    node.add_prop(key, list2str)

        annotated_newick = tree.write(props=avail_props, 
        parser=utils.get_internal_parser(treeparser),format_root_node=True)
        
        # check if alignment pro is there
        if 'alignment' in prop2type.keys():
            alignment_annotation = True
        else:
            alignment_annotation = False

        # check if domain pro is there
        if 'dom_arq' in prop2type.keys():
            domain_annotation = True
        else:
            domain_annotation = False

        # check taxonomic annotation
        taxonomic_props = ['rank', 'sci_name', 'taxid', 'named_lineage']
        
        all_present = all(prop in prop2type.keys() for prop in taxonomic_props)
        if all_present:
            taxonomic_annotation = True
        else:
            taxonomic_annotation = False
        
        # Store the processed data
        trees[treename] = {
            'tree': tree_data,
            'treeparser': treeparser,
            #'columns': columns,
            #'metadata': job_args.get("metadata"),
            'node_props': node_props,
            "updated_tree": '',
            'annotated_tree': annotated_newick,
            'prop2type': prop2type,
            'layouts': [],
            'taxonomic_annotation': taxonomic_annotation,
            'rank_list':[],
            'alignment_annotation': alignment_annotation,
            'domain_annotation': domain_annotation,
        }
    else:
        # Process metadata
        metadata_options = {}
        columns = {}

        if metadata_file_list:
            metadata_dict, node_props, columns, prop2type = parse_csv(metadata_file_list, delimiter=separator)
        else:
            metadata_dict = {}
            node_props = []
            prop2type = {}
        
        node_props = list(prop2type.keys())

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
        if emapper_file_path:
            emapper_metadata_dict, emapper_node_props, emapper_columns = parse_emapper_annotations(emapper_file_path)
            metadata_dict = utils.merge_dictionaries(metadata_dict, emapper_metadata_dict)
            node_props.extend(emapper_node_props)
            columns.update(emapper_columns)
            prop2type.update({
                'seed_ortholog': str,
                'evalue': float,
                'score': float,
                'eggNOG_OGs': list,
                'max_annot_lvl': str,
                'COG_category': str,
                'Description': str,
                'Preferred_name': str,
                'GOs': list,
                'EC':str,
                'KEGG_ko': list,
                'KEGG_Pathway': list,
                'KEGG_Module': list,
                'KEGG_Reaction':list,
                'KEGG_rclass':list,
                'BRITE':list,
                'KEGG_TC':list,
                'CAZy':list,
                'BiGG_Reaction':list,
                'PFAMs':list
            })
            
        emapper_options = {
            "emapper_pfam": pfam_file_path if pfam_file_path else None,
            "emapper_mode": True if emapper_file_path else False
        }

        # Run annotation
        threads = 6
        annotated_tree, prop2type = run_tree_annotate(
            tree,
            **metadata_options,
            **emapper_options,
            **taxonomic_options,
            **analytic_options,
            **alignment_options,
            column2method=column2method,
            threads=threads
        )

        # Process Matrix
        node_props_array = []
        if matrix_file_list:
            tmp2filename = {os.path.basename(value): key for key, value in matrix_file_paths.items()}
            array_dict = parse_tsv_to_array(matrix_file_list, delimiter=matrix_separator)
            filename2array = {}

            for key, value in array_dict.items():
                filename2array[tmp2filename[key]] = value
            
            annotated_tree = run_array_annotate(annotated_tree, filename2array, column2method=column2method)
            # update prop2type
            for filename in filename2array.keys():
                prop2type[filename] = list
                prop2type[utils.add_suffix(filename, 'avg')] = list
                prop2type[utils.add_suffix(filename, 'max')] = list
                prop2type[utils.add_suffix(filename, 'min')] = list
                prop2type[utils.add_suffix(filename, 'sum')] = list
                prop2type[utils.add_suffix(filename, 'std')] = list

            node_props_array = list(filename2array.keys())

        if job_args.get("taxon_column"):
            rank_list = sorted(list(set(utils.tree_prop_array(annotated_tree, 'rank'))))
        else:
            rank_list = []

        # Post-processing of annotated tree properties
        list_keys = [key for key, value in prop2type.items() if value == list]

        list_sep = '||'
        for node in annotated_tree.traverse():
            for key in list_keys:
                if node.props.get(key):
                    cont2str = list(map(str, node.props.get(key)))
                    list2str = list_sep.join(cont2str)
                    node.add_prop(key, list2str)
        
        # Name the nodes
        annotated_tree = name_nodes(annotated_tree)

        #avail_props = [key for key in prop2type.keys() if key not in ['name', 'dist', 'support']]
        avail_props = list(prop2type.keys())

        annotated_newick = annotated_tree.write(props=avail_props, format_root_node=True)

        # Add node properties for display
        node_props = metadata_options.get('node_props', [])
        if node_props_array:
            node_props.extend(node_props_array)

        # check if tree has support values
        sample_node = annotated_tree.children[0]
        if 'support' in sample_node.props:
            node_props.extend(['name', 'dist', 'support'])
        else:
            node_props.extend(['name', 'dist'])
        if job_args.get("taxon_column"):
            taxonomic_props = ['rank', 'sci_name', 'taxid', 'evoltype', 'dup_sp', 'dup_percent', 'lca']
            node_props.extend(taxonomic_props)
        
        # Store the processed data
        trees[treename] = {
            'tree': tree_data,
            'treeparser': treeparser,
            #'columns': columns,
            #'metadata': job_args.get("metadata"),
            'node_props': node_props,
            "updated_tree": "",
            'annotated_tree': annotated_newick,
            'prop2type': prop2type,
            'layouts': [],
            'taxonomic_annotation': True if job_args.get("taxon_column") else False,
            'rank_list':rank_list,
            'alignment_annotation': True if alignment_file_path else False,
            'domain_annotation': True if pfam_file_path else False,
        }

    # Cleanup temporary alignment file after usage
    if alignment_file_path and os.path.exists(alignment_file_path):
        os.remove(alignment_file_path)

    # Mark job as complete
    job_status[treename] = "complete"
        

    if job_status[treename] == "failed":
        print(f"Error processing job {treename}: {e}")
    
    # Remove stored chunks for the completed job
    if treename in uploaded_chunks:
        del uploaded_chunks[treename]

# Route to serve static files like CSS and JS
@app.route('/static/<filepath:path>')
def server_static(filepath):
    return static_file(filepath, root='./static')

@app.route('/')
def upload_tree():
    return template('upload_tree')

@app.route('/upload_chunk', method='POST')
def upload_chunk():
    # Determine the file type based on the request
    chunk = request.files.get('treeFile') or request.files.get('metadataFile') or request.files.get('alignmentFile') or request.files.get('pfamFile') or request.files.get('emapperFile') or request.files.get('matrixFile')
    chunk_index = int(request.forms.get("chunkIndex"))
    total_chunks = int(request.forms.get("totalChunks"))
    treename = request.forms.get("treename")
    file_id = request.forms.get("fileId", None)  # Use provided fileId or generate one

    # Identify the file type
    if 'treeFile' in request.files:
        file_type = "tree"
    elif 'metadataFile' in request.files:
        file_type = "metadata"
    elif 'matrixFile' in request.files:
        file_type = "matrix"
    elif 'alignmentFile' in request.files:
        file_type = "alignment"
    elif 'emapperFile' in request.files:
        file_type = "emapper"
    elif 'pfamFile' in request.files:
        file_type = "pfam"
    else:
        return "Unknown file type", 400

    # Initialize storage
    if treename not in uploaded_chunks:
        uploaded_chunks[treename] = {
            "tree": {}, 
            "metadata": {}, 
            "matrix": {},
            "alignment": {}, 
            "pfam": {},
            "emapper": {}
        }

    # Handle metadata specifically for multiple files
    if file_type == "metadata" or file_type == "matrix":
        if file_id not in uploaded_chunks[treename][file_type]:
            uploaded_chunks[treename][file_type][file_id] = []
        uploaded_chunks[treename][file_type][file_id].append((chunk_index, chunk.file.read()))

        # Assemble the file when all chunks are received
        if len(uploaded_chunks[treename][file_type][file_id]) == total_chunks:
            uploaded_chunks[treename][file_type][file_id].sort()
            with NamedTemporaryFile(delete=False) as assembled_file:
                for _, part in uploaded_chunks[treename][file_type][file_id]:
                    assembled_file.write(part)
                uploaded_chunks[treename][file_type][file_id] = assembled_file.name
    else:
        # Handle other file types
        if "chunks" not in uploaded_chunks[treename][file_type]:
            uploaded_chunks[treename][file_type]["chunks"] = []
        uploaded_chunks[treename][file_type]["chunks"].append((chunk_index, chunk.file.read()))

        if len(uploaded_chunks[treename][file_type]["chunks"]) == total_chunks:
            uploaded_chunks[treename][file_type]["chunks"].sort()
            with NamedTemporaryFile(delete=False) as assembled_file:
                for _, part in uploaded_chunks[treename][file_type]["chunks"]:
                    assembled_file.write(part)
                uploaded_chunks[treename][file_type] = assembled_file.name

    return "Chunk received"

@app.route('/job_status/<job_id>')
def job_status_check(job_id):
    status = job_status.get(job_id, "not_found")
    response.content_type = 'application/json'
    return {"status": status}


# @app.route('/show_tree/<treename>')
# def show_tree(treename):
#     # Fetch and display the processed tree, if available
#     tree_info = trees.get(treename)
#     if tree_info:
#         return f"Tree '{treename}' processed successfully! Displaying results."
#     return "Processing still in progress or failed."

@app.route('/tree/<treename>/result')
def show_tree(treename):
    # Fetch the tree data by its name
    tree_info = trees.get(treename)
    if tree_info:
        return template('tree_details', treename=treename, tree_info=tree_info)
    else:
        return f"Tree '{treename}' not found."

@app.route('/tree/<treename>')
def tree_status(treename):
    """Serves the job_running.html template to show job status."""
    return template('job_running', treename=treename, job_id=treename)

@app.route('/check_job_status')
def check_job_status():
    """API endpoint to check the status of a given job."""
    job_id = request.query.get('job_id')
    status = job_status.get(job_id, "not_found")
    return status

@app.route('/explore_tree/<treename>', method=['GET', 'POST', 'PUT'])
def explore_tree(treename):
    tree_info = trees.get(treename)
    
    if not tree_info:
        return f"Tree '{treename}' not found."

    # Retrieve and initialize layouts, properties, and tree
    current_layouts = tree_info.get('layouts', [])
    color_config = {}
    layout_manager = {}

    if current_layouts:
        layout_manager = {layout.name: layout for layout in current_layouts}
    
    prop2type = tree_info['prop2type']
    current_props = sorted(list(tree_info['prop2type'].keys()))

    if tree_info['updated_tree']:
        t = Tree(tree_info['updated_tree'])
    else:
        t = Tree(tree_info['annotated_tree'])
    
    
    # Default configuration settings
    default_configs = {
        "level": 1,
        "column_width": 70,
        "padding_x": 1,
        "padding_y": 0,
        "color_config": color_config,
        "internal_num_rep": 'avg'
    }
    tree_info['default_configs'] = default_configs

    # Retrieve or initialize layouts_metadata in tree_info
    if 'layouts_metadata' not in tree_info:
        tree_info['layouts_metadata'] = []

    layouts_metadata = tree_info['layouts_metadata']  # Reference to the persisted metadata

    # Process POST request
    if request.method == 'POST':
        # Check if the reset action is triggered
        if request.forms.get('reset_tree'):
            print("Resetting tree...")
            # Reset the tree by assigning `annotated_tree` to `updated_tree`
            tree_info['updated_tree'] = None
            tree_info['layouts'] = []  # Clear applied layouts
            current_layouts = []  # Clear layouts
            layouts_metadata.clear()  # Clear metadata
        else:
            layers_data = request.forms.get('layers')
            level = int(request.forms.get('level', default_configs['level']))
            column_width = int(request.forms.get('column_width', default_configs['column_width']))
            padding_x = int(request.forms.get('padding_x', default_configs['padding_x']))
            padding_y = int(request.forms.get('padding_y', default_configs['padding_y']))
            internal_num_rep = request.forms.get('internal_num_rep', default_configs['internal_num_rep'])
            color_config = default_configs.get('color_config', {})

            if layers_data:
                layers = json.loads(layers_data)
                existing_layouts = {layout['layout_name']: layout for layout in layouts_metadata}
                updated_metadata = []
                
                for layer in layers:
                    # Query queryType
                    query_type = layer.get('queryType', '')
                    if query_type == 'rank_limit':
                        rank_selection = layer.get('rankSelection')
                        t = utils.taxatree_prune(t, rank_limit=rank_selection)
                        tree_info['updated_tree'] = t.write(props=current_props, format_root_node=True)
                    elif query_type == 'prune':
                        # prune tree by condition 
                        query_box = layer.get('query', '')
                        query_strings = convert_query_string(query_box)
                        t = utils.conditional_prune(t, query_strings, prop2type)
                        tree_info['updated_tree'] = t.write(props=current_props, format_root_node=True)
                    
                    # Process each layer individually without altering its structure
                    current_layouts, current_props, level, color_config = process_layer(
                        t, layer, tree_info, current_layouts, current_props, level,
                        column_width, padding_x, padding_y, color_config, internal_num_rep, default_paired_color, 
                    )
                    
                    
                    # Update layouts_metadata after process_layer
                    #layouts_metadata.clear()  # Reset to avoid duplicates or outdated data
                    # Update layouts_metadata without clearing it
                    for layout in current_layouts:
                        layout_name = layout.name
                        
                        if layout_name in existing_layouts:
                            # Update existing layout
                            existing_layout = existing_layouts[layout_name]
                            existing_layout['config'].update({
                                "level": getattr(layout, 'column', level),
                                "column_width": getattr(layout, 'column_width', default_configs['column_width']),
                                "padding_x": getattr(layout, 'padding_x', default_configs['padding_x']),
                                "padding_y": getattr(layout, 'padding_y', default_configs['padding_y']),
                                "internal_num_rep": getattr(layout, 'internal_num_rep', default_configs['internal_num_rep']),
                            })
                            existing_layout['layer'].update(layer.get('layer', {}))
                            updated_metadata.append(existing_layout)
                        else:
                            layout_prefix = layout.name.split('_')[0].lower() # get the layout prefix 
                            # taxonomic layout
                            
                            if layout_prefix.startswith('taxa'):
                                layout_config = {
                                    "layout_name": layout.name,
                                    "applied_props": [],
                                    "config": {
                                        "level": getattr(layout, 'column', level),
                                        "column_width": getattr(layout, 'column_width', default_configs['column_width']),
                                        "padding_x": getattr(layout, 'padding_x', default_configs['padding_x']),
                                        # "padding_y": getattr(layout, 'padding_y', default_configs['padding_y']),
                                        # "internal_num_rep": getattr(layout, 'internal_num_rep', default_configs['internal_num_rep']),
                                    },
                                    "layer": {}
                                }
                                
                                updated_metadata.append(layout_config)
                                layout_manager[layout.name] = layout
                            
                            # alignment layout
                            elif layout_prefix in ['alignment', 'domain']:
                                layout_config = {
                                    "layout_name": layout.name,
                                    "applied_props": [layout_prefix],
                                    "config": {
                                        # "level": getattr(layout, 'column', level),
                                        # "column_width": getattr(layout, 'column_width', default_configs['column_width']),
                                        # "padding_x": getattr(layout, 'padding_x', default_configs['padding_x']),
                                        # "padding_y": getattr(layout, 'padding_y', default_configs['padding_y']),
                                        # "internal_num_rep": getattr(layout, 'internal_num_rep', default_configs['internal_num_rep']),
                                        # "color_config": getattr(layout, 'color_config', {})
                                    },
                                    "layer": {}
                                }
                                if layout_prefix == "alignment":
                                    # basic 
                                    layout_config = {
                                        "layout_name": layout.name,  # Retrieve layout name from processed layouts
                                        "applied_props": [layout_prefix],  # Props linked to this layout
                                        "config": {
                                            # "level": getattr(layout, 'column', level),
                                            #"column_width": getattr(layout, 'column_width', default_configs['column_width']),
                                        },
                                        "layer": {}
                                    }
                                    layout_config['layer']['algStart'] = layer.get('algStart', '')
                                    layout_config['layer']['algEnd'] = layer.get('algEnd', '')

                                updated_metadata.append(layout_config)
                                layout_manager[layout.name] = layout
                            
                            # profiling layout
                            elif layout_prefix == 'numerical-matrix':
                                # Append with applied props and full config
                                name = layout.name
                                applied_props = layout.matrix_props
                                layout_config = {
                                    "layout_name": name,  # Retrieve layout name from processed layouts
                                    "applied_props": applied_props,  # Props linked to this layout
                                    "config": {
                                        "level": getattr(layout, 'column', level),
                                        "column_width": getattr(layout, 'column_width', default_configs['column_width']),
                                        "padding_x": getattr(layout, 'padding_x', default_configs['padding_x']),
                                        "padding_y": getattr(layout, 'padding_y', default_configs['padding_y']),
                                        "internal_num_rep": getattr(layout, 'internal_num_rep', default_configs['internal_num_rep']),
                                        #"color_config": getattr(layout, 'color_config', {})
                                    },
                                    "layer": {}
                                }
                                
                                minval, maxval = layout.value_range
                                layout_config['layer']['maxVal'] = maxval
                                layout_config['layer']['minVal'] = minval
                                layout_config['layer']["colorMin"] = layer.get('colorMin', '#0000ff')
                                layout_config['layer']["colorMid"] = layer.get('colorMid', '#ffffff')
                                layout_config['layer']["colorMax"] = layer.get('colorMax', '#ff0000')

                                updated_metadata.append(layout_config)
                                layout_manager[layout.name] = layout

                            elif layout_prefix in ['profiling', 'categorical-matrix', 'binary-matrix']:
                                # Append with applied props and full config
                                name = layout.name
                                applied_props = layout.matrix_props
                                layout_config = {
                                    "layout_name": name,  # Retrieve layout name from processed layouts
                                    "applied_props": applied_props,  # Props linked to this layout
                                    "config": {
                                        "level": getattr(layout, 'column', level),
                                        "column_width": getattr(layout, 'column_width', default_configs['column_width']),
                                        "padding_x": getattr(layout, 'padding_x', default_configs['padding_x']),
                                        "padding_y": getattr(layout, 'padding_y', default_configs['padding_y']),
                                        # "internal_num_rep": getattr(layout, 'internal_num_rep', default_configs['internal_num_rep']),
                                        #"color_config": getattr(layout, 'color_config', {})
                                    },
                                    "layer": {}
                                }
                                updated_metadata.append(layout_config)
                                layout_manager[layout.name] = layout
                            
                            # auto query layout
                            elif layout_prefix in ['collapsed-by', 'highlighted-by']:
                                # Append with applied props and full config
                                name = layout.name
                                conditions = list(layout.color2conditions.values())
                                layout_config = {
                                    "layout_name": name,  # Retrieve layout name from processed layouts
                                    "applied_props": conditions,  # Props linked to this layout
                                    "config": {
                                        # "level": getattr(layout, 'column', level),
                                        # "column_width": getattr(layout, 'column_width', default_configs['column_width']),
                                        # "padding_x": getattr(layout, 'padding_x', default_configs['padding_x']),
                                        # "padding_y": getattr(layout, 'padding_y', default_configs['padding_y']),
                                        # "internal_num_rep": getattr(layout, 'internal_num_rep', default_configs['internal_num_rep']),
                                        #"color_config": getattr(layout, 'color_config', {})
                                    },
                                    "layer": {}
                                }
                                updated_metadata.append(layout_config)
                                layout_manager[layout.name] = layout
                                
                            else:
                                
                                # Append with applied props and full config
                                name = layout.name
                                applied_props = layout.prop
                                
                                # categorical layout
                                if layout_prefix in categorical_prefix:
                                    # basic 
                                    color_config = color_config.get(applied_props, {})
                                    layout_config = {
                                        "layout_name": name,  # Retrieve layout name from processed layouts
                                        "applied_props": [applied_props],  # Props linked to this layout
                                        "config": {
                                            "level": getattr(layout, 'column', level),
                                            "column_width": getattr(layout, 'column_width', default_configs['column_width']),
                                            "padding_x": getattr(layout, 'padding_x', default_configs['padding_x']),
                                            "padding_y": getattr(layout, 'padding_y', default_configs['padding_y']),
                                            #"color_config": color_config
                                        },
                                        "layer": {}
                                    }
                                    if layout_prefix == 'textbranch':
                                        layout_config['layer']['textPosition'] = layer.get('textPosition', 'branch_bottom')
                                        if layer.get('textunicolorColor'):
                                            layout_config['layer']['istextUnicolor'] = 'True'
                                        else:
                                            layout_config['layer']['istextUnicolor'] = 'False'
                                        layout_config['layer']['textColorScheme'] = layer.get('textColorScheme')
                                        layout_config['layer']['textunicolorColor'] = layer.get('textunicolorColor')
                                    
                                    else:
                                        layout_config['layer']['categoricalColorscheme'] = layer.get('categoricalColorscheme', 'default')   
                                
                                    # numerical layout
                                
                                elif layout_prefix in ['circlenode', 'squarenode', 'trianglenode']:
                                    # basic 
                                    color_config = color_config.get(applied_props, {})
                                    layout_config = {
                                        "layout_name": name,  # Retrieve layout name from processed layouts
                                        "applied_props": [applied_props],  # Props linked to this layout
                                        "config": {
                                            "level": getattr(layout, 'column', level),
                                            "column_width": getattr(layout, 'column_width', default_configs['column_width']),
                                            "padding_x": getattr(layout, 'padding_x', default_configs['padding_x']),
                                            "padding_y": getattr(layout, 'padding_y', default_configs['padding_y']),
                                            #"color_config": color_config
                                        },
                                        "layer": {}
                                    }
                                    layout_config['layer']['categoricalColorscheme'] = layer.get('categoricalColorscheme', 'default')
                                    layout_config['layer']['symbolOption'] = layer.get('symbolOption', 'circle')
                                    layout_config['layer']['symbolSize'] = float(layer.get('symbolSize', 5))
                                    layout_config['layer']['fgopacity'] = float(layer.get('fgopacity', 0.8))

                                elif layout_prefix in numerical_prefix:
                                    if layout_prefix == 'barplot':
                                        # basic 
                                        layout_config = {
                                            "layout_name": name,  # Retrieve layout name from processed layouts
                                            "applied_props": [applied_props],  # Props linked to this layout
                                            "config": {
                                                "level": getattr(layout, 'column', level),
                                                #"column_width": getattr(layout, 'column_width', default_configs['column_width']),
                                                "padding_x": getattr(layout, 'padding_x', default_configs['padding_x']),
                                                "padding_y": getattr(layout, 'padding_y', default_configs['padding_y']),
                                                "internal_num_rep": getattr(layout, 'internal_num_rep', default_configs['internal_num_rep']),
                                                #"color_config": color_config.get(applied_props, {})
                                            },
                                            "layer": {}
                                        }
                                        layout_config['layer']['barplotWidth'] = int(layer.get('barplotWidth', 200))
                                        layout_config['layer']['barplotRange'] = layout.size_range[-1]
                                        layout_config['layer']['barplotColorOption'] = layer.get('barplotColorOption', 'same')
                                        
                                        # layout_config['layer']['barplotScaleProperty'] = layer.get('barplotScaleProperty', '')
                                        # layout_config['layer']['barplotColorScheme'] = layer.get('barplotColorScheme', 'default')
                                        # layout_config['layer']['barplotCustomColor'] = layer.get('barplotCustomColor', '#ff0000')

                                        if layout.color is not None:
                                            layout_config['layer']['barplotColor'] = layout.color
                                            layout_config['layer']['barplotFillProperty'] = None

                                        if layout.colors is not None:
                                            layout_config['layer']['barplotFillProperty'] = layer.get('barplotFillProperty', None)
                                            #layout_config['layer']['barplotColorDict'] = layout.colors
                                            layout_config['layer']['barplotColor'] = None
                                        
                                    elif layout_prefix == 'heatmap':
                                        # basic 
                                        layout_config = {
                                            "layout_name": name,  # Retrieve layout name from processed layouts
                                            "applied_props": [applied_props],  # Props linked to this layout
                                            "config": {
                                                "level": getattr(layout, 'column', level),
                                                "column_width": getattr(layout, 'column_width', default_configs['column_width']),
                                                "padding_x": getattr(layout, 'padding_x', default_configs['padding_x']),
                                                "padding_y": getattr(layout, 'padding_y', default_configs['padding_y']),
                                                "internal_num_rep": getattr(layout, 'internal_num_rep', default_configs['internal_num_rep']),
                                                #"color_config": color_config.get(applied_props, {})
                                            },
                                            "layer": {}
                                        }
                                        minval, maxval = layout.value_range
                                        layout_config['layer']['maxVal'] = maxval
                                        layout_config['layer']['minVal'] = minval
                                        layout_config['layer']['colorMin'] = layer.get('colorMin', '#0000ff')
                                        layout_config['layer']['colorMid'] = layer.get('colorMid', '#ffffff')
                                        layout_config['layer']['colorMax'] = layer.get('colorMax', '#ff0000')
                                    
                                    elif layout_prefix == 'branchscore':
                                        # basic 
                                        layout_config = {
                                            "layout_name": name,  # Retrieve layout name from processed layouts
                                            "applied_props": [applied_props],  # Props linked to this layout
                                            "config": {
                                                # "level": getattr(layout, 'column', level),
                                                # "column_width": getattr(layout, 'column_width', default_configs['column_width']),
                                                # "padding_x": getattr(layout, 'padding_x', default_configs['padding_x']),
                                                # "padding_y": getattr(layout, 'padding_y', default_configs['padding_y']),
                                                "internal_num_rep": getattr(layout, 'internal_num_rep', default_configs['internal_num_rep']),
                                                #"color_config": color_config.get(applied_props, {})
                                            },
                                            "layer": {}
                                        }
                                        minval, maxval = layout.value_range
                                        layout_config['layer']['maxVal'] = maxval
                                        layout_config['layer']['minVal'] = minval
                                        layout_config['layer']['colorMin'] = layer.get('colorMin', '#0000ff')
                                        layout_config['layer']['colorMid'] = layer.get('colorMid', '#ffffff')
                                        layout_config['layer']['colorMax'] = layer.get('colorMax', '#ff0000')
                                    elif layout_prefix == 'numerical-bubble':
                                        # basic 
                                        layout_config = {
                                            "layout_name": name,  # Retrieve layout name from processed layouts
                                            "applied_props": [applied_props],  # Props linked to this layout
                                            "config": {
                                                "level": getattr(layout, 'column', level),
                                                # "column_width": getattr(layout, 'column_width', default_configs['column_width']),
                                                # "padding_x": getattr(layout, 'padding_x', default_configs['padding_x']),
                                                # "padding_y": getattr(layout, 'padding_y', default_configs['padding_y']),
                                                "internal_num_rep": getattr(layout, 'internal_num_rep', default_configs['internal_num_rep']),
                                                #"color_config": color_config.get(applied_props, {})
                                            },
                                            "layer": {}
                                        }
                                        minval, maxval = layout.bubble_rage
                                        layout_config['layer']['maxVal'] = maxval
                                        layout_config['layer']['minVal'] = minval
                                        layout_config['layer']['colorMin'] = layer.get('colorMin', '#0000ff')
                                        layout_config['layer']['colorMid'] = layer.get('colorMid', '#ffffff')
                                        layout_config['layer']['colorMax'] = layer.get('colorMax', '#ff0000')
                                    else:
                                        layout_config = {
                                            "layout_name": name,  # Retrieve layout name from processed layouts
                                            "applied_props": [applied_props],  # Props linked to this layout
                                            "config": {
                                                "level": getattr(layout, 'column', level),
                                                "column_width": getattr(layout, 'column_width', default_configs['column_width']),
                                                "padding_x": getattr(layout, 'padding_x', default_configs['padding_x']),
                                                "padding_y": getattr(layout, 'padding_y', default_configs['padding_y']),
                                                "internal_num_rep": getattr(layout, 'internal_num_rep', default_configs['internal_num_rep']),
                                                #"color_config": color_config.get(applied_props, {})
                                            },
                                            "layer": {}
                                        }
                                        layout_config['layer']["maxVal"] = layer.get('maxVal', '')
                                        layout_config['layer']["minVal"] = layer.get('minVal', '')
                                        layout_config['layer']["colorMin"] = layer.get('colorMin', '#0000ff')
                                        layout_config['layer']["colorMid"] = layer.get('colorMid', '#ffffff')
                                        layout_config['layer']["colorMax"] = layer.get('colorMax', '#ff0000')

                                # binary layout
                                elif layout_prefix in binary_prefix:
                                    # basic 
                                    layout_config = {
                                        "layout_name": name,  # Retrieve layout name from processed layouts
                                        "applied_props": [applied_props],  # Props linked to this layout
                                        "config": {
                                            "level": getattr(layout, 'column', level),
                                            "column_width": getattr(layout, 'column_width', default_configs['column_width']),
                                            "padding_x": getattr(layout, 'padding_x', default_configs['padding_x']),
                                            "padding_y": getattr(layout, 'padding_y', default_configs['padding_y']),
                                            #"color_config": color_config.get(applied_props, {})
                                        },
                                        "layer": {}
                                    }
                                    layout_config['layer']['aggregateOption'] = layer.get('aggregateOption', 'gradient')
                                    layout_config['layer']['selectedColor'] = layout.color
                                
                                # alignment
                                else:
                                    # basic 
                                    layout_config = {
                                        "layout_name": name,  # Retrieve layout name from processed layouts
                                        "applied_props": [applied_props],  # Props linked to this layout
                                        "config": {
                                            # "level": getattr(layout, 'column', level),
                                            # "column_width": getattr(layout, 'column_width', default_configs['column_width']),
                                        },
                                        "layer": {}
                                    }
                                updated_metadata.append(layout_config)
                            
                                #layouts_metadata.append(layout_config)
                                layout_manager[layout.name] = layout

                        layouts_metadata.clear()
                        layouts_metadata.extend(updated_metadata)
    
    
    # Process PUT request for updating layouts 
    if request.method == 'PUT':
        if request.forms.get('updated_metadata'):
            try:
                # layouts_metadata.clear()
                updated_metadata = json.loads(request.forms.get('updated_metadata'))
                tree_info['layouts_metadata'] = updated_metadata
                current_layouts = []
                
                for layout_meta in updated_metadata:
                    layout = layout_manager.get(layout_meta['layout_name'])
                    layout_prefix = layout_meta['layout_name'].split('_')[0].lower()
                    if layout:
                        # for categorical
                        if layout_prefix in categorical_prefix:
                            prop = layout_meta['applied_props'][0]
                            # reset color config
                            if layout_prefix  == 'piechart':
                                original_prop = prop.split('_counter')[0] # get the original prop
                                categorical_color_scheme = layout_meta['layer'].get('categoricalColorscheme', 'default')
                                prop_values = sorted(list(set(utils.tree_prop_array(t, original_prop))))
                                paired_color = get_colormap_hex_colors(categorical_color_scheme, len(prop_values))

                                color_config[prop] = {}
                                color_config[prop]['value2color'] = utils.assign_color_to_values(prop_values, paired_color)
                                color_config[prop]['detail2color'] = {}
                                
                                # change directly in layout
                                layout.column = layout_meta['config']['level']
                                layout.width = layout_meta['config']['column_width']
                                layout.padding_x = layout_meta['config']['padding_x']
                                layout.padding_y = layout_meta['config']['padding_y']
                                layout.color_dict = color_config.get(prop).get('value2color')
                            elif layout_prefix == 'textbranch':

                                layout.position = layout_meta['layer'].get('textPosition')

                                if layout_meta['layer'].get('istextUnicolor') == 'True':
                                    layout.text_color = layout_meta['layer'].get('textunicolorColor')
                                    layout.color_dict = {}
                                else:
                                    layout.text_color = None
                                    text_color_scheme = layout_meta['layer'].get('textColorScheme', None)
                                    prop_values = sorted(list(set(utils.tree_prop_array(t, prop))))
                                    paired_color = get_colormap_hex_colors(text_color_scheme, len(prop_values))
                                    color_config[prop] = {}
                                    color_config[prop]['value2color'] = utils.assign_color_to_values(prop_values, paired_color)
                                    color_config[prop]['detail2color'] = {}
                                    color_dict = color_config.get(prop).get('value2color')
                                    layout.color_dict = color_dict
                                
                                # change directly in layout
                                layout.column = layout_meta['config']['level']
                                layout.width = layout_meta['config']['column_width']
                                layout.padding_x = layout_meta['config']['padding_x']
                                layout.padding_y = layout_meta['config']['padding_y']
                                
                            else:
                                categorical_color_scheme = layout_meta['layer'].get('categoricalColorscheme', 'default')
                                prop_values = sorted(list(set(utils.tree_prop_array(t, prop))))
                                paired_color = get_colormap_hex_colors(categorical_color_scheme, len(prop_values))
                                color_config[prop] = {}
                                color_config[prop]['value2color'] = utils.assign_color_to_values(prop_values, paired_color)
                                color_config[prop]['detail2color'] = {}
                                
                                # change directly in layout
                                layout.column = layout_meta['config']['level']
                                layout.width = layout_meta['config']['column_width']
                                layout.padding_x = layout_meta['config']['padding_x']
                                layout.padding_y = layout_meta['config']['padding_y']
                                layout.color_dict = color_config.get(prop).get('value2color')

                        
                        elif layout_prefix in ['circlenode', 'squarenode', 'trianglenode']:
                            prop = layout_meta['applied_props'][0]
                            layout.symbol = layout_meta['layer'].get('symbolOption', 'circle')
                            layout.symbol_size = float(layout_meta['layer'].get('symbolSize', 5))
                            layout.fgopacity = float(layout_meta['layer'].get('fgopacity', 0.8))
                            
                            categorical_color_scheme = layout_meta['layer'].get('categoricalColorscheme', 'default')
                            prop_values = sorted(list(set(utils.tree_prop_array(t, prop))))
                            paired_color = get_colormap_hex_colors(categorical_color_scheme, len(prop_values))
                            color_config[prop] = {}
                            color_config[prop]['value2color'] = utils.assign_color_to_values(prop_values, paired_color)
                            color_config[prop]['detail2color'] = {}
                            
                            # change directly in layout
                            layout.column = layout_meta['config']['level']
                            layout.width = layout_meta['config']['column_width']
                            layout.padding_x = layout_meta['config']['padding_x']
                            layout.padding_y = layout_meta['config']['padding_y']
                            layout.color_dict = color_config.get(prop).get('value2color')
                        # for binary
                        elif layout_prefix in binary_prefix:
                            prop = layout_meta['applied_props'][0]
                            selected_color = layout_meta['layer'].get('selectedColor', '#ff0000')
                            aggregate_option = layout_meta['layer'].get('aggregateOption', 'gradient')
                            
                            if aggregate_option == 'gradient':
                                aggregate = False
                            else:
                                aggregate = True
                            
                            # change directly in layout
                            layout.aggregate = aggregate
                            layout.color = selected_color
                            layout.width = layout_meta['config']['column_width']
                            layout.padding_x = layout_meta['config']['padding_x']
                            layout.padding_y = layout_meta['config']['padding_y']
                        
                        # for numerical
                        elif layout_prefix in numerical_prefix:
                            prop = layout_meta['applied_props'][0]
                            if layout_prefix == 'numerical-matrix':
                                maxval = layout_meta['layer'].get('maxVal', '')                            
                                minval = layout_meta['layer'].get('minVal', '')
                                
                                color_min = layout_meta['layer'].get('colorMin', '#0000ff')
                                color_mid = layout_meta['layer'].get('colorMid', '#ffffff')
                                color_max = layout_meta['layer'].get('colorMax', '#ff0000')
                                selected_props = layout_meta['applied_props']
                                for index, prop in enumerate(selected_props):
                                    color_config[prop] = {}
                                    color_config[prop]['value2color'] = {}
                                    color_config[prop]['detail2color'] = {}
                                    color_config[prop]['detail2color']['color_max'] = (color_max, maxval)
                                    color_config[prop]['detail2color']['color_mid'] = (color_mid, '')
                                    color_config[prop]['detail2color']['color_min'] = (color_min, minval)
                                matrix, minval, maxval, value2color, results_list, list_props, single_props = tree_plot.numerical2matrix(t, 
                                selected_props, count_negative=True, 
                                internal_num_rep=layout_meta['config']['internal_num_rep'], 
                                color_config=color_config, norm_method='min-max', prop2type=tree_info['prop2type'])
                                
                                # TODO add list_props at this moment
                                if list_props:
                                    pass

                                if single_props:
                                    index_map = {value: idx for idx, value in enumerate(selected_props)}
                                    sorted_single_props = sorted(single_props, key=lambda x: index_map[x])
                                    layout.matrix = matrix
                                    layout.matrix_props = sorted_single_props
                                    layout.value_color = value2color
                                    layout.value_range = [minval, maxval]
                                    layout.column = layout_meta['config']['level']
                                    layout.poswidth = layout_meta['config']['column_width']
                            else:
                                maxval = layout_meta['layer'].get('maxVal', '')                            
                                minval = layout_meta['layer'].get('minVal', '')
                                
                                color_min = layout_meta['layer'].get('colorMin', '#0000ff')
                                color_mid = layout_meta['layer'].get('colorMid', '#ffffff')
                                color_max = layout_meta['layer'].get('colorMax', '#ff0000')
                                color_config[prop] = {}
                                color_config[prop]['value2color'] = {}
                                color_config[prop]['detail2color'] = {}
                                color_config[prop]['detail2color']['color_max'] = (color_max, maxval)
                                color_config[prop]['detail2color']['color_mid'] = (color_mid, '')
                                color_config[prop]['detail2color']['color_min'] = (color_min, minval)
                                
                            if layout_prefix == 'barplot':
                                barplot_width = layout_meta['layer'].get('barplotWidth', 200)
                                barplot_width = int(barplot_width)
                                barplot_scale = layout_meta['layer'].get('barplotScaleProperty', '')
                                barplot_range = layout_meta['layer'].get('barplotRange')
                                
                                if barplot_range == '':
                                    barplot_range = layout.size_range[1]
                                else:
                                    barplot_range = float(barplot_range)
                                    
                                size_range = [0, barplot_range]
                                barplot_color_option = layout_meta['layer'].get('barplotColorOption', 'same')
                                
                                # barplot_custom_color = layout_meta['layer'].get('barplotCustomColor', '#ff0000')
                                
                                if barplot_color_option == 'same':
                                    layout.colors = None
                                    layout.color_prop = None
                                    layout.color = layout_meta['layer'].get('barplotColor', None)


                                if barplot_color_option == 'colorby':
                                    layout.color = None
                                    barplot_colorby = layout_meta['layer'].get('barplotFillProperty', '')
                                    if barplot_colorby != '' or barplot_colorby is not None:
                                        prop_values = sorted(list(set(utils.tree_prop_array(t, barplot_colorby))))
                                        color_dict = utils.assign_color_to_values(prop_values, default_paired_color)

                                    layout.colors = color_dict
                                    layout.color_prop = barplot_colorby

                                layout.size_range = size_range
                                layout.width = barplot_width
                                layout.padding_x = layout_meta['config']['padding_x']
                                layout.padding_y = layout_meta['config']['padding_y']
                                layout.internal_rep = layout_meta['config']['internal_num_rep']
                                layout.column = layout_meta['config']['level']

                            if layout_prefix.startswith('branchscore'):
                                layouts = tree_plot.get_branchscore_layouts(
                                    t, [prop], prop2type=tree_info['prop2type'],  
                                    internal_rep=layout_meta['config']['internal_num_rep'],
                                    color_config=color_config
                                )
                                layout = layouts[0] 
                                layout.internal_num_rep = layout_meta['config']['internal_num_rep']
                            
                            if layout_prefix == 'numerical-bubble':
                                # Convert maxval and minval to floats if they exist
                                if maxval is not None and maxval != '':
                                    maxval = float(maxval)
                                else:
                                    maxval = None  # Reset to None if invalid
                                if minval is not None and minval != '':
                                    minval = float(minval)
                                else:
                                    minval = None  # Reset to None if invalid

                                # Ensure bubble_range includes 0 if valid
                                if minval is not None and maxval is not None:
                                    bubble_range = [minval, maxval]
                                else:
                                    bubble_range = []

                                bubble_layouts, level, _ = tree_plot.get_numerical_bubble_layouts(
                                    t, [prop], layout_meta['config']['level'], 
                                    tree_info['prop2type'], internal_rep=layout_meta['config']['internal_num_rep'],
                                    bubble_range=bubble_range, color_config=color_config
                                )
                                layout = bubble_layouts[0]
                                layout.width = layout_meta['config']['column_width']
                                layout.padding_x = layout_meta['config']['padding_x']
                                layout.padding_y = layout_meta['config']['padding_y']
                                layout.internal_num_rep = layout_meta['config']['internal_num_rep']

                            if layout_meta['layout_name'].lower().startswith('heatmap') and layout_meta['layout_name'].lower().endswith('min-max'):
                                layouts, _ = tree_plot.get_heatmap_layouts(
                                    t, [prop], layout_meta['config']['level'], 
                                    column_width=layout_meta['config']['column_width'],
                                    padding_x=layout_meta['config']['padding_x'], padding_y=layout_meta['config']['padding_y'], 
                                    internal_rep=layout_meta['config']['internal_num_rep'],
                                    color_config=color_config, norm_method='min-max'
                                )
                                layout = layouts[0]
                                layout.width = layout_meta['config']['column_width']
                                layout.padding_x = layout_meta['config']['padding_x']
                                layout.padding_y = layout_meta['config']['padding_y']
                                layout.internal_num_rep = layout_meta['config']['internal_num_rep']
                                
                            if layout_meta['layout_name'].lower().startswith('heatmap') and layout_meta['layout_name'].lower().endswith('mean'):
                                layouts, _ = tree_plot.get_heatmap_layouts(
                                    t, [prop], layout_meta['config']['level'],
                                    column_width=layout_meta['config']['column_width'],
                                    padding_x=layout_meta['config']['padding_x'], padding_y=layout_meta['config']['padding_y'], 
                                    internal_rep=layout_meta['config']['internal_num_rep'],
                                    color_config=color_config, norm_method='mean'
                                )
                                layout = layouts[0]
                                layout.width = layout_meta['config']['column_width']
                                layout.padding_x = layout_meta['config']['padding_x']
                                layout.padding_y = layout_meta['config']['padding_y']
                                layout.internal_num_rep = layout_meta['config']['internal_num_rep']
                                
                            if layout_meta['layout_name'].lower().startswith('heatmap') and layout_meta['layout_name'].lower().endswith('zscore'):
                                layouts, _ = tree_plot.get_heatmap_layouts(
                                    t, [prop], layout_meta['config']['level'],
                                    column_width=layout_meta['config']['column_width'],
                                    padding_x=layout_meta['config']['padding_x'], padding_y=layout_meta['config']['padding_y'], 
                                    internal_rep=layout_meta['config']['internal_num_rep'],
                                    color_config=color_config, norm_method='zscore'
                                )
                                layout = layouts[0]
                                layout.width = layout_meta['config']['column_width']
                                layout.padding_x = layout_meta['config']['padding_x']
                                layout.padding_y = layout_meta['config']['padding_y']
                                layout.internal_num_rep = layout_meta['config']['internal_num_rep']
                            
                            
                        current_layouts.append(layout)

                tree_info['layouts'] = current_layouts

                start_explore_thread(t, treename, current_layouts, current_props)
                #return "Layouts metadata updated successfully."
            except json.JSONDecodeError:
                return "Invalid metadata format.", 400
    # Start the ete exploration thread

    if request.method == 'GET':
        start_explore_thread(t, treename, current_layouts, current_props)

    # Before rendering the template, convert to JSON
    #layouts_json = json.dumps(tree_info['layouts'])
    layouts_metadata_json = json.dumps(layouts_metadata)

    # Render template
    return template(
        'explore_tree_v3',
        treename=treename,
        tree_info=tree_info,
        selected_props=current_props,
        color_schemes=continuous_colormaps,
        layouts_metadata=layouts_metadata,
        layouts_metadata_json = json.dumps(layouts_metadata)
    )

def get_colormap_hex_colors(colormap_name, num_colors):
    if colormap_name == 'default':
        return default_paired_color
    else:
        cmap = plt.get_cmap(colormap_name)
        colors = [mcolors.to_hex(cmap(i / (num_colors - 1))) for i in range(num_colors)]
        return colors

def process_layer(t, layer, tree_info, current_layouts, current_props, level, column_width, padding_x, padding_y, color_config, internal_num_rep, paired_color):
    """
    Processes a single layer, applying user selections and updating layouts, properties, and level.
    """
    layout2colorconfg = {}
    selected_props = layer.get('props', [])
    selected_layout = layer.get('layout', '')
    query_type = layer.get('queryType', '')
    rank_selection = layer.get('rankSelection', '')
    query_box = layer.get('query', '')

    prop2type = tree_info['prop2type']

    # Process color configuration for the current layer
    # categorical settings
    categorical_color_scheme = layer.get('categoricalColorscheme', 'default')

    # numerical settings
    maxval = layer.get('maxVal', '') # this should be automatically calculated
    minval = layer.get('minVal', '') # this should be automatically calculated
    color_min = layer.get('colorMin', '#0000ff')
    color_mid = layer.get('colorMid', '#ffffff')
    color_max = layer.get('colorMax', '#ff0000')

    # barplot settings
    barplot_width = layer.get('barplotWidth', 200)
    barplot_scale = layer.get('barplotScaleProperty', '')
    barplot_range = layer.get('barplotRange', '')
    barplot_color_option = layer.get('barplotColorOption', 'same')
    barplot_color_scheme = layer.get('barplotColorScheme', 'default')
    barplot_colorby = layer.get('barplotFillProperty', None)
    
    # binary settings
    if selected_layout == 'binary-layout':
        same_color = layer.get('isbinaryUnicolor', True)

        if not same_color:
            bianry_color_scheme = layer.get('binaryColorScheme', 'default')
        else:
            unicolorColor = layer.get('binaryunicolorColor', '#ff0000')

        aggregate_option = layer.get('aggregateOption', 'gradient')
        
        if aggregate_option == 'gradient':
            aggregate = False
        else:
            aggregate = True

    # alignment settings
    alg_start = layer.get('algStart', '')
    alg_end = layer.get('algEnd', '')
    
    # Apply selected layout based on type directly within this function
    if selected_layout in categorical_layout_list:
        for index, prop in enumerate(selected_props):
            prop_values = sorted(list(set(utils.tree_prop_array(t, prop))))
            paired_color = get_colormap_hex_colors(categorical_color_scheme, len(prop_values))
            color_config[prop] = {}
            color_config[prop]['value2color'] = utils.assign_color_to_values(prop_values, paired_color)
            color_config[prop]['detail2color'] = {}
        
        if selected_layout == 'nodesymbol-layout':
            for prop in selected_props:
                color_dict = color_config.get(prop)['value2color']
                symbol = layer.get('symbolOption', 'circle')
                symbol_size = layer.get('symbolSize', 5)
                if symbol_size:
                    symbol_size = float(symbol_size)
                fgopacity = layer.get('fgopacity', 0.8)
                
                layout = layouts.text_layouts.LayoutSymbolNode(f'{symbol}Node_{prop}', prop=prop,
                    column=level, symbol=symbol, symbol_color=None, color_dict=color_dict,
                    symbol_size=symbol_size, 
                    padding_x=padding_x, padding_y=padding_y, fgopacity=fgopacity,
                    scale=True, legend=True, active=True
                )
                level +=1 
                current_layouts.append(layout)
            
        else:
            level, current_layouts = apply_categorical_layouts(t, selected_layout, selected_props, tree_info, current_layouts, level, column_width, padding_x, padding_y, color_config)
        
    elif selected_layout == 'textbranch-layout':
        # text branch

        text_position = layer.get('textPosition', 'branch_bottom')
        text_color_scheme = layer.get('textColorScheme', None)

        text_color = layer.get('textunicolorColor', None)

        if text_color:
            for prop in selected_props:
                layout = layouts.text_layouts.LayoutTextbranch(name='TextBranch_'+prop, 
                column=level, text_color=text_color, color_dict={}, prop=prop, 
                position=text_position, width=column_width, 
                padding_x=padding_x, padding_y=padding_y)
                level +=1 
                current_layouts.append(layout)
        else:
            for index, prop in enumerate(selected_props):
                prop_values = sorted(list(set(utils.tree_prop_array(t, prop))))
                paired_color = get_colormap_hex_colors(text_color_scheme, len(prop_values))
                color_config[prop] = {}
                color_config[prop]['value2color'] = utils.assign_color_to_values(prop_values, paired_color)
                color_config[prop]['detail2color'] = {}
            for prop in selected_props:
                color_dict = color_config.get(prop).get('value2color')
                layout = layouts.text_layouts.LayoutTextbranch(name='TextBranch_'+prop, 
                column=level, text_color=None, color_dict=color_dict, prop=prop, 
                position=text_position, width=column_width, 
                padding_x=padding_x, padding_y=padding_y)
                level +=1 
                current_layouts.append(layout)
    
    elif selected_layout == 'binary-layout':
        for index, prop in enumerate(selected_props):
            prop_values = utils.tree_prop_array(t, prop, leaf_only=True)
            max_count = utils.find_bool_representations(prop_values)
            prop_values = sorted(set(prop_values))
            color_config[prop] = {}
            color_config[prop]['value2color'] = {}
            if same_color:
                color_config[prop]['value2color']['True'] = unicolorColor
            else:
                color_config[prop]['value2color']['True'] = paired_color[index]

        binary_layouts, level, _ = tree_plot.get_binary_layouts(
        t, selected_props, level, prop2type=tree_info['prop2type'],
        column_width=column_width, reverse=False, padding_x=padding_x, padding_y=padding_y,
        color_config=color_config, same_color=same_color, aggregate=aggregate)
        current_layouts.extend(binary_layouts)

    elif selected_layout == 'branchscore-layout':
        for index, prop in enumerate(selected_props):
            color_config[prop] = {}
            color_config[prop]['value2color'] = {}
            color_config[prop]['detail2color'] = {}
            color_config[prop]['detail2color']['color_max'] = (color_max, maxval)
            color_config[prop]['detail2color']['color_mid'] = (color_mid, '')
            color_config[prop]['detail2color']['color_min'] = (color_min, minval)
        branchscore_layouts = tree_plot.get_branchscore_layouts(
            t, selected_props, prop2type=tree_info['prop2type'], 
            padding_x=padding_x, padding_y=padding_y, internal_rep=internal_num_rep,
            color_config=color_config
        )
        current_layouts.extend(branchscore_layouts)
    elif selected_layout == 'barplot-layout':
        for index, prop in enumerate(selected_props):
            if barplot_color_option == 'colorby':
                color_config = None
            elif barplot_color_option == 'custom':
                barplot_custom_color = layer.get('barplotCustomColor', '#ff0000')
                color_config[prop] = {}
                color_config[prop]['value2color'] = {}
                color_config[prop]['detail2color'] = {}
                color_config[prop]['detail2color']['barplot_color'] = (barplot_custom_color, '')
            else:
                paired_color = get_colormap_hex_colors(barplot_color_scheme, len(selected_props))
                color_config[prop] = {}
                color_config[prop]['value2color'] = {}
                color_config[prop]['detail2color'] = {}
                color_config[prop]['detail2color']['barplot_color'] = (paired_color[index], '')
        if barplot_range != '':
            barplot_range = float(barplot_range)
        if barplot_width != '':
            barplot_width = int(barplot_width)
        
        barplot_layouts, level, _ = tree_plot.get_barplot_layouts(
            t, selected_props, level, tree_info['prop2type'], 
            column_width=barplot_width, padding_x=padding_x, padding_y=padding_y,
            internal_rep=internal_num_rep, anchor_column=barplot_scale, color_config=color_config,
            barplot_colorby=barplot_colorby, max_range=barplot_range)
        current_layouts.extend(barplot_layouts)

    elif selected_layout == 'numerical-bubble-layout':
        for index, prop in enumerate(selected_props):
            color_config[prop] = {}
            color_config[prop]['value2color'] = {}
            color_config[prop]['detail2color'] = {}
            color_config[prop]['detail2color']['color_max'] = (color_max, maxval)
            color_config[prop]['detail2color']['color_mid'] = (color_mid, '')
            color_config[prop]['detail2color']['color_min'] = (color_min, minval)
        # Convert maxval and minval to floats if they exist
        if maxval is not None and maxval != '':
            maxval = float(maxval)
        else:
            maxval = None
        if minval is not None and minval != '':
            minval = float(minval)
        else:
            minval = None

        if minval is not None and maxval is not None:
            bubble_range=[minval, maxval]
        else:
            bubble_range = []
        bubble_layouts, level, _ = tree_plot.get_numerical_bubble_layouts(
            t, selected_props, level, tree_info['prop2type'], internal_rep=internal_num_rep,
            bubble_range=bubble_range, color_config=color_config
        )
        current_layouts.extend(bubble_layouts)
    elif selected_layout == 'heatmap-layout':
        for index, prop in enumerate(selected_props):
            color_config[prop] = {}
            color_config[prop]['value2color'] = {}
            color_config[prop]['detail2color'] = {}
            color_config[prop]['detail2color']['color_max'] = (color_max, maxval)
            color_config[prop]['detail2color']['color_mid'] = (color_mid, '')
            color_config[prop]['detail2color']['color_min'] = (color_min, minval)
        heatmap_layouts, level = tree_plot.get_heatmap_layouts(
            t, selected_props, level, column_width=column_width, padding_x=padding_x, padding_y=padding_y, 
            internal_rep=internal_num_rep, color_config=color_config, norm_method='min-max',
            global_scaling=True
        )
        current_layouts.extend(heatmap_layouts)
    elif selected_layout == 'heatmap-mean-layout':
        for index, prop in enumerate(selected_props):
            color_config[prop] = {}
            color_config[prop]['value2color'] = {}
            color_config[prop]['detail2color'] = {}
            color_config[prop]['detail2color']['color_max'] = (color_max, maxval)
            color_config[prop]['detail2color']['color_mid'] = (color_mid, '')
            color_config[prop]['detail2color']['color_min'] = (color_min, minval)
        heatmap_layouts, level = tree_plot.get_heatmap_layouts(t, selected_props, level,
            column_width=column_width, padding_x=padding_x, padding_y=padding_y, 
            internal_rep=internal_num_rep, color_config=color_config, norm_method='mean',
            global_scaling=True)
        current_layouts.extend(heatmap_layouts)
    elif selected_layout == 'heatmap-zscore-layout':
        for index, prop in enumerate(selected_props):
            color_config[prop] = {}
            color_config[prop]['value2color'] = {}
            color_config[prop]['detail2color'] = {}
            color_config[prop]['detail2color']['color_max'] = (color_max, maxval)
            color_config[prop]['detail2color']['color_mid'] = (color_mid, '')
            color_config[prop]['detail2color']['color_min'] = (color_min, minval)
        heatmap_layouts, level = tree_plot.get_heatmap_layouts(t, selected_props, level,
            column_width=column_width, padding_x=padding_x, padding_y=padding_y, 
            internal_rep=internal_num_rep, color_config=color_config, norm_method='zscore',
            global_scaling=True)
        current_layouts.extend(heatmap_layouts)
    elif selected_layout == 'numerical-matrix-layout':
        for index, prop in enumerate(selected_props):
            color_config[prop] = {}
            color_config[prop]['value2color'] = {}
            color_config[prop]['detail2color'] = {}
            color_config[prop]['detail2color']['color_max'] = (color_max, maxval)
            color_config[prop]['detail2color']['color_mid'] = (color_mid, '')
            color_config[prop]['detail2color']['color_min'] = (color_min, minval)
        matrix, minval, maxval, value2color, results_list, list_props, single_props = tree_plot.numerical2matrix(t, 
        selected_props, count_negative=True, internal_num_rep=internal_num_rep, 
        color_config=color_config, norm_method='min-max', prop2type=prop2type)
        if list_props:
            index_map = {value: idx for idx, value in enumerate(selected_props)}
            sorted_list_props = sorted(list_props, key=lambda x: index_map[x])
            for list_prop in sorted_list_props:
                matrix, minval, maxval, value2color = results_list[list_prop]
                matrix_layout = tree_plot.profile_layouts.LayoutPropsMatrixOld(name=f"Numerical-matrix_{list_prop}", 
                    matrix=matrix, matrix_type='numerical', matrix_props=[list_prop], is_list=True, 
                    value_color=value2color, value_range=[minval, maxval], column=level,
                    poswidth=column_width)
                level += 1
                current_layouts.append(matrix_layout)
        if single_props:
            index_map = {value: idx for idx, value in enumerate(selected_props)}
            sorted_single_props = sorted(single_props, key=lambda x: index_map[x])
            matrix_layout = tree_plot.profile_layouts.LayoutPropsMatrixOld(name=f"Numerical-matrix_{sorted_single_props}", 
                matrix=matrix, matrix_type='numerical', matrix_props=sorted_single_props, is_list=False, 
                value_color=value2color, value_range=[minval, maxval], column=level,
                poswidth=column_width)
            level += 1
            current_layouts.append(matrix_layout)
    elif selected_layout in ['acr-discrete-layout', 'acr-continuous-layout', 'ls-layout']:
        for index, prop in enumerate(selected_props):
            prop_values = sorted(list(set(utils.tree_prop_array(t, prop))))
            paired_color = get_colormap_hex_colors(categorical_color_scheme, len(prop_values))
            color_config[prop] = {}
            color_config[prop]['value2color'] = utils.assign_color_to_values(prop_values, paired_color)
            color_config[prop]['detail2color'] = {}
        level, current_layouts = apply_analytic_layouts(t, selected_layout, selected_props, tree_info, current_layouts, level, column_width, padding_x, padding_y, color_config)
    elif selected_layout in ['taxoncollapse-layout', 'taxonclade-layout', 'taxonrectangle-layout']:
        level, current_layouts = apply_taxonomic_layouts(t, selected_layout, selected_props, tree_info, current_layouts, level, column_width, padding_x, padding_y, paired_color)
    elif selected_layout == 'alignment-layout':
        if alg_start != '':
            alg_start = int(alg_start)
        if alg_end != '':
            alg_end = int(alg_end)

        if alg_start != '' and alg_end != '':
            window = [alg_start, alg_end]
        else:
            window = []

        lengh = len(max(utils.tree_prop_array(t, 'alignment'), key=len))
        aln_layout = layouts.seq_layouts.LayoutAlignment(
            name='Alignment_layout', alignment_prop='alignment', 
            column=level, scale_range=lengh,
            window=window, summarize_inner_nodes=True)
        current_layouts.append(aln_layout)
    elif selected_layout == 'domain-layout':
        domain_layout = layouts.seq_layouts.LayoutDomain(name="Domain_layout", prop='dom_arq')
        current_layouts.append(domain_layout)
        
    # Process query actions for the current layer
    current_layouts = apply_queries(query_type, query_box, current_layouts, paired_color, tree_info, level)
    
    current_props.extend(selected_props)

    return current_layouts, current_props, level, color_config

def apply_categorical_layouts(t, selected_layout, selected_props, tree_info, current_layouts, level, column_width, padding_x, padding_y, color_config):
    """
    Applies categorical layouts such as rectangle, label, colorbranch, bubble, piechart, background, and profiling.
    """
    prop2type = tree_info['prop2type']
    
    # Apply specific categorical layout configurations
    if selected_layout == 'rectangle-layout':
        rectangle_layouts, level, _ = tree_plot.get_rectangle_layouts(
            t, selected_props, level, prop2type=prop2type, column_width=column_width,
            padding_x=padding_x, padding_y=padding_y, color_config=color_config
        )
        current_layouts.extend(rectangle_layouts)

    elif selected_layout == 'label-layout':
        label_layouts, level, _ = tree_plot.get_label_layouts(
            t, selected_props, level, prop2type=prop2type, column_width=column_width,
            padding_x=padding_x, padding_y=padding_y, color_config=color_config
        )
        current_layouts.extend(label_layouts)

    elif selected_layout == 'colorbranch-layout':
        colorbranch_layouts, level, _ = tree_plot.get_colorbranch_layouts(
            t, selected_props, level, prop2type=prop2type, column_width=column_width,
            padding_x=padding_x, padding_y=padding_y, color_config=color_config
        )
        current_layouts.extend(colorbranch_layouts)

    elif selected_layout == 'categorical-bubble-layout':
        bubble_layouts, level, _ = tree_plot.get_categorical_bubble_layouts(
            t, selected_props, level, prop2type=prop2type, column_width=column_width,
            padding_x=padding_x, padding_y=padding_y, color_config=color_config
        )
        current_layouts.extend(bubble_layouts)

    elif selected_layout == 'piechart-layout':
        piechart_layouts = tree_plot.get_piechart_layouts(
            t, selected_props, prop2type=prop2type, padding_x=padding_x,
            padding_y=padding_y, color_config=color_config
        )
        current_layouts.extend(piechart_layouts)

    elif selected_layout == 'background-layout':
        background_layouts, level, _ = tree_plot.get_background_layouts(
            t, selected_props, level, prop2type=prop2type
        )
        current_layouts.extend(background_layouts)

    elif selected_layout == 'profiling-layout':
        for profiling_prop in selected_props:
            
            matrix, value2color, all_profiling_values = tree_plot.multiple2matrix(
                t, profiling_prop, prop2type=prop2type, color_config=color_config,
            )
            
            matrix_layout = tree_plot.profile_layouts.LayoutPropsMatrixBinary(
                name=f"Profiling_{profiling_prop}", matrix=matrix, matrix_props=all_profiling_values,
                value_range=[0, 1], value_color=value2color, column=level, poswidth=column_width
            )
            current_layouts.append(matrix_layout)
            level += 1

    elif selected_layout == 'categorical-matrix-layout':
        # Draws matrix array for categorical properties
        matrix, value2color = tree_plot.categorical2matrix(
            t, selected_props, color_config=color_config
        )
        matrix_layout = tree_plot.profile_layouts.LayoutPropsMatrixOld(
            name=f"Categorical-matrix_{selected_props}", matrix=matrix, matrix_type='categorical',
            matrix_props=selected_props, value_color=value2color, column=level, poswidth=column_width
        )
        current_layouts.append(matrix_layout)
        level += 1

    return level, current_layouts

def apply_analytic_layouts(t, selected_layout, selected_props, tree_info, current_layouts, level, column_width, padding_x, padding_y, color_config):
    """
    Applies analytic layouts like acr-discrete, acr-continuous, and ls-layout.
    """
    if selected_layout == 'acr-discrete-layout':
        acr_layouts, level, _ = tree_plot.get_acr_discrete_layouts(
            t, selected_props, level, prop2type=tree_info['prop2type'],
            column_width=column_width, padding_x=padding_x, padding_y=padding_y,
            color_config=color_config
        )
        current_layouts.extend(acr_layouts)
    elif selected_layout == 'acr-continuous-layout':
        acr_layouts = tree_plot.get_acr_continuous_layouts(
            t, selected_props, level, prop2type=tree_info['prop2type'],
            padding_x=padding_x, padding_y=padding_y
        )
        current_layouts.extend(acr_layouts)
    elif selected_layout == 'ls-layout':
        ls_layouts, ls_props = tree_plot.get_ls_layouts(
            t, selected_props, level, prop2type=tree_info['prop2type'],
            padding_x=padding_x, padding_y=padding_y, color_config=color_config
        )
        current_layouts.extend(ls_layouts)
    return level, current_layouts

def apply_taxonomic_layouts(t, selected_layout, selected_props, tree_info, current_layouts, level, column_width, padding_x, padding_y, paired_color):
    """
    Applies taxonomic layouts such as taxoncollapse, taxonclade, and taxonrectangle.
    """
    taxon_color_dict = {}
    taxa_layouts = []
    rank2values = defaultdict(list)

    # Traverse tree nodes to build rank2values dictionary
    for n in t.traverse():
        if n.props.get('lca'):
            lca_dict = utils.string_to_dict(n.props.get('lca'))
            
            for rank, sci_name in lca_dict.items():
                rank2values[rank].append(sci_name)

        # Add rank information if available
        current_rank = n.props.get('rank')
        if current_rank and current_rank != 'Unknown':
            rank2values[current_rank].append(n.props.get('sci_name', ''))

    # Assign colors based on unique values for each rank
    for rank, values in sorted(rank2values.items()):
        unique_values = list(set(values))
        color_dict = utils.assign_color_to_values(unique_values, paired_color)

        # Determine layout type and create layout instances accordingly
        if selected_layout == 'taxoncollapse-layout':
            taxa_layout = layouts.taxon_layouts.TaxaCollapse(
                name="TaxaCollapse_" + rank,
                rank=rank,
                rect_width=column_width,
                color_dict=color_dict,
                column=level
            )
            taxa_layouts.append(taxa_layout)
        elif selected_layout == 'taxonclade-layout':
            taxa_layout = layouts.taxon_layouts.TaxaClade(
                name='TaxaClade_' + rank,
                level=level,
                rank=rank,
                color_dict=color_dict
            )
            taxa_layouts.append(taxa_layout)
        elif selected_layout == 'taxonrectangle-layout':
            taxa_layout = layouts.taxon_layouts.TaxaRectangular(
                name="TaxaRect_" + rank,
                rank=rank,
                rect_width=column_width,
                color_dict=color_dict,
                column=level
            )
            taxa_layouts.append(taxa_layout)

        # Update taxon_color_dict for LayoutSciName and LayoutEvolEvents
        taxon_color_dict[rank] = color_dict

    # Add scientific name and evolutionary events layouts
    taxa_layouts.append(
        layouts.taxon_layouts.LayoutSciName(
            name='Taxa_Scientific name',
            color_dict=taxon_color_dict
        )
    )
    taxa_layouts.append(
        layouts.taxon_layouts.LayoutEvolEvents(
            name='Taxa_Evolutionary events',
            prop="evoltype",
            speciation_color="blue",
            duplication_color="red",
            node_size=3,
            legend=True
        )
    )

    # Extend current layouts with generated taxonomic layouts
    current_layouts.extend(taxa_layouts)
    level += 1

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
            name=f'Highlighted-by_{condition}', color2conditions=color2conditions, 
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
            name=f'Collapsed-by_{condition}', color2conditions=color2conditions, 
            column=level, prop2type=tree_info['prop2type']
        )
        current_layouts.append(c_layout)
    
    return current_layouts

def start_explore_thread(t, treename, current_layouts, current_props):
    """
    Starts the ete exploration in a separate thread.
    """

    def explore():
        t.explore(name=treename, layouts=current_layouts, port=5051, open_browser=False, include_props=current_props)
    
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

if __name__ == "__main__":
    start_server()