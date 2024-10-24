from bottle import route, run, request, redirect, template, static_file
import threading

from ete4 import Tree
from treeprofiler.tree_annotate import run_tree_annotate  # or other functions you need

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

    # Store the tree data (in memory or file/db)
    trees[treename] = {
        'tree': tree_data,
        'metadata': metadata
    }

    # here start the tree annotation process
    #tree, eteformat_flag = utils.validate_tree(args.tree, args.input_type, args.internal)
    # load the tree
    # parse the metadata
    # run the tree_annotate function

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


# Route to annotate the tree using treeprofiler
@route('/annotate_tree/<treename>', method=['POST'])
def annotate_tree(treename):
    tree_info = trees.get(treename)
    if not tree_info:
        return "Tree not found"

    # Extract metadata and call run_tree_annotate here with the relevant options
    annotated_tree = run_tree_annotate(tree_info['tree'], metadata=tree_info['metadata'])
    
    # Save the annotated tree (e.g., Newick format, TSV, etc.)
    # You could also store it in memory or generate download links
    
    return "Tree annotation complete!"

@route('/explore_tree/<treename>')
def explore_tree(treename):
    tree_info = trees.get(treename)
    if tree_info:
        # Convert Newick format tree string into an ETE4 Tree object
        t = Tree(tree_info['tree'])

        # Run the ETE tree explorer in a separate thread
        def start_explore():
            t.explore(name=treename, port=5050, open_browser=True)

        explorer_thread = threading.Thread(target=start_explore)
        explorer_thread.start()

        return f"ETE4 explorer launched for tree '{treename}'. It should open in a new browser tab."

    return f"Tree '{treename}' not found."

run(host='localhost', port=8080)

