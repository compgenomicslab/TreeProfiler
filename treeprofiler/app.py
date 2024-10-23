from bottle import route, run, request, redirect, template
import threading

from ete4 import Tree

# In-memory storage for simplicity (you can replace this with a database or file storage)
trees = {}

@route('/')
def upload_tree():
    return template('upload_tree')

@route('/upload', method='POST')
def do_upload():
    treename = request.forms.get('treename')
    tree_data = request.forms.get('tree')
    aln_data = request.forms.get('aln')

    # Store the tree data (in memory or file/db)
    trees[treename] = {
        'tree': tree_data,
        'alignment': aln_data
    }

    # Redirect to the /tree/<treename> route after uploading
    return redirect(f'/tree/{treename}')

@route('/tree/<treename>')
def show_tree(treename):
    # Fetch the tree data by its name
    tree_info = trees.get(treename)
    if tree_info:
        return template('tree_details', treename=treename, tree=tree_info['tree'], aln=tree_info['alignment'])
    else:
        return f"Tree '{treename}' not found."

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

