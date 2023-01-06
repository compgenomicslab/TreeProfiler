from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.action_chains import ActionChains
import time
import os
from ete4 import Tree, PhyloTree
from multiprocessing import Process

DRIVER = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'chromedriver')

def tree_session(tree, layouts, port):
    #t = Tree(TREEFILE, format=1)
    tree.explore(tree_name='example', layouts=layouts, show_leaf_name=True, 
                show_branch_length=True, show_branch_support=True, port=port,
                custom_api={}, custom_route={}) 

def snap_tree(port):
    url = 'http://127.0.0.1:{}/'.format(str(port))

    options = Options()
    options.headless = True

    driver = webdriver.Chrome(executable_path=DRIVER, options=options)
    driver.get(url)
    time.sleep(0.5)

    actions = ActionChains(driver)
    actions.send_keys('d')
    actions.perform()
    time.sleep(0.5)
    driver.quit()

def plot(tree, layouts, port, plot_file):
    p1 = Process(target=tree_session, args=(tree, layouts, port,))
    p2 = Process(target=snap_tree, args=(port,))

    p1.start()
    time.sleep(1)
    p2.start()
    time.sleep(2)
    p1.terminate()
    p1.join()
#main(annotated_tree, layouts, port, plot_file)