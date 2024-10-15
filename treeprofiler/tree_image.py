from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.action_chains import ActionChains
from webdriver_manager.chrome import ChromeDriverManager
from pyvirtualdisplay import Display  # Import pyvirtualdisplay to create a virtual X server
import time
import os
from ete4 import Tree
from multiprocessing import Process
import subprocess
import re
import shutil

# Custom directory to store downloaded drivers
DRIVER_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'drivers')

# Step 1: Get the installed Chrome browser version
def get_chrome_version():
    try:
        version_output = subprocess.check_output(["google-chrome", "--version"], text=True)
        version = re.search(r"(\d+\.\d+\.\d+\.\d+)", version_output).group(1)
        return version
    except Exception as e:
        print(f"Error checking Chrome version: {e}")
        return None

# Step 2: Get the ChromeDriver version stored in the custom path
def get_chromedriver_version(driver_path):
    try:
        version_output = subprocess.check_output([driver_path, "--version"], text=True)
        version = re.search(r"(\d+\.\d+\.\d+\.\d+)", version_output).group(1)
        return version
    except Exception as e:
        print(f"Error checking ChromeDriver version: {e}")
        return None

# Step 3: Check if there's already a driver in the custom path that matches the installed Chrome version
def find_existing_driver(custom_path, chrome_version):
    if not os.path.exists(custom_path):
        return None
    
    for root, dirs, files in os.walk(custom_path):
        for file in files:
            if file == "chromedriver":  # Adjust for your OS (e.g., chromedriver.exe on Windows)
                driver_path = os.path.join(root, file)
                driver_version = get_chromedriver_version(driver_path)
                if driver_version and driver_version.split('.')[0] == chrome_version.split('.')[0]:
                    print(f"Found matching ChromeDriver (v{driver_version}) in {driver_path}")
                    return driver_path
    return None

# Step 4: Download the ChromeDriver and move it to the custom path
def download_and_move_driver(custom_path):
    # Download the driver using ChromeDriverManager
    downloaded_driver = ChromeDriverManager().install()

    if not os.path.exists(custom_path):
        os.makedirs(custom_path)

    driver_filename = "chromedriver"  # Adjust for your OS (e.g., chromedriver.exe on Windows)
    destination = os.path.join(custom_path, driver_filename)

    # Copy the downloaded driver to the custom path
    shutil.copy(downloaded_driver, destination)
    return destination

# Step 5: Set up the ChromeDriver in headless mode
# def setup_driver(custom_path):
#     # Check if Chrome is available in /usr/bin
#     chrome_binary = shutil.which("google-chrome")
    
#     if chrome_binary:
#         print(f"Google Chrome found at {chrome_binary}")
#         return chrome_binary, True  # Return the Chrome binary and flag it as found
#     else:
#         print("Google Chrome not found in /usr/bin. Using ChromeDriver instead.")
#         chrome_version = get_chrome_version()
#         if chrome_version:
#             print(f"Installed Chrome version: {chrome_version}")
#         else:
#             print("Could not determine Chrome version. Downloading ChromeDriver automatically.")
        
#         # Check for existing ChromeDriver in the custom path
#         existing_driver = find_existing_driver(custom_path, chrome_version)
        
#         if existing_driver:
#             print("Using existing ChromeDriver.")
#             driver_path = existing_driver
#         else:
#             print("No matching ChromeDriver found. Downloading new driver.")
#             driver_path = download_and_move_driver(custom_path)

#         return driver_path, False  # Return the ChromeDriver path and flag it as ChromeDriver

def setup_driver(custom_path):
    chrome_version = get_chrome_version()
    if chrome_version:
        print(f"Installed Chrome version: {chrome_version}")
    else:
        print("Could not determine Chrome version. Downloading ChromeDriver automatically.")
    
    # Check for existing ChromeDriver in the custom path
    existing_driver = find_existing_driver(custom_path, chrome_version)
    
    if existing_driver:
        print("Using existing ChromeDriver.")
        driver_path = existing_driver
    else:
        print("No matching ChromeDriver found. Downloading new driver.")
        driver_path = download_and_move_driver(custom_path)

    return driver_path, False  # Return the ChromeDriver path and flag it as ChromeDriver

# Step 6: Updated tree_session and snap_tree functions with Xvfb virtual display
def tree_session(tree, layouts, port):
    # Start the ETE tree exploration session
    tree.explore(layouts=layouts, keep_server=True, open_browser=False, compress=False,  
                show_leaf_name=True, show_branch_length=True, show_branch_support=True, port=port)

def snap_tree(port, plot_file, driver_path, is_chrome_binary):
    # Start a virtual display (for headless environments)
    display = Display(visible=0, size=(1920, 1080))  # Set larger display size
    display.start()
    print("Virtual display started")
    
    try:
        url = 'http://127.0.0.1:{}/'.format(str(port))
        print(f"Opening URL: {url}")
        
        options = Options()
        options.headless = False  # Disable headless mode for debugging, change to True later
        options.add_argument("--no-sandbox")
        options.add_argument("--disable-dev-shm-usage")
        options.add_argument("--disable-gpu")
        #options.add_argument("--remote-debugging-port=9222")  # Required for remote debugging
        options.add_argument("--disable-software-rasterizer")  # Fix for headless environments
        options.add_argument("--remote-debugging-port=9230")

        # Specify the directory where the file will be downloaded
        default_download_dir = os.getcwd()  # Set to current working directory

        # Set Chrome preferences for automatic download to the specific directory
        prefs = {
            "download.default_directory": default_download_dir,  # Set download directory
            "download.prompt_for_download": False,  # Disable download confirmation dialog
            "download.directory_upgrade": True,  # Automatically upgrade download directory
            "safebrowsing.enabled": True  # Enable safe browsing
        }
        options.add_experimental_option("prefs", prefs)

        # If we're using Chrome, set the binary location
        if is_chrome_binary:
            options.binary_location = driver_path

        # Set up Chrome with the driver or binary
        print("Setting up Chrome browser")
        driver = webdriver.Chrome(service=Service(driver_path), options=options)
        driver.get(url)
        print("Page loaded. Current URL:", driver.current_url)

        # Perform actions (ensure the page is focused before sending the 'd' key)
        actions = ActionChains(driver)
        actions.click().send_keys('d').perform()  # Ensure page has focus before sending key
        print("Sent 'd' key to the page")
        
        # Wait for the file to download
        default_download_dir = os.getcwd()  # Specify the directory where the file will be downloaded
        file_name = 'tree-1.svg'  # Specify the name of the downloaded file
        file_path = os.path.join(default_download_dir, file_name)
        print(f"Looking for file at: {file_path}")

        # Wait for the file to appear in the directory
        wait_time = 0
        while not os.path.exists(file_path):
            print("Waiting for the file to download...")
            time.sleep(1)
            wait_time += 1
            if wait_time > 30:  # Timeout after 30 seconds
                print("Timeout waiting for file download.")
                break

        # Check if the file was downloaded and move it
        if os.path.exists(file_path):
            os.rename(file_path, plot_file)
            print(f"File downloaded and moved to: {plot_file}")
        else:
            print("File was not downloaded.")

        driver.quit()

    except Exception as e:
        print(f"Error occurred: {e}")
    
    finally:
        print("Stopping virtual display")
        display.stop()  # Stop the virtual display when done

# Step 7: Get image function using multiprocessing
def get_image(tree, layouts, port, plot_file):
    driver_path, is_chrome_binary = setup_driver(DRIVER_DIR)  # Set up ChromeDriver or Chrome binary
    #driver_path = "/home/deng/Projects/metatree_drawer/test_selenium/drivers/chromedriver"
    #is_chrome_binary = False
    p1 = Process(target=tree_session, args=(tree, layouts, port,))
    p2 = Process(target=snap_tree, args=(port, plot_file, driver_path, is_chrome_binary,))
    
    p1.start()
    time.sleep(4)  # Give the tree session some time to initialize
    p2.start()
    time.sleep(4)  # Wait for the snapshot process
    p1.kill()
    p2.kill()