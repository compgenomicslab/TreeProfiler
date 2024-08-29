import os
import glob

def replace_static_path(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        content = file.read()

    # Replace occurrences of '_static' with 'static'
    content = content.replace('_static', 'static')

    with open(file_path, 'w', encoding='utf-8') as file:
        file.write(content)

def process_html_files(directory):
    html_files = glob.glob(os.path.join(directory, '**/*.html'), recursive=True)
    for html_file in html_files:
        replace_static_path(html_file)

if __name__ == "__main__":
    html_dir = os.path.join('build', 'html')  # Adjust this if your build directory is different
    process_html_files(html_dir)
