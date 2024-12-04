import glob
import os
import shutil
import jupytext

def call_jupytext(md_files, _contents_path, _contents_cache_path):
    is_cached = os.path.exists(_contents_cache_path)
    if not is_cached:
        for md_file in md_files:
            basename = os.path.splitext(os.path.basename(md_file))[0]
            output_name = os.path.join(_contents_path, f"{basename}.ipynb")
            nb = jupytext.read(md_file)
            jupytext.write(jupytext.read(md_file), output_name, fmt="ipynb", version=4)
        return True

    is_dirty = False

    for md_file in md_files:
        basename = os.path.splitext(os.path.basename(md_file))[0]
        output_path = os.path.join(_contents_path, f"{basename}.ipynb")
        cached_output_path = os.path.join(_contents_cache_path, f"{basename}.ipynb")
        cmd_execution_time = os.stat(cached_output_path).st_mtime
        md_file_modification_time = os.stat(md_file).st_mtime
        if cmd_execution_time < md_file_modification_time:
            nb = jupytext.read(md_file)
            jupytext.write(nb, output_path, fmt="ipynb", version=4)
            is_dirty = True
        else:
            shutil.copyfile(
                cached_output_path,
                os.path.join(_contents_path, f"{basename}.ipynb")
            )

    return is_dirty

if __name__ == '__main__':
    _contents_cache_path = os.path.join("build", "_contents")
    _contents_path = os.path.join("source", "_contents")

    os.makedirs(os.path.expanduser(_contents_path), exist_ok=True)

    md_files = glob.glob("source/tutorial/stats/*.md")
    is_dirty = call_jupytext(
        md_files,
        _contents_path,
        _contents_cache_path
    )

    if is_dirty:
        os.makedirs(os.path.expanduser(_contents_cache_path), exist_ok=True)
        shutil.copytree(_contents_path, _contents_cache_path, dirs_exist_ok=True)
