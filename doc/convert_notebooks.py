import glob
import os
import shutil

def call_jupytext(md_files, _contents_path, _contents_cache_path):
    is_cached = os.path.exists(_contents_cache_path)
    if not is_cached:
        for md_file in md_files:
            basename = os.path.splitext(os.path.basename(md_file))[0]
            output_name = os.path.join(_contents_path, f"{basename}.ipynb")
            os.system(f"python3 -m jupytext --output {output_name} {md_file}")
        return True

    is_dirty = False
    for md_file in md_files:
        basename = os.path.splitext(os.path.basename(md_file))[0]
        output_path = os.path.join(_contents_path, f"{basename}.ipynb")
        cached_output_path = os.path.join(_contents_cache_path, f"{basename}.ipynb")
        cmd_execution_time = os.path.getctime(cached_output_path)
        md_file_modification_time = os.path.getmtime(md_file)
        if cmd_execution_time <= md_file_modification_time:
            os.system(f"python3 -m jupytext --output {output_path} {md_file}")
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
