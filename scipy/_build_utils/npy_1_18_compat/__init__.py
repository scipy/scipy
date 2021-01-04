import pathlib

def get_include() -> str:
    return str(pathlib.Path(__file__).parent)
