from typing import Union, Optional
import numpy as np

class cKDTree:
    def __init__(self,
                 data: np.ndarray,
                 leafsize: int = 16,
                 compact_nodes: bool = True,
                 copy_data: bool = False,
                 balanced_tree: bool = True,
                 boxsize: Optional[Union[np.ndarray, float]] = None) -> None: ...

    def query_pairs(self,
                    r: float,
                    p: float = 2.,
                    eps: float = 0,
                    output_type: str = 'set') -> Union[set, np.ndarray]: ...
