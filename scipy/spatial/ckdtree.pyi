from typing import Union, Optional
import numpy as np

class cKDTree:
    def __init__(self,
                 data: np.ndarray,
                 leafsize: int = ...,
                 compact_nodes: bool = ...,
                 copy_data: bool = ...,
                 balanced_tree: bool = ...,
                 boxsize: Optional[Union[np.ndarray, float]] = ...) -> None: ...

    def query_pairs(self,
                    r: float,
                    p: float = ...,
                    eps: float = ...,
                    output_type: str = ...) -> Union[set, np.ndarray]: ...
