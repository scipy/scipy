import _zsuperlu
_memory_dict = _zsuperlu._memory_dict
del _zsuperlu
from _zsuperlu import *
from _ssuperlu import *
from _dsuperlu import *
from _csuperlu import *

