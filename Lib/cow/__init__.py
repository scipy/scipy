#
# cow - Cluster Of Workstations
#

import os as _os
execfile(_os.path.join(__path__[0],'pre___init__.py'),globals(),locals())

import sync_cluster
from cow import *

Mrange = range(1,17)    
import cow # cow stands for "Cluster Of Workstations"
server_list = []
server_list.append((sync_cluster.host,10000))

local = cow.machine_cluster(server_list)
