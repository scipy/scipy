
from twisted.application import service
from twisted.application import internet

from nevow import appserver
import sys
sys.path.append('.')
from livedocs import LiveDocs

application = service.Application("livedocs", uid=65534, gid=65534)
import os

#  A Hack to set the place for temporary variables correctly
os.environ['HOME'] = '/home/nobody'

internet.TCPServer(
    81,
    appserver.NevowSite(
        LiveDocs()
    )
).setServiceParent(application)

