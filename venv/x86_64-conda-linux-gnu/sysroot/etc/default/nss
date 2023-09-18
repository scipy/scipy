# /etc/default/nss
# This file can theoretically contain a bunch of customization variables
# for Name Service Switch in the GNU C library.  For now there are only
# four variables:
#
# NETID_AUTHORITATIVE
#   If set to TRUE, the initgroups() function will accept the information
#   from the netid.byname NIS map as authoritative.  This can speed up the
#   function significantly if the group.byname map is large.  The content
#   of the netid.byname map is used AS IS.  The system administrator has
#   to make sure it is correctly generated.
#NETID_AUTHORITATIVE=TRUE
#
# SERVICES_AUTHORITATIVE
#   If set to TRUE, the getservbyname{,_r}() function will assume
#   services.byservicename NIS map exists and is authoritative, particularly
#   that it contains both keys with /proto and without /proto for both
#   primary service names and service aliases.  The system administrator
#   has to make sure it is correctly generated.
#SERVICES_AUTHORITATIVE=TRUE
#
# SETENT_BATCH_READ
#  If set to TRUE, various setXXent() functions will read the entire
#  database at once and then hand out the requests one by one from
#  memory with every getXXent() call.  Otherwise each getXXent() call
#  might result into a network communication with the server to get
#  the next entry.
SETENT_BATCH_READ=TRUE
#
# ADJUNCT_AS_SHADOW
#  If set to TRUE, the passwd routines in the NIS NSS module will not
#  use the passwd.adjunct.byname tables to fill in the password data
#  in the passwd structure.  This is a security problem if the NIS
#  server cannot be trusted to send the passwd.adjuct table only to
#  privileged clients.  Instead the passwd.adjunct.byname table is
#  used to synthesize the shadow.byname table if it does not exist.
#ADJUNCT_AS_SHADOW=TRUE
