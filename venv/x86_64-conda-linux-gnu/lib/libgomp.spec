# This spec file is read by gcc when linking.  It is used to specify the
# standard libraries we need in order to link with libgomp.
*link_gomp: -lgomp %{static: -ldl }
