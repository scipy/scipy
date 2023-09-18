#! /bin/sh

if test -n "${xml_catalog_files_libxml2:-}"; then
    export XML_CATALOG_FILES="${xml_catalog_files_libxml2}"
else
    unset XML_CATALOG_FILES
fi
unset xml_catalog_files_libxml2
