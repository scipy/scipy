#! /bin/sh

if test -n "${XML_CATALOG_FILES:-}"; then
    xml_catalog_files_libxml2="${XML_CATALOG_FILES}"
    XML_CATALOG_FILES="${XML_CATALOG_FILES} "
else
    xml_catalog_files_libxml2=""
    XML_CATALOG_FILES=""
fi


# Replace space with '%20'; equivalent to
# conda_catalog_files=${CONDA_PREFIX// /%20}, except trailing space is
# ignored.
conda_catalog_files=""
ifs_libxml2="${IFS}"
IFS=" "
rem="${CONDA_PREFIX}"
for pre in ${rem}; do
    while test "${rem#"${pre}"}" = "${rem}"; do
	conda_catalog_files="${conda_catalog_files}%20"
	rem=${rem#" "}
    done
    conda_catalog_files="${conda_catalog_files}${pre}"
    rem=${rem#"${pre}"}
done
IFS="${ifs_libxml2}"

conda_catalog_files="file://${conda_catalog_files}/etc/xml/catalog file:///etc/xml/catalog"
export XML_CATALOG_FILES="${XML_CATALOG_FILES}${conda_catalog_files}"
unset conda_catalog_files ifs_libxml2 rem
