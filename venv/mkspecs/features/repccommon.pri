# Detect repc when QtRO is installed into non-Qt prefix
cmd = $${QT.remoteobjects.bins}/repc
contains(QMAKE_HOST.os, Windows) {
    cmd = $$system_path($${cmd}.exe)
}
exists($$cmd): QT_TOOL.repc.binary = $$cmd

# qtPrepareTool honors QT_TOOL.repc.binary if set
qtPrepareTool(QMAKE_REPC, repc)

REPC_INCLUDEPATHES = $$QT.remoteobjects.includes
for (path, REPC_INCLUDEPATHES) {
    REPC_INCLUDEPATH += -I $$path
}

isEmpty(QMAKE_MOD_REPC):QMAKE_MOD_REPC = rep_

repc_TYPE = $$upper($$repc_type)

load(moc)

groups =
for(entry, REPC_$$repc_TYPE) {
    files = $$eval($${entry}.files)
    isEmpty(files) {
        files = $$entry
        group = repc_$${repc_type}
    } else {
        group = $${entry}_repc_$${repc_type}
    }
    groups *= $$group

    input_list = $$upper($$group)_LIST
    for(subent, $$list($$unique(files))) {
        $$input_list += $$subent

        # Add directory of *.rep file to include path
        file_path = $$_PRO_FILE_PWD_/$$subent
        INCLUDEPATH *= $$dirname(file_path)
    }
}

for(group, groups) {
    GROUP = $$upper($$group)
    input_list = $${GROUP}_LIST

    $${group}_header.output  = $$QMAKE_MOD_REPC${QMAKE_FILE_BASE}_$${repc_type}.h
    $${group}_header.commands = $$QMAKE_REPC $$repc_option $$REPC_INCLUDEPATH ${QMAKE_FILE_NAME} ${QMAKE_FILE_OUT}
    $${group}_header.depends = ${QMAKE_FILE_NAME} $$QT_TOOL.repc.binary
    $${group}_header.variable_out = $${GROUP}_HEADERS
    $${group}_header.input = $$input_list

    $${group}_moc.commands = $$moc_header.commands $$REPC_INCLUDEPATH
    $${group}_moc.output = $$moc_header.output
    $${group}_moc.input = $${GROUP}_HEADERS
    $${group}_moc.variable_out = GENERATED_SOURCES
    !contains(TEMPLATE, vc.*): \
        $${group}_moc.name = $$moc_header.name

    QMAKE_EXTRA_COMPILERS += $${group}_header $${group}_moc
}
