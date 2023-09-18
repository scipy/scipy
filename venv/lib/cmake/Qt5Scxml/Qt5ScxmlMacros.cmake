#
# Copyright (C) 2016 Klaralvdalens Datakonsult AB (KDAB).
# Contact: https://www.qt.io/licensing/
#
# This file is part of the QtScxml module of the Qt Toolkit.
#
# $QT_BEGIN_LICENSE:LGPL$
# Commercial License Usage
# Licensees holding valid commercial Qt licenses may use this file in
# accordance with the commercial license agreement provided with the
# Software or, alternatively, in accordance with the terms contained in
# a written agreement between you and The Qt Company. For licensing terms
# and conditions see https://www.qt.io/terms-conditions. For further
# information use the contact form at https://www.qt.io/contact-us.
#
# GNU Lesser General Public License Usage
# Alternatively, this file may be used under the terms of the GNU Lesser
# General Public License version 3 as published by the Free Software
# Foundation and appearing in the file LICENSE.LGPL3 included in the
# packaging of this file. Please review the following information to
# ensure the GNU Lesser General Public License version 3 requirements
# will be met: https://www.gnu.org/licenses/lgpl-3.0.html.
#
# GNU General Public License Usage
# Alternatively, this file may be used under the terms of the GNU
# General Public License version 2.0 or (at your option) the GNU General
# Public license version 3 or any later version approved by the KDE Free
# Qt Foundation. The licenses are as published by the Free Software
# Foundation and appearing in the file LICENSE.GPL2 and LICENSE.GPL3
# included in the packaging of this file. Please review the following
# information to ensure the GNU General Public License requirements will
# be met: https://www.gnu.org/licenses/gpl-2.0.html and
# https://www.gnu.org/licenses/gpl-3.0.html.
#
# $QT_END_LICENSE$

if(NOT Qt5Scxml_QSCXMLC_EXECUTABLE)
    message(FATAL_ERROR "qscxmlc executable not found -- Check installation.")
endif()

# qt5_add_statecharts(outfiles inputfile ... )

function(qt5_add_statecharts outfiles)
    set(options)
    set(oneValueArgs)
    set(multiValueArgs OPTIONS)

    cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    set(scxml_files ${ARGS_UNPARSED_ARGUMENTS})

    foreach(it ${scxml_files})
        get_filename_component(outfilename ${it} NAME_WE)
        get_filename_component(infile ${it} ABSOLUTE)
        set(outfile ${CMAKE_CURRENT_BINARY_DIR}/${outfilename})
        set(outfile_cpp ${CMAKE_CURRENT_BINARY_DIR}/${outfilename}.cpp)
        set(outfile_h ${CMAKE_CURRENT_BINARY_DIR}/${outfilename}.h)

        add_custom_command(OUTPUT ${outfile_cpp} ${outfile_h}
                           COMMAND ${Qt5Scxml_QSCXMLC_EXECUTABLE}
                           ARGS ${ARGS_OPTIONS} --output ${outfile} ${infile}
                           MAIN_DEPENDENCY ${infile}
                           VERBATIM)
        list(APPEND ${outfiles} ${outfile_cpp})
    endforeach()
    set_source_files_properties(${outfiles} PROPERTIES SKIP_AUTOMOC TRUE)
    set(${outfiles} ${${outfiles}} PARENT_SCOPE)
endfunction()

if(NOT QT_NO_CREATE_VERSIONLESS_FUNCTIONS)
    function(qt_add_statecharts outfiles)
        if(QT_DEFAULT_MAJOR_VERSION EQUAL 5)
            qt5_add_statecharts("${outfiles}" ${ARGN})
        elseif(QT_DEFAULT_MAJOR_VERSION EQUAL 6)
            qt6_add_statecharts("${outfiles}" ${ARGN})
        endif()
        set("${outfiles}" "${${outfiles}}" PARENT_SCOPE)
    endfunction()
endif()
