/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Data Visualization module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QDATAVISUALIZATIONGLOBAL_H
#define QDATAVISUALIZATIONGLOBAL_H

#include <QtCore/qglobal.h>

#define QT_DATAVISUALIZATION_VERSION_STR QT_VERSION_STR
/*
   QT_DATAVISUALIZATION_VERSION is (major << 16) + (minor << 8) + patch.
*/
#define QT_DATAVISUALIZATION_VERSION QT_VERSION
/*
   can be used like #if (QT_DATAVISUALIZATION_VERSION >= QT_DATAVISUALIZATION_VERSION_CHECK(1, 0, 0))
*/
#define QT_DATAVISUALIZATION_VERSION_CHECK(major, minor, patch) ((major<<16)|(minor<<8)|(patch))

#ifndef QT_STATIC
#  if defined(QT_BUILD_DATAVISUALIZATION_LIB)
#    define QT_DATAVISUALIZATION_EXPORT Q_DECL_EXPORT
#  else
#    define QT_DATAVISUALIZATION_EXPORT Q_DECL_IMPORT
#  endif
#else
#  define QT_DATAVISUALIZATION_EXPORT
#endif

#ifndef Q_CLANG_QDOC
#define QT_BEGIN_NAMESPACE_DATAVISUALIZATION namespace QtDataVisualization {
#define QT_END_NAMESPACE_DATAVISUALIZATION }
#else /* Let documentation be generated with the standard Qt namespace */
#define QT_BEGIN_NAMESPACE_DATAVISUALIZATION QT_BEGIN_NAMESPACE
#define QT_END_NAMESPACE_DATAVISUALIZATION QT_END_NAMESPACE
#endif // Q_CLANG_QDOC

#endif
