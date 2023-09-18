/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtXmlPatterns module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl-3.0.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or (at your option) the GNU General
** Public license version 3 or any later version approved by the KDE Free
** Qt Foundation. The licenses are as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL2 and LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-2.0.html and
** https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#ifndef Patternist_Debug_H
#define Patternist_Debug_H

#include <QtDebug>
/**
 * @file
 * @short Contains macros for debugging.
 */

QT_BEGIN_NAMESPACE

/**
 * @short Enables detailed parser debug output.
 *
 * If this macro is defined, @em a @em lot of debugging information will be outputted.
 * This is all the state transitions, token shifting, and rule reductions that
 * the parser do.
 *
 * This is automatically disabled if @c QT_NO_DEBUG is defined.
 */
#define Patternist_DEBUG_PARSER

/**
 * @short Enables debug printing statements.
 *
 * Patternist does not use qDebug(), but pDebug() instead. It only output
 * if this define is defined.
 *
 * It is automatically disabled if @c QT_NO_DEBUG is defined.
 */
#define Patternist_DEBUG

#undef Patternist_DEBUG // disable it for now

#ifdef QT_NO_DEBUG
#   undef Patternist_DEBUG_PARSER
#   undef Patternist_DEBUG
#endif

namespace QPatternist
{
#ifdef Patternist_DEBUG
    inline QDebug pDebug()
    {
        return qDebug();
    }
#else
    inline QNoDebug pDebug()
    {
        return QNoDebug();
    }
#endif
}

QT_END_NAMESPACE
#endif
