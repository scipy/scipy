/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtQml module of the Qt Toolkit.
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

#ifndef QTQMLGLOBAL_P_H
#define QTQMLGLOBAL_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <QtCore/private/qglobal_p.h>
#include <QtQml/qtqmlglobal.h>
#include <QtQml/qqmlextensionplugin.h>
#ifndef QT_QML_BOOTSTRAPPED
#  include <QtQml/private/qtqml-config_p.h>
#endif
#include <private/qqmlapiversion_p.h>

#define Q_QML_PRIVATE_EXPORT Q_QML_EXPORT

void Q_QML_PRIVATE_EXPORT qml_register_types_QtQml();
GHS_KEEP_REFERENCE(qml_register_types_QtQml);

#if !defined(QT_QMLDEVTOOLS_LIB) && !defined(QT_BUILD_QMLDEVTOOLS_LIB)
#  define Q_QML_AUTOTEST_EXPORT Q_AUTOTEST_EXPORT
#else
#  define Q_QML_AUTOTEST_EXPORT
#endif

// When doing macOS universal builds, JIT needs to be disabled for the ARM slice.
// Because both arm and x86_64 slices are built in one clang frontend invocation
// we need this hack to ensure each backend invocation sees the correct value
// of the feature definition.

// Unset dummy value
#undef QT_QML_JIT_SUPPORTED_IMPL
// Compute per-arch value and save in extra define
#if QT_CONFIG(qml_jit) && !(defined(Q_OS_MACOS) && defined(Q_PROCESSOR_ARM))
#  define QT_QML_JIT_SUPPORTED_IMPL 1
#else
#  define QT_QML_JIT_SUPPORTED_IMPL 0
#endif
// Unset original feature value
#undef QT_FEATURE_qml_jit
// Set new value based on previous computation
#if QT_QML_JIT_SUPPORTED_IMPL
#  define QT_FEATURE_qml_jit 1
#else
#  define QT_FEATURE_qml_jit -1
#endif

#endif // QTQMLGLOBAL_P_H
