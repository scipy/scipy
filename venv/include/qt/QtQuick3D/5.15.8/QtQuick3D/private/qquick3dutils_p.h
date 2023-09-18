/****************************************************************************
**
** Copyright (C) 2019 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of Qt Quick 3D.
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

#ifndef QQUICK3DUTILS_P_H
#define QQUICK3DUTILS_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API. It exists purely as an
// implementation detail. This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <type_traits>

#include <QtCore/QtGlobal>

QT_BEGIN_NAMESPACE

// Assigns 'updated' to 'orig' and returns true if they are different
template<typename T, typename std::enable_if<!std::is_floating_point<T>::value, int>::type = 0>
bool qUpdateIfNeeded(T &orig, T updated)
{
    if (orig == updated)
        return false;
    orig = updated;
    return true;
}

// Assigns 'updated' to 'orig' and returns true if they are different, compared with qFuzzyCompare
template <typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
bool qUpdateIfNeeded(T &orig, T updated)
{
    if (qFuzzyCompare(orig, updated))
        return false;
    orig = updated;
    return true;
}

QT_END_NAMESPACE

#endif // QQUICK3DUTILS_P_H
