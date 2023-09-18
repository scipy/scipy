/****************************************************************************
**
** Copyright (C) 2020 Klarälvdalens Datakonsult AB, a KDAB Group company, info@kdab.com, author Giuseppe D'Angelo <giuseppe.dangelo@kdab.com>
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtQuick module of the Qt Toolkit.
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

#ifndef QQUICKPARTICLEFLATSET_P_H
#define QQUICKPARTICLEFLATSET_P_H

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

#include <QtGlobal>

#include <vector>
#include <algorithm>
#include <iterator>

QT_BEGIN_NAMESPACE

// Minimal API, just for the consumption of Qt Quick Particles.
// For extra safety, it's in a private namespace

namespace QtQuickParticlesPrivate {

template <typename T>
class QFlatSet
{
public:
    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;
    using value_type = typename std::vector<T>::value_type;
    using size_type = int;

    iterator find(const T &t)
    {
        return std::find(begin(), end(), t);
    }

    const_iterator find(const T &t) const
    {
        return std::find(begin(), end(), t);
    }

    bool contains(const T &t) const
    {
        return find(t) != end();
    }

    void clear()
    {
        m_data.clear();
    }

    void reserve(int capacity)
    {
        m_data.reserve(capacity);
    }

    iterator insert(const T &t)
    {
        auto i = find(t);
        if (i != end())
            return i;
        T copy = t;
        m_data.push_back(std::move(copy));
        return std::prev(m_data.end());
    }

    iterator insert(T &&t)
    {
        auto i = find(t);
        if (i != end())
            return i;
        m_data.push_back(std::move(t));
        return std::prev(m_data.end());
    }

    size_type remove(const T &t)
    {
        auto i = std::find(m_data.begin(), m_data.end(), t);
        if (i != m_data.end()) {
            m_data.erase(i);
            return 1;
        }
        return 0;
    }

    iterator operator<<(const T &t)
    {
        return insert(t);
    }

    iterator operator<<(T &&t)
    {
        return insert(std::move(t));
    }

    iterator begin() { return m_data.begin(); }
    const_iterator begin() const { return m_data.begin(); }
    const_iterator cbegin() const { return m_data.cbegin(); }

    iterator end() { return m_data.end(); }
    const_iterator end() const { return m_data.end(); }
    const_iterator cend() const { return m_data.cend(); }

private:
    std::vector<T> m_data;
};

} // namespace QtQuickParticlesPrivate

QT_END_NAMESPACE

#endif // QQUICKPARTICLEFLATSET_P_H
