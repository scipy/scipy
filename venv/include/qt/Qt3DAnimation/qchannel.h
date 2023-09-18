/****************************************************************************
**
** Copyright (C) 2017 Klaralvdalens Datakonsult AB (KDAB).
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
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

#ifndef QT3DANIMATION_QCHANNEL_H
#define QT3DANIMATION_QCHANNEL_H

#include <QtCore/qscopedpointer.h>
#include <QtCore/qstring.h>
#include <Qt3DAnimation/qt3danimation_global.h>
#include <Qt3DAnimation/qchannelcomponent.h>

QT_BEGIN_NAMESPACE

namespace Qt3DAnimation {

class QChannelPrivate;

class Q_3DANIMATIONSHARED_EXPORT QChannel
{
public:
    QChannel();
    explicit QChannel(const QString &name);
    QChannel(const QChannel &);
    QChannel &operator=(const QChannel &);
    ~QChannel();

    void setName(const QString &name);
    QString name() const;

    void setJointIndex(int jointIndex);
    int jointIndex() const;

    int channelComponentCount() const;
    void appendChannelComponent(const QChannelComponent &component);
    void insertChannelComponent(int index, const QChannelComponent &component);
    void removeChannelComponent(int index);
    void clearChannelComponents();

    // Iterator API
    typedef const QChannelComponent *const_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    const_iterator begin()  const Q_DECL_NOTHROW;
    const_iterator cbegin() const Q_DECL_NOTHROW { return begin(); }
    const_iterator end()    const Q_DECL_NOTHROW;
    const_iterator cend()   const Q_DECL_NOTHROW { return end(); }
    const_reverse_iterator rbegin()  const Q_DECL_NOTHROW { return const_reverse_iterator(end()); }
    const_reverse_iterator crbegin() const Q_DECL_NOTHROW { return rbegin(); }
    const_reverse_iterator rend()    const Q_DECL_NOTHROW { return const_reverse_iterator(begin()); }
    const_reverse_iterator crend()   const Q_DECL_NOTHROW { return rend(); }

    friend bool operator==(const QChannel &,
                           const QChannel &) Q_DECL_NOTHROW;
    friend bool operator!=(const QChannel &,
                           const QChannel &) Q_DECL_NOTHROW;

private:
    QScopedPointer<QChannelPrivate> d;
};

} // namespace Qt3DAnimation

QT_END_NAMESPACE

#endif // QT3DANIMATION_QCHANNEL_H
