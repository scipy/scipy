/****************************************************************************
**
** Copyright (C) 2014 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DCORE_QSCENECHANGE_H
#define QT3DCORE_QSCENECHANGE_H

#include <Qt3DCore/qnodeid.h>
#include <Qt3DCore/qt3dcore_global.h>
#include <QtCore/QSharedPointer>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

enum Q3D_DECL_DEPRECATED ChangeFlag {
    NodeCreated             = 1 << 0,
    NodeDeleted             = 1 << 1,
    PropertyUpdated         = 1 << 2,
    PropertyValueAdded      = 1 << 3,
    PropertyValueRemoved    = 1 << 4,
    ComponentAdded          = 1 << 5,
    ComponentRemoved        = 1 << 6,
    CommandRequested        = 1 << 7,
    CallbackTriggered       = 1 << 8,
    AllChanges              = 0xFFFFFFFF
};

QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
Q_DECLARE_FLAGS(ChangeFlags, ChangeFlag)
Q_DECLARE_OPERATORS_FOR_FLAGS(ChangeFlags)
QT_WARNING_POP

class QNode;
class QSceneChangePrivate;

class Q_3DCORESHARED_EXPORT QSceneChange
{
public:
    enum DeliveryFlag {
        BackendNodes = 0x0001,
        Nodes = 0x0010,
        DeliverToAll = BackendNodes | Nodes
    };
    Q_DECLARE_FLAGS(DeliveryFlags, DeliveryFlag)

    virtual ~QSceneChange();

    QT_WARNING_PUSH
    QT_WARNING_DISABLE_DEPRECATED
    ChangeFlag type() const Q_DECL_NOTHROW;
    QT_WARNING_POP

    void setDeliveryFlags(DeliveryFlags flags) Q_DECL_NOTHROW;
    DeliveryFlags deliveryFlags() const Q_DECL_NOTHROW;

    QNodeId subjectId() const Q_DECL_NOTHROW;

protected:
    Q_DECLARE_PRIVATE(QSceneChange)
    QT_WARNING_PUSH
    QT_WARNING_DISABLE_DEPRECATED
    Q3D_DECL_DEPRECATED explicit QSceneChange(ChangeFlag type, QNodeId subjectId);
    Q3D_DECL_DEPRECATED explicit QSceneChange(QSceneChangePrivate &dd,
                                              ChangeFlag type, QNodeId subjectId);
    QT_WARNING_POP
    QSceneChangePrivate *d_ptr;

private:
    Q_DISABLE_COPY(QSceneChange)
};

Q_DECLARE_OPERATORS_FOR_FLAGS(QSceneChange::DeliveryFlags)

typedef QSharedPointer<QSceneChange> QSceneChangePtr;

} // namespace Qt3DCore

QT_END_NAMESPACE

#endif // QT3DCORE_QSCENECHANGE_H
