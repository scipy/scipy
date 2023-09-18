/****************************************************************************
**
** Copyright (C) 2015 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DRENDER_QPICKLINEEVENT_H
#define QT3DRENDER_QPICKLINEEVENT_H

#include <Qt3DRender/qpickevent.h>
#include <QtCore/qsharedpointer.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {
class QPickLineEventPrivate;

class Q_3DRENDERSHARED_EXPORT QPickLineEvent : public QPickEvent
{
    Q_OBJECT
    Q_PROPERTY(uint edgeIndex READ edgeIndex CONSTANT)
    Q_PROPERTY(uint vertex1Index READ vertex1Index CONSTANT)
    Q_PROPERTY(uint vertex2Index READ vertex2Index CONSTANT)
public:
    QPickLineEvent();
    QPickLineEvent(const QPointF &position, const QVector3D& worldIntersection, const QVector3D& localIntersection, float distance,
                   uint edgeIndex, uint vertex1Index, uint vertex2Index, Buttons button, int buttons, int modifiers);
    ~QPickLineEvent();

public:
    uint edgeIndex() const;
    uint vertex1Index() const;
    uint vertex2Index() const;

private:
    Q_DECLARE_PRIVATE(QPickLineEvent)
};

typedef QSharedPointer<QPickLineEvent> QPickLineEventPtr;

} // Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_QPICKLINEEVENT_H
