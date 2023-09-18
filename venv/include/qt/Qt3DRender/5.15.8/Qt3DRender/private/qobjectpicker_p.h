/****************************************************************************
**
** Copyright (C) 2016 Klaralvdalens Datakonsult AB (KDAB).
** Contact: http://www.qt-project.org/legal
**
** This file is part of the Qt3D module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL3$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see http://www.qt.io/terms-conditions. For further
** information use the contact form at http://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPLv3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or later as published by the Free
** Software Foundation and appearing in the file LICENSE.GPL included in
** the packaging of this file. Please review the following information to
** ensure the GNU General Public License version 2.0 requirements will be
** met: http://www.gnu.org/licenses/gpl-2.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QT3DRENDER_QOBJECTPICKER_P_H
#define QT3DRENDER_QOBJECTPICKER_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of other Qt classes.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <Qt3DCore/private/qcomponent_p.h>
#include <Qt3DCore/qnodeid.h>
#include <Qt3DRender/private/qt3drender_global_p.h>
#include <Qt3DRender/private/objectpicker_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

/*!
    \internal
*/
class Q_3DRENDERSHARED_PRIVATE_EXPORT QObjectPickerPrivate : public Qt3DCore::QComponentPrivate
{
public:
    QObjectPickerPrivate()
        : QComponentPrivate()
        , m_hoverEnabled(false)
        , m_dragEnabled(false)
        , m_pressed(false)
        , m_containsMouse(false)
        , m_acceptedLastPressedEvent(true)
        , m_priority(0)
    {
        m_shareable = false;
    }

    Q_DECLARE_PUBLIC(QObjectPicker)
    bool m_hoverEnabled;
    bool m_dragEnabled;
    bool m_pressed;
    bool m_containsMouse;
    bool m_acceptedLastPressedEvent;
    int m_priority;

    enum EventType {
        Pressed,
        Released,
        Clicked,
        Moved
    };

    void propagateEvent(QPickEvent *event, EventType type);

    void pressedEvent(QPickEvent *event);
    void clickedEvent(QPickEvent *event);
    void movedEvent(QPickEvent *event);
    void releasedEvent(QPickEvent *event);

    void setPressed(bool pressed);
    void setContainsMouse(bool containsMouse);
};

struct QObjectPickerData
{
    bool hoverEnabled;
    bool dragEnabled;
    int priority;
};

struct QObjectPickerEvent
{
    QPickEventPtr event;
    Qt3DCore::QNodeId viewportNodeId;
};

} // namespace Qt3DRender

QT_END_NAMESPACE

Q_DECLARE_METATYPE(Qt3DRender::QObjectPickerEvent);

#endif // QT3DRENDER_QOBJECTPICKER_P_H
