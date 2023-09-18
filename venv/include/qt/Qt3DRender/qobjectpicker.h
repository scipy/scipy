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

#ifndef QT3DRENDER_QOBJECTPICKER_H
#define QT3DRENDER_QOBJECTPICKER_H

#include <Qt3DCore/qcomponent.h>
#include <Qt3DRender/qt3drender_global.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QAttribute;
class QObjectPickerPrivate;
class QPickEvent;

class Q_3DRENDERSHARED_EXPORT QObjectPicker : public Qt3DCore::QComponent
{
    Q_OBJECT
    Q_PROPERTY(bool hoverEnabled READ isHoverEnabled WRITE setHoverEnabled NOTIFY hoverEnabledChanged)
    Q_PROPERTY(bool dragEnabled READ isDragEnabled WRITE setDragEnabled NOTIFY dragEnabledChanged)
    Q_PROPERTY(bool pressed READ isPressed NOTIFY pressedChanged)
    Q_PROPERTY(bool containsMouse READ containsMouse NOTIFY containsMouseChanged)
    Q_PROPERTY(int priority READ priority WRITE setPriority NOTIFY priorityChanged REVISION 13)

public:
    explicit QObjectPicker(QNode *parent = nullptr);
    ~QObjectPicker();

    bool isHoverEnabled() const;
    bool isDragEnabled() const;

    bool containsMouse() const;
    bool isPressed() const;

    int priority() const;

public Q_SLOTS:
    void setHoverEnabled(bool hoverEnabled);
    void setDragEnabled(bool dragEnabled);
    Q_REVISION(13) void setPriority(int priority);

protected:
    // TODO Unused remove in Qt6
    void sceneChangeEvent(const Qt3DCore::QSceneChangePtr &change) override;

Q_SIGNALS:
    void pressed(Qt3DRender::QPickEvent *pick);
    void released(Qt3DRender::QPickEvent *pick);
    void clicked(Qt3DRender::QPickEvent *pick);
    void moved(Qt3DRender::QPickEvent *pick);
    void entered();
    void exited();
    void hoverEnabledChanged(bool hoverEnabled);
    void dragEnabledChanged(bool dragEnabled);
    void pressedChanged(bool pressed);
    void containsMouseChanged(bool containsMouse);
    Q_REVISION(13) void priorityChanged(int priority);

private:
    Q_DECLARE_PRIVATE(QObjectPicker)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

} // Qt3D

QT_END_NAMESPACE

#endif // QT3DRENDER_QOBJECTPICKER_H
