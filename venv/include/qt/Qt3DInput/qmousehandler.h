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

#ifndef QT3DINPUT_QMOUSEHANDLER_H
#define QT3DINPUT_QMOUSEHANDLER_H

#include <Qt3DInput/qt3dinput_global.h>
#include <Qt3DInput/qmouseevent.h>
#include <Qt3DCore/qcomponent.h>

QT_BEGIN_NAMESPACE

namespace Qt3DInput {

class QMouseDevice;
class QMouseHandlerPrivate;

class Q_3DINPUTSHARED_EXPORT QMouseHandler : public Qt3DCore::QComponent
{
    Q_OBJECT
    Q_PROPERTY(Qt3DInput::QMouseDevice *sourceDevice READ sourceDevice WRITE setSourceDevice NOTIFY sourceDeviceChanged)
    Q_PROPERTY(bool containsMouse READ containsMouse NOTIFY containsMouseChanged)

public:
    explicit QMouseHandler(QNode *parent = nullptr);
    ~QMouseHandler();

    QMouseDevice *sourceDevice() const;
    bool containsMouse() const;

public Q_SLOTS:
    void setSourceDevice(QMouseDevice *mouseDevice);

Q_SIGNALS:
    void sourceDeviceChanged(QMouseDevice *mouseDevice);
    void containsMouseChanged(bool containsMouse);

    void clicked(Qt3DInput::QMouseEvent *mouse);
    void doubleClicked(Qt3DInput::QMouseEvent *mouse);

    void entered();
    void exited();

    void pressed(Qt3DInput::QMouseEvent *mouse);
    void released(Qt3DInput::QMouseEvent *mouse);
    void pressAndHold(Qt3DInput::QMouseEvent *mouse);
    void positionChanged(Qt3DInput::QMouseEvent *mouse);

#if QT_CONFIG(wheelevent)
    void wheel(Qt3DInput::QWheelEvent *wheel);
#endif

protected:
    // TODO Unused remove in Qt6
    void sceneChangeEvent(const Qt3DCore::QSceneChangePtr &change) override;
    void setContainsMouse(bool contains);

private:
    Q_DECLARE_PRIVATE(QMouseHandler)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

} // namespace Qt3DInput

QT_END_NAMESPACE

#endif // QT3DINPUT_QMOUSEHANDLER_H
