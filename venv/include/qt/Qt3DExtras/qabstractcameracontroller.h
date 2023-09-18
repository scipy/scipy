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

#ifndef QT3DEXTRAS_QABSTRACTCAMERACONTROLLER_H
#define QT3DEXTRAS_QABSTRACTCAMERACONTROLLER_H

#include <Qt3DCore/QEntity>
#include <Qt3DExtras/qt3dextras_global.h>

QT_BEGIN_NAMESPACE

namespace Qt3DInput {
class QKeyboardDevice;
class QMouseDevice;
}

namespace Qt3DRender {
class QCamera;
}

namespace Qt3DExtras {

class QAbstractCameraControllerPrivate;

class Q_3DEXTRASSHARED_EXPORT QAbstractCameraController : public Qt3DCore::QEntity
{
    Q_OBJECT
    Q_PROPERTY(Qt3DRender::QCamera *camera READ camera WRITE setCamera NOTIFY cameraChanged)
    Q_PROPERTY(float linearSpeed READ linearSpeed WRITE setLinearSpeed NOTIFY linearSpeedChanged)
    Q_PROPERTY(float lookSpeed READ lookSpeed WRITE setLookSpeed NOTIFY lookSpeedChanged)
    Q_PROPERTY(float acceleration READ acceleration WRITE setAcceleration NOTIFY accelerationChanged)
    Q_PROPERTY(float deceleration READ deceleration WRITE setDeceleration NOTIFY decelerationChanged)

public:
    ~QAbstractCameraController();

    Qt3DRender::QCamera *camera() const;
    float linearSpeed() const;
    float lookSpeed() const;

    float acceleration() const;
    float deceleration() const;

    void setCamera(Qt3DRender::QCamera *camera);
    void setLinearSpeed(float linearSpeed);
    void setLookSpeed(float lookSpeed);

    void setAcceleration(float acceleration);
    void setDeceleration(float deceleration);

Q_SIGNALS:
    void cameraChanged();
    void linearSpeedChanged();
    void lookSpeedChanged();

    void accelerationChanged(float acceleration);
    void decelerationChanged(float deceleration);

protected:
    explicit QAbstractCameraController(Qt3DCore::QNode *parent = nullptr);
    QAbstractCameraController(QAbstractCameraControllerPrivate &dd, Qt3DCore::QNode *parent = nullptr);

    Qt3DInput::QKeyboardDevice *keyboardDevice() const;
    Qt3DInput::QMouseDevice *mouseDevice() const;

public:
    struct InputState
    {
        float rxAxisValue;
        float ryAxisValue;
        float txAxisValue;
        float tyAxisValue;
        float tzAxisValue;

        bool leftMouseButtonActive;
        bool middleMouseButtonActive;
        bool rightMouseButtonActive;

        bool altKeyActive;
        bool shiftKeyActive;
    };

private:
    virtual void moveCamera(const InputState &state, float dt) = 0;

private:
    Q_DECLARE_PRIVATE(QAbstractCameraController)
};

}   // Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DEXTRAS_QABSTRACTCAMERACONTROLLER_H
