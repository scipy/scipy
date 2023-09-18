/****************************************************************************
**
** Copyright (C) 2014 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DEXTRAS_QABSTRACTCAMERACONTROLLER_P_H
#define QT3DEXTRAS_QABSTRACTCAMERACONTROLLER_P_H

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

#include <Qt3DExtras/qabstractcameracontroller.h>
#include <QtGui/QVector3D>

#include <Qt3DCore/private/qentity_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {
class QCamera;
}

namespace Qt3DLogic {
class QFrameAction;
}

namespace Qt3DInput {
class QKeyboardDevice;
class QMouseDevice;
class QLogicalDevice;
class QAction;
class QActionInput;
class QAxis;
class QAnalogAxisInput;
class QButtonAxisInput;
class QAxisActionHandler;
}

namespace Qt3DExtras {

class QAbstractCameraControllerPrivate : public Qt3DCore::QEntityPrivate
{
public:
    QAbstractCameraControllerPrivate();

    void init();
    void applyInputAccelerations();

    Qt3DRender::QCamera *m_camera;

    Qt3DInput::QAction *m_leftMouseButtonAction;
    Qt3DInput::QAction *m_middleMouseButtonAction;
    Qt3DInput::QAction *m_rightMouseButtonAction;
    Qt3DInput::QAction *m_altButtonAction;
    Qt3DInput::QAction *m_shiftButtonAction;
    Qt3DInput::QAction *m_escapeButtonAction;

    Qt3DInput::QAxis *m_rxAxis;
    Qt3DInput::QAxis *m_ryAxis;
    Qt3DInput::QAxis *m_txAxis;
    Qt3DInput::QAxis *m_tyAxis;
    Qt3DInput::QAxis *m_tzAxis;

    Qt3DInput::QActionInput *m_leftMouseButtonInput;
    Qt3DInput::QActionInput *m_middleMouseButtonInput;
    Qt3DInput::QActionInput *m_rightMouseButtonInput;
    Qt3DInput::QActionInput *m_altButtonInput;
    Qt3DInput::QActionInput *m_shiftButtonInput;
    Qt3DInput::QActionInput *m_escapeButtonInput;

    Qt3DInput::QAnalogAxisInput *m_mouseRxInput;
    Qt3DInput::QAnalogAxisInput *m_mouseRyInput;
    Qt3DInput::QAnalogAxisInput *m_mouseTzXInput;
    Qt3DInput::QAnalogAxisInput *m_mouseTzYInput;
    Qt3DInput::QButtonAxisInput *m_keyboardTxPosInput;
    Qt3DInput::QButtonAxisInput *m_keyboardTyPosInput;
    Qt3DInput::QButtonAxisInput *m_keyboardTzPosInput;
    Qt3DInput::QButtonAxisInput *m_keyboardTxNegInput;
    Qt3DInput::QButtonAxisInput *m_keyboardTyNegInput;
    Qt3DInput::QButtonAxisInput *m_keyboardTzNegInput;

    Qt3DInput::QKeyboardDevice *m_keyboardDevice;
    Qt3DInput::QMouseDevice *m_mouseDevice;

    Qt3DInput::QLogicalDevice *m_logicalDevice;

    Qt3DLogic::QFrameAction *m_frameAction;

    float m_linearSpeed;
    float m_lookSpeed;

    float m_acceleration;
    float m_deceleration;

    QVector3D m_sceneUp;

    Q_DECLARE_PUBLIC(QAbstractCameraController)
};

} // namespace Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DEXTRAS_QABSTRACTCAMERACONTROLLER_P_H
