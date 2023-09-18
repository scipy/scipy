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

#ifndef QT3DINPUT_INPUT_MOUSEDEVICE_H
#define QT3DINPUT_INPUT_MOUSEDEVICE_H

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

#include <Qt3DInput/QMouseEvent>

#include <Qt3DInput/private/qabstractphysicaldevicebackendnode_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DInput {

class QInputAspect;

namespace Input {

class InputHandler;

class Q_AUTOTEST_EXPORT MouseDevice : public Qt3DInput::QAbstractPhysicalDeviceBackendNode
{
public:
    struct MouseState {

        MouseState()
            : xAxis(0.0f)
            , yAxis(0.0f)
            , wXAxis(0.0f)
            , wYAxis(0.0f)
            , leftPressed(false)
            , rightPressed(false)
            , centerPressed(false)
        {}

        float xAxis;
        float yAxis;
        float wXAxis;
        float wYAxis;
        bool leftPressed;
        bool rightPressed;
        bool centerPressed;
    };

    MouseDevice();
    ~MouseDevice();

    void setInputHandler(InputHandler *handler);
    InputHandler *inputHandler() const;

    float axisValue(int axisIdentifier) const override;
    bool isButtonPressed(int buttonIdentifier) const override;

    void updateMouseEvents(const QList<QT_PREPEND_NAMESPACE(QMouseEvent)> &events);
#if QT_CONFIG(wheelevent)
    void updateWheelEvents(const QList<QT_PREPEND_NAMESPACE(QWheelEvent)> &events);
#endif

    MouseState mouseState() const;
    QPointF previousPos() const;
    bool wasPressed() const;
    float sensitivity() const;
    bool updateAxesContinuously() const;

    void syncFromFrontEnd(const Qt3DCore::QNode *frontEnd, bool firstTime) override;

private:
    InputHandler *m_inputHandler;

    MouseState m_mouseState;
    QPointF m_previousPos;
    bool m_wasPressed;
    float m_sensitivity;
    bool m_updateAxesContinuously;
};

class MouseDeviceFunctor : public Qt3DCore::QBackendNodeMapper
{
public:
    explicit MouseDeviceFunctor(Qt3DInput::QInputAspect *inputAspect, InputHandler *handler);

    Qt3DCore::QBackendNode *create(const Qt3DCore::QNodeCreatedChangeBasePtr &change) const override;
    Qt3DCore::QBackendNode *get(Qt3DCore::QNodeId id) const override;
    void destroy(Qt3DCore::QNodeId id) const override;

private:
    QInputAspect *m_inputAspect;
    InputHandler *m_handler;
};

} // namespace Input
} // namespace Qt3DInput

QT_END_NAMESPACE

#endif // MOUSEDEVICE_H
