/****************************************************************************
**
** Copyright (C) 2015 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the Qt Gamepad module
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

#ifndef QGAMEPADMANAGER_H
#define QGAMEPADMANAGER_H

#include <QtCore/QObject>
#include <QtCore/QMap>
#include <QtGamepad/qtgamepadglobal.h>

QT_BEGIN_NAMESPACE

class QGamepadBackend;
class QGamepad;
class QGamepadManagerPrivate;

class Q_GAMEPAD_EXPORT QGamepadManager : public QObject
{
    Q_OBJECT
    Q_FLAGS(GamepadButton GamepadButtons)
    Q_FLAGS(GamepadAxis GamepadAxes)
    Q_PROPERTY(QList<int> connectedGamepads READ connectedGamepads NOTIFY connectedGamepadsChanged)

public:
    enum GamepadButton {
        ButtonInvalid = -1,
        ButtonA = 0,
        ButtonB,
        ButtonX,
        ButtonY,
        ButtonL1,
        ButtonR1,
        ButtonL2,
        ButtonR2,
        ButtonSelect,
        ButtonStart,
        ButtonL3,
        ButtonR3,
        ButtonUp,
        ButtonDown,
        ButtonRight,
        ButtonLeft,
        ButtonCenter,
        ButtonGuide
    };
    Q_DECLARE_FLAGS(GamepadButtons, GamepadButton)

    enum GamepadAxis {
        AxisInvalid = -1,
        AxisLeftX = 0,
        AxisLeftY,
        AxisRightX,
        AxisRightY
    };
    Q_DECLARE_FLAGS(GamepadAxes, GamepadAxis)

    static QGamepadManager* instance();

    bool isGamepadConnected(int deviceId) const;
    QString gamepadName(int deviceId) const;
    const QList<int> connectedGamepads() const;

public Q_SLOTS:
    bool isConfigurationNeeded(int deviceId) const;
    bool configureButton(int deviceId, GamepadButton button);
    bool configureAxis(int deviceId, GamepadAxis axis);
    bool setCancelConfigureButton(int deviceId, GamepadButton button);
    void resetConfiguration(int deviceId);
    void setSettingsFile(const QString &file);

Q_SIGNALS:
    void connectedGamepadsChanged();
    void gamepadConnected(int deviceId);
    void gamepadNameChanged(int deviceId, const QString &name);
    void gamepadDisconnected(int deviceId);
    void gamepadAxisEvent(int deviceId, QGamepadManager::GamepadAxis axis, double value);
    void gamepadButtonPressEvent(int deviceId, QGamepadManager::GamepadButton button, double value);
    void gamepadButtonReleaseEvent(int deviceId, QGamepadManager::GamepadButton button);
    void buttonConfigured(int deviceId, QGamepadManager::GamepadButton button);
    void axisConfigured(int deviceId, QGamepadManager::GamepadAxis axis);
    void configurationCanceled(int deviceId);

private:
    QGamepadManager();
    ~QGamepadManager();

    Q_DECLARE_PRIVATE(QGamepadManager)
    Q_DISABLE_COPY(QGamepadManager)
    Q_PRIVATE_SLOT(d_func(), void _q_forwardGamepadConnected(int))
    Q_PRIVATE_SLOT(d_func(), void _q_forwardGamepadNameChanged(int, const QString&))
    Q_PRIVATE_SLOT(d_func(), void _q_forwardGamepadDisconnected(int))
    Q_PRIVATE_SLOT(d_func(), void _q_forwardGamepadAxisEvent(int, QGamepadManager::GamepadAxis, double))
    Q_PRIVATE_SLOT(d_func(), void _q_forwardGamepadButtonPressEvent(int, QGamepadManager::GamepadButton, double))
    Q_PRIVATE_SLOT(d_func(), void _q_forwardGamepadButtonReleaseEvent(int, QGamepadManager::GamepadButton))
};

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QGamepadManager::GamepadButton)
Q_DECLARE_METATYPE(QGamepadManager::GamepadAxis)

#endif // QGAMEPADMANAGER_H
