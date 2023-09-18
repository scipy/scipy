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

#ifndef QGAMEPAD_H
#define QGAMEPAD_H

#include <QtCore/QObject>
#include <QtGamepad/qtgamepadglobal.h>
#include <QtGamepad/QGamepadManager>

QT_BEGIN_NAMESPACE

class QGamepadPrivate;

class Q_GAMEPAD_EXPORT QGamepad : public QObject
{
    Q_OBJECT
    Q_PROPERTY(int deviceId READ deviceId WRITE setDeviceId NOTIFY deviceIdChanged)
    Q_PROPERTY(bool connected READ isConnected NOTIFY connectedChanged)
    Q_PROPERTY(QString name READ name NOTIFY nameChanged)
    Q_PROPERTY(double axisLeftX READ axisLeftX NOTIFY axisLeftXChanged)
    Q_PROPERTY(double axisLeftY READ axisLeftY NOTIFY axisLeftYChanged)
    Q_PROPERTY(double axisRightX READ axisRightX NOTIFY axisRightXChanged)
    Q_PROPERTY(double axisRightY READ axisRightY NOTIFY axisRightYChanged)
    Q_PROPERTY(bool buttonA READ buttonA NOTIFY buttonAChanged)
    Q_PROPERTY(bool buttonB READ buttonB NOTIFY buttonBChanged)
    Q_PROPERTY(bool buttonX READ buttonX NOTIFY buttonXChanged)
    Q_PROPERTY(bool buttonY READ buttonY NOTIFY buttonYChanged)
    Q_PROPERTY(bool buttonL1 READ buttonL1 NOTIFY buttonL1Changed)
    Q_PROPERTY(bool buttonR1 READ buttonR1 NOTIFY buttonR1Changed)
    Q_PROPERTY(double buttonL2 READ buttonL2 NOTIFY buttonL2Changed)
    Q_PROPERTY(double buttonR2 READ buttonR2 NOTIFY buttonR2Changed)
    Q_PROPERTY(bool buttonSelect READ buttonSelect NOTIFY buttonSelectChanged)
    Q_PROPERTY(bool buttonStart READ buttonStart NOTIFY buttonStartChanged)
    Q_PROPERTY(bool buttonL3 READ buttonL3 NOTIFY buttonL3Changed)
    Q_PROPERTY(bool buttonR3 READ buttonR3 NOTIFY buttonR3Changed)
    Q_PROPERTY(bool buttonUp READ buttonUp NOTIFY buttonUpChanged)
    Q_PROPERTY(bool buttonDown READ buttonDown NOTIFY buttonDownChanged)
    Q_PROPERTY(bool buttonLeft READ buttonLeft NOTIFY buttonLeftChanged)
    Q_PROPERTY(bool buttonRight READ buttonRight NOTIFY buttonRightChanged)
    Q_PROPERTY(bool buttonCenter READ buttonCenter NOTIFY buttonCenterChanged)
    Q_PROPERTY(bool buttonGuide READ buttonGuide NOTIFY buttonGuideChanged)
public:
    explicit QGamepad(int deviceId = 0, QObject *parent = nullptr);
    ~QGamepad();

    int deviceId() const;

    bool isConnected() const;

    QString name() const;

    double axisLeftX() const;
    double axisLeftY() const;
    double axisRightX() const;
    double axisRightY() const;
    bool buttonA() const;
    bool buttonB() const;
    bool buttonX() const;
    bool buttonY() const;
    bool buttonL1() const;
    bool buttonR1() const;
    double buttonL2() const;
    double buttonR2() const;
    bool buttonSelect() const;
    bool buttonStart() const;
    bool buttonL3() const;
    bool buttonR3() const;
    bool buttonUp() const;
    bool buttonDown() const;
    bool buttonLeft() const;
    bool buttonRight() const;
    bool buttonCenter() const;
    bool buttonGuide() const;

Q_SIGNALS:

    void deviceIdChanged(int value);
    void connectedChanged(bool value);
    void nameChanged(QString value);
    void axisLeftXChanged(double value);
    void axisLeftYChanged(double value);
    void axisRightXChanged(double value);
    void axisRightYChanged(double value);
    void buttonAChanged(bool value);
    void buttonBChanged(bool value);
    void buttonXChanged(bool value);
    void buttonYChanged(bool value);
    void buttonL1Changed(bool value);
    void buttonR1Changed(bool value);
    void buttonL2Changed(double value);
    void buttonR2Changed(double value);
    void buttonSelectChanged(bool value);
    void buttonStartChanged(bool value);
    void buttonL3Changed(bool value);
    void buttonR3Changed(bool value);
    void buttonUpChanged(bool value);
    void buttonDownChanged(bool value);
    void buttonLeftChanged(bool value);
    void buttonRightChanged(bool value);
    void buttonCenterChanged(bool value);
    void buttonGuideChanged(bool value);

public Q_SLOTS:

    void setDeviceId(int number);

private:
    Q_DECLARE_PRIVATE(QGamepad)
    Q_DISABLE_COPY(QGamepad)
    Q_PRIVATE_SLOT(d_func(), void _q_handleGamepadConnected(int))
    Q_PRIVATE_SLOT(d_func(), void _q_handleGamepadNameChanged(int, const QString &))
    Q_PRIVATE_SLOT(d_func(), void _q_handleGamepadDisconnected(int))
    Q_PRIVATE_SLOT(d_func(), void _q_handleGamepadAxisEvent(int, QGamepadManager::GamepadAxis, double))
    Q_PRIVATE_SLOT(d_func(), void _q_handleGamepadButtonPressEvent(int, QGamepadManager::GamepadButton, double))
    Q_PRIVATE_SLOT(d_func(), void _q_handleGamepadButtonReleaseEvent(int, QGamepadManager::GamepadButton))
};

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QGamepad*)

#endif // QGAMEPAD_H
