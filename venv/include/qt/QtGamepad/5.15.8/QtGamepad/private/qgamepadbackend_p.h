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

#ifndef QGAMEPADBACKEND_P_H
#define QGAMEPADBACKEND_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <QtGamepad/qtgamepadglobal.h>

#include <QtCore/QObject>
#include <QtGamepad/QGamepadManager>

QT_BEGIN_NAMESPACE

class Q_GAMEPAD_EXPORT QGamepadBackend : public QObject
{
    Q_OBJECT
public:
    template <typename T>
    struct AxisInfo {
        AxisInfo(T minValue = 0, T maxValue = 1, QGamepadManager::GamepadAxis gamepadAxis = QGamepadManager::AxisInvalid)
            : minValue(minValue)
            , maxValue(maxValue)
            , gamepadAxis(gamepadAxis)
        {}

        virtual double normalized(T value) const
        {
            if (gamepadAxis == QGamepadManager::AxisInvalid && minValue >= 0) {
                //case 0.0 - 1.0
                return double(value - minValue) / double(maxValue - minValue);
            } else {
                //case -1.0 - 1.0
                return 2.0 * double(value - minValue) / double(maxValue - minValue) - 1.0;
            }
        }
        T minValue;
        T maxValue;
        QGamepadManager::GamepadAxis gamepadAxis;
    };

public:
    explicit QGamepadBackend(QObject *parent = nullptr);
    virtual bool isConfigurationNeeded(int deviceId);
    virtual void resetConfiguration(int deviceId);
    virtual bool configureButton(int deviceId, QGamepadManager::GamepadButton button);
    virtual bool configureAxis(int deviceId, QGamepadManager::GamepadAxis axis);
    virtual bool setCancelConfigureButton(int deviceId, QGamepadManager::GamepadButton button);
    virtual void setSettingsFile(const QString &file);
    virtual void saveSettings(int productId, const QVariant &value);
    virtual QVariant readSettings(int productId);

public Q_SLOTS:
    virtual bool start();
    virtual void stop();

Q_SIGNALS:
    void buttonConfigured(int deviceId, QGamepadManager::GamepadButton button);
    void axisConfigured(int deviceId, QGamepadManager::GamepadAxis axis);
    void configurationCanceled(int deviceId);
    void gamepadAdded(int deviceId);
    void gamepadNamed(int deviceId, const QString &name);
    void gamepadRemoved(int deviceId);
    void gamepadAxisMoved(int deviceId, QGamepadManager::GamepadAxis axis, double value);
    void gamepadButtonPressed(int deviceId, QGamepadManager::GamepadButton button, double value);
    void gamepadButtonReleased(int deviceId, QGamepadManager::GamepadButton button);

protected:
    QString m_settingsFilePath;
};

QT_END_NAMESPACE

#endif // QGAMEPADBACKEND_P_H
