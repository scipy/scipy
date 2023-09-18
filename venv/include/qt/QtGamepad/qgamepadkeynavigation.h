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

#ifndef QGAMEPADKEYNAVIGATION_H
#define QGAMEPADKEYNAVIGATION_H

#include <QtCore/QObject>
#include <QtCore/QMap>
#include <QtGamepad/qtgamepadglobal.h>

#include <QtGamepad/QGamepadManager>

QT_BEGIN_NAMESPACE

class QKeyEvent;
class QGamepad;
class QGamepadKeyNavigationPrivate;
class Q_GAMEPAD_EXPORT QGamepadKeyNavigation : public QObject
{
    Q_OBJECT
    Q_PROPERTY(bool active READ active WRITE setActive NOTIFY activeChanged)
    Q_PROPERTY(QGamepad *gamepad READ gamepad WRITE setGamepad NOTIFY gamepadChanged)
    Q_PROPERTY(Qt::Key upKey READ upKey WRITE setUpKey NOTIFY upKeyChanged)
    Q_PROPERTY(Qt::Key downKey READ downKey WRITE setDownKey NOTIFY downKeyChanged)
    Q_PROPERTY(Qt::Key leftKey READ leftKey WRITE setLeftKey NOTIFY leftKeyChanged)
    Q_PROPERTY(Qt::Key rightKey READ rightKey WRITE setRightKey NOTIFY rightKeyChanged)
    Q_PROPERTY(Qt::Key buttonAKey READ buttonAKey WRITE setButtonAKey NOTIFY buttonAKeyChanged)
    Q_PROPERTY(Qt::Key buttonBKey READ buttonBKey WRITE setButtonBKey NOTIFY buttonBKeyChanged)
    Q_PROPERTY(Qt::Key buttonXKey READ buttonXKey WRITE setButtonXKey NOTIFY buttonXKeyChanged)
    Q_PROPERTY(Qt::Key buttonYKey READ buttonYKey WRITE setButtonYKey NOTIFY buttonYKeyChanged)
    Q_PROPERTY(Qt::Key buttonSelectKey READ buttonSelectKey WRITE setButtonSelectKey NOTIFY buttonSelectKeyChanged)
    Q_PROPERTY(Qt::Key buttonStartKey READ buttonStartKey WRITE setButtonStartKey NOTIFY buttonStartKeyChanged)
    Q_PROPERTY(Qt::Key buttonGuideKey READ buttonGuideKey WRITE setButtonGuideKey NOTIFY buttonGuideKeyChanged)
    Q_PROPERTY(Qt::Key buttonL1Key READ buttonL1Key WRITE setButtonL1Key NOTIFY buttonL1KeyChanged)
    Q_PROPERTY(Qt::Key buttonR1Key READ buttonR1Key WRITE setButtonR1Key NOTIFY buttonR1KeyChanged)
    Q_PROPERTY(Qt::Key buttonL2Key READ buttonL2Key WRITE setButtonL2Key NOTIFY buttonL2KeyChanged)
    Q_PROPERTY(Qt::Key buttonR2Key READ buttonR2Key WRITE setButtonR2Key NOTIFY buttonR2KeyChanged)
    Q_PROPERTY(Qt::Key buttonL3Key READ buttonL3Key WRITE setButtonL3Key NOTIFY buttonL3KeyChanged)
    Q_PROPERTY(Qt::Key buttonR3Key READ buttonR3Key WRITE setButtonR3Key NOTIFY buttonR3KeyChanged)
public:
    explicit QGamepadKeyNavigation(QObject *parent = nullptr);

    bool active() const;
    QGamepad *gamepad() const;

    Qt::Key upKey() const;
    Qt::Key downKey() const;
    Qt::Key leftKey() const;
    Qt::Key rightKey() const;
    Qt::Key buttonAKey() const;
    Qt::Key buttonBKey() const;
    Qt::Key buttonXKey() const;
    Qt::Key buttonYKey() const;
    Qt::Key buttonSelectKey() const;
    Qt::Key buttonStartKey() const;
    Qt::Key buttonGuideKey() const;
    Qt::Key buttonL1Key() const;
    Qt::Key buttonR1Key() const;
    Qt::Key buttonL2Key() const;
    Qt::Key buttonR2Key() const;
    Qt::Key buttonL3Key() const;
    Qt::Key buttonR3Key() const;

Q_SIGNALS:
    void activeChanged(bool isActive);
    void gamepadChanged(QGamepad *gamepad);

    void upKeyChanged(Qt::Key key);
    void downKeyChanged(Qt::Key key);
    void leftKeyChanged(Qt::Key key);
    void rightKeyChanged(Qt::Key key);
    void buttonAKeyChanged(Qt::Key key);
    void buttonBKeyChanged(Qt::Key key);
    void buttonXKeyChanged(Qt::Key key);
    void buttonYKeyChanged(Qt::Key key);
    void buttonSelectKeyChanged(Qt::Key key);
    void buttonStartKeyChanged(Qt::Key key);
    void buttonGuideKeyChanged(Qt::Key key);
    void buttonL1KeyChanged(Qt::Key key);
    void buttonR1KeyChanged(Qt::Key key);
    void buttonL2KeyChanged(Qt::Key key);
    void buttonR2KeyChanged(Qt::Key key);
    void buttonL3KeyChanged(Qt::Key key);
    void buttonR3KeyChanged(Qt::Key key);

public Q_SLOTS:
    void setActive(bool isActive);
    void setGamepad(QGamepad *gamepad);

    void setUpKey(Qt::Key key);
    void setDownKey(Qt::Key key);
    void setLeftKey(Qt::Key key);
    void setRightKey(Qt::Key key);
    void setButtonAKey(Qt::Key key);
    void setButtonBKey(Qt::Key key);
    void setButtonXKey(Qt::Key key);
    void setButtonYKey(Qt::Key key);
    void setButtonSelectKey(Qt::Key key);
    void setButtonStartKey(Qt::Key key);
    void setButtonGuideKey(Qt::Key key);
    void setButtonL1Key(Qt::Key key);
    void setButtonR1Key(Qt::Key key);
    void setButtonL2Key(Qt::Key key);
    void setButtonR2Key(Qt::Key key);
    void setButtonL3Key(Qt::Key key);
    void setButtonR3Key(Qt::Key key);

private:
    Q_DECLARE_PRIVATE(QGamepadKeyNavigation)
    Q_DISABLE_COPY(QGamepadKeyNavigation)
    Q_PRIVATE_SLOT(d_func(), void _q_processGamepadButtonPressEvent(int, QGamepadManager::GamepadButton, double))
    Q_PRIVATE_SLOT(d_func(), void _q_processGamepadButtonReleaseEvent(int, QGamepadManager::GamepadButton))
};

QT_END_NAMESPACE

#endif // QGAMEPADKEYNAVIGATION_H
