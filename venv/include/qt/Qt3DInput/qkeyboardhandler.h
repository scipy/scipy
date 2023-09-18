/****************************************************************************
**
** Copyright (C) 2014 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DINPUT_INPUT_QKEYBOARDHANDLER_H
#define QT3DINPUT_INPUT_QKEYBOARDHANDLER_H

#include <Qt3DInput/qt3dinput_global.h>
#include <Qt3DInput/qkeyevent.h>
#include <Qt3DCore/qcomponent.h>

QT_BEGIN_NAMESPACE

namespace Qt3DInput {

class QKeyboardDevice;
class QKeyboardHandlerPrivate;

class Q_3DINPUTSHARED_EXPORT QKeyboardHandler : public Qt3DCore::QComponent
{
    Q_OBJECT
    Q_PROPERTY(Qt3DInput::QKeyboardDevice *sourceDevice READ sourceDevice WRITE setSourceDevice NOTIFY sourceDeviceChanged)
    Q_PROPERTY(bool focus READ focus WRITE setFocus NOTIFY focusChanged)
public:
    explicit QKeyboardHandler(QNode *parent = nullptr);
    ~QKeyboardHandler();

    Qt3DInput::QKeyboardDevice *sourceDevice() const;
    bool focus() const;

public Q_SLOTS:
    void setSourceDevice(Qt3DInput::QKeyboardDevice *keyboardDevice);
    void setFocus(bool focus);

Q_SIGNALS:
    void sourceDeviceChanged(QKeyboardDevice *keyboardDevice);
    void focusChanged(bool focus);

    void digit0Pressed(Qt3DInput::QKeyEvent *event);
    void digit1Pressed(Qt3DInput::QKeyEvent *event);
    void digit2Pressed(Qt3DInput::QKeyEvent *event);
    void digit3Pressed(Qt3DInput::QKeyEvent *event);
    void digit4Pressed(Qt3DInput::QKeyEvent *event);
    void digit5Pressed(Qt3DInput::QKeyEvent *event);
    void digit6Pressed(Qt3DInput::QKeyEvent *event);
    void digit7Pressed(Qt3DInput::QKeyEvent *event);
    void digit8Pressed(Qt3DInput::QKeyEvent *event);
    void digit9Pressed(Qt3DInput::QKeyEvent *event);

    void leftPressed(Qt3DInput::QKeyEvent *event);
    void rightPressed(Qt3DInput::QKeyEvent *event);
    void upPressed(Qt3DInput::QKeyEvent *event);
    void downPressed(Qt3DInput::QKeyEvent *event);
    void tabPressed(Qt3DInput::QKeyEvent *event);
    void backtabPressed(Qt3DInput::QKeyEvent *event);

    void asteriskPressed(Qt3DInput::QKeyEvent *event);
    void numberSignPressed(Qt3DInput::QKeyEvent *event);
    void escapePressed(Qt3DInput::QKeyEvent *event);
    void returnPressed(Qt3DInput::QKeyEvent *event);
    void enterPressed(Qt3DInput::QKeyEvent *event);
    void deletePressed(Qt3DInput::QKeyEvent *event);
    void spacePressed(Qt3DInput::QKeyEvent *event);
    void backPressed(Qt3DInput::QKeyEvent *event);
    void cancelPressed(Qt3DInput::QKeyEvent *event);
    void selectPressed(Qt3DInput::QKeyEvent *event);
    void yesPressed(Qt3DInput::QKeyEvent *event);
    void noPressed(Qt3DInput::QKeyEvent *event);
    void context1Pressed(Qt3DInput::QKeyEvent *event);
    void context2Pressed(Qt3DInput::QKeyEvent *event);
    void context3Pressed(Qt3DInput::QKeyEvent *event);
    void context4Pressed(Qt3DInput::QKeyEvent *event);
    void callPressed(Qt3DInput::QKeyEvent *event);
    void hangupPressed(Qt3DInput::QKeyEvent *event);
    void flipPressed(Qt3DInput::QKeyEvent *event);
    void menuPressed(Qt3DInput::QKeyEvent *event);
    void volumeUpPressed(Qt3DInput::QKeyEvent *event);
    void volumeDownPressed(Qt3DInput::QKeyEvent *event);

    void pressed(Qt3DInput::QKeyEvent *event);
    void released(Qt3DInput::QKeyEvent *event);

protected:
    // TODO Unused remove in Qt6
    void sceneChangeEvent(const Qt3DCore::QSceneChangePtr &change) override;

private:
    Q_DECLARE_PRIVATE(QKeyboardHandler)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

} // namespace Qt3DInput

QT_END_NAMESPACE

#endif // QT3DINPUT_INPUT_QKEYBOARDHANDLER_H
