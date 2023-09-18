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

#ifndef QT3DINPUT_INPUT_QKEYBOARDDEVICE_H
#define QT3DINPUT_INPUT_QKEYBOARDDEVICE_H

#include <Qt3DInput/QKeyboardHandler>
#include <Qt3DInput/qt3dinput_global.h>
#include <Qt3DInput/qabstractphysicaldevice.h>

QT_BEGIN_NAMESPACE

namespace Qt3DInput {

class QKeyboardDevicePrivate;
class QKeyboardHandler;

class Q_3DINPUTSHARED_EXPORT QKeyboardDevice : public Qt3DInput::QAbstractPhysicalDevice
{
    Q_OBJECT
    Q_PROPERTY(Qt3DInput::QKeyboardHandler *activeInput READ activeInput NOTIFY activeInputChanged)

public:
    explicit QKeyboardDevice(QNode *parent = nullptr);
    ~QKeyboardDevice();

    QKeyboardHandler *activeInput() const;

    int axisCount() const final;
    int buttonCount() const final;
    QStringList axisNames() const final;
    QStringList buttonNames() const final;
    int axisIdentifier(const QString &name) const final;
    int buttonIdentifier(const QString &name) const final;

protected:
    explicit QKeyboardDevice(QKeyboardDevicePrivate &dd, QNode *parent = nullptr);
    // TODO Unused remove in Qt6
    void sceneChangeEvent(const Qt3DCore::QSceneChangePtr &change) override;

Q_SIGNALS:
    void activeInputChanged(QKeyboardHandler *activeInput);

private:
    Q_DECLARE_PRIVATE(QKeyboardDevice)
    void setActiveInput(QKeyboardHandler *activeInput);
};

} // namespace Qt3DInput

QT_END_NAMESPACE

#endif // QT3DINPUT_INPUT_QKEYBOARDDEVICE_H
