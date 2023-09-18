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

#ifndef QT3DINPUT_QLOGICALDEVICE_H
#define QT3DINPUT_QLOGICALDEVICE_H

#include <Qt3DInput/qt3dinput_global.h>
#include <Qt3DCore/qcomponent.h>

QT_BEGIN_NAMESPACE

namespace Qt3DInput {

class QLogicalDevicePrivate;
class QAction;
class QAxis;

class Q_3DINPUTSHARED_EXPORT QLogicalDevice : public Qt3DCore::QComponent
{
    Q_OBJECT
public:
    explicit QLogicalDevice(Qt3DCore::QNode *parent = nullptr);
    ~QLogicalDevice();

    void addAction(QAction *action);
    void removeAction(QAction *action);
    QVector<QAction *> actions() const;

    void addAxis(QAxis *axis);
    void removeAxis(QAxis *axis);
    QVector<QAxis *> axes() const;

private:
    Q_DECLARE_PRIVATE(QLogicalDevice)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

} // Qt3DInput

QT_END_NAMESPACE

#endif // QT3DINPUT_QLOGICALDEVICE_H
