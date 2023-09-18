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

#ifndef QT3DINPUT_INPUT_QUICK_QUICK3DLOGICALDEVICE_H
#define QT3DINPUT_INPUT_QUICK_QUICK3DLOGICALDEVICE_H

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

#include <Qt3DInput/QAxis>
#include <Qt3DInput/QAction>
#include <Qt3DInput/QLogicalDevice>
#include <QtQml/QQmlListProperty>

#include <Qt3DQuickInput/private/qt3dquickinput_global_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DInput {
namespace Input {
namespace Quick {

class Q_3DQUICKINPUTSHARED_PRIVATE_EXPORT Quick3DLogicalDevice  : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QQmlListProperty<Qt3DInput::QAxis> axes READ qmlAxes CONSTANT)
    Q_PROPERTY(QQmlListProperty<Qt3DInput::QAction> actions READ qmlActions CONSTANT)
public:
    explicit Quick3DLogicalDevice(QObject *parent = nullptr);

    inline QLogicalDevice *parentLogicalDevice() const { return qobject_cast<QLogicalDevice *>(parent()); }
    QQmlListProperty<QAxis> qmlAxes();
    QQmlListProperty<QAction> qmlActions();

private:
    static void appendAxis(QQmlListProperty<QAxis> *list, QAxis *axes);
    static QAxis *axisAt(QQmlListProperty<QAxis> *list, int index);
    static int axesCount(QQmlListProperty<QAxis> *list);
    static void clearAxes(QQmlListProperty<QAxis> *list);

    static void appendAction(QQmlListProperty<QAction> *list, QAction *action);
    static QAction *actionAt(QQmlListProperty<QAction> *list, int index);
    static int actionCount(QQmlListProperty<QAction> *list);
    static void clearActions(QQmlListProperty<QAction> *list);
};

} // namespace Quick
} // namespace Input
} // namespace Qt3DInput

QT_END_NAMESPACE

#endif // QT3DINPUT_INPUT_QUICK_QUICK3DLOGICALDEVICE_H
