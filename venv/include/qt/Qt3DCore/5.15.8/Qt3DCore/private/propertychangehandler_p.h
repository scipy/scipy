/****************************************************************************
**
** Copyright (C) 2014 Klarälvdalens Datakonsult AB, a KDAB Group company, info@kdab.com, author Milian Wolff <milian.wolff@kdab.com>
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtWebChannel module of the Qt Toolkit.
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

#ifndef SIGNALHANDLER_H
#define SIGNALHANDLER_H

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

#include <Qt3DCore/qt3dcore_global.h>
#include <QtCore/QDebug>
#include <QtCore/QHash>
#include <QtCore/QMetaMethod>
#include <QtCore/QObject>
#include <QtCore/QVector>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

class Q_3DCORESHARED_EXPORT PropertyChangeHandlerBase : public QObject
{
    Q_OBJECT
public:
    PropertyChangeHandlerBase(QObject *parent = nullptr);

    /**
     * Connect to the change signal of @p property in @p object.
     */
    void connectToPropertyChange(const QObject *object, int propertyIndex);

    /**
     * Disconnect from the change signal of @p property in @p object.
     */
    void disconnectFromPropertyChange(const QObject *object, int propertyIndex);
};

/**
 * The property change handler is similar to QSignalSpy, but geared towards the usecase of Qt3D.
 *
 * It allows connecting to any number of property change signals of the receiver object and forwards
 * the signal invocations to the Receiver by calling its propertyChanged function.
 */
template<class Receiver>
class PropertyChangeHandler : public PropertyChangeHandlerBase
{
public:
    PropertyChangeHandler(Receiver *receiver, QObject *parent = nullptr);

    /**
     * @internal
     *
     * Custom implementation of qt_metacall which calls propertyChanged() on the receiver for
     * connected signals.
     */
    int qt_metacall(QMetaObject::Call call, int methodId, void **args) override;

private:
    Receiver *m_receiver;
};

template<class Receiver>
PropertyChangeHandler<Receiver>::PropertyChangeHandler(Receiver *receiver, QObject *parent)
    : PropertyChangeHandlerBase(parent)
    , m_receiver(receiver)
{
}

template<class Receiver>
int PropertyChangeHandler<Receiver>::qt_metacall(QMetaObject::Call call, int methodId, void **args)
{
    methodId = QObject::qt_metacall(call, methodId, args);
    if (methodId < 0)
        return methodId;

    if (call == QMetaObject::InvokeMetaMethod) {
        m_receiver->propertyChanged(methodId);
        return -1;
    }
    return methodId;
}

}

QT_END_NAMESPACE

#endif // SIGNALHANDLER_H
