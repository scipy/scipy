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

#ifndef QT3D_QUICK_QQMLASPECTENGINE_H
#define QT3D_QUICK_QQMLASPECTENGINE_H

#include <Qt3DCore/QAspectEngine>
#include <Qt3DQuick/qt3dquick_global.h>
#include <QtQml/QQmlEngine>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {
namespace Quick {

class QQmlAspectEnginePrivate;

class Q_3DQUICKSHARED_EXPORT QQmlAspectEngine: public QObject
{
    Q_OBJECT
public:
    enum Status { Null, Ready, Loading, Error };

    explicit QQmlAspectEngine(QObject *parent = nullptr);
    ~QQmlAspectEngine();

    Status status() const;
    void setSource(const QUrl &source);

    QQmlEngine *qmlEngine() const;
    QAspectEngine *aspectEngine() const;

Q_SIGNALS:
    void statusChanged(Status);
    void sceneCreated(QObject* rootObject);

private:
    Q_DECLARE_PRIVATE(QQmlAspectEngine)
    Q_PRIVATE_SLOT(d_func(), void _q_continueExecute())
};

} // namespace Quick
} // namespace Qt3DCore

QT_END_NAMESPACE

#endif // QT3D_QUICK_QQMLASPECTENGINE_H
