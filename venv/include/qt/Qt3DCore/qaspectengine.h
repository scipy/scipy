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

#ifndef QT3DCORE_QASPECTENGINE_H
#define QT3DCORE_QASPECTENGINE_H

#include <Qt3DCore/qt3dcore_global.h>
#include <QtCore/QObject>
#include <QtCore/QVector>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

class QAbstractAspect;
class QAspectThread;
class QAspectEnginePrivate;
class QEntity;
class QNode;

typedef QSharedPointer<QEntity> QEntityPtr;

class Q_3DCORESHARED_EXPORT QAspectEngine : public QObject
{
    Q_OBJECT
public:
    enum RunMode {
        Manual = 0,
        Automatic
    };
    Q_ENUM(RunMode)

    explicit QAspectEngine(QObject *parent = nullptr);
    ~QAspectEngine();

    void setRootEntity(QEntityPtr root);
    QEntityPtr rootEntity() const;

    void setRunMode(RunMode mode);
    RunMode runMode() const;

    void registerAspect(QAbstractAspect *aspect);
    void registerAspect(const QString &name);
    void unregisterAspect(QAbstractAspect *aspect);
    void unregisterAspect(const QString &name);

    QVector<QAbstractAspect*> aspects() const;

    QVariant executeCommand(const QString &command);

    void processFrame();

private:
    Q_DECLARE_PRIVATE(QAspectEngine)
};

} // namespace Qt3DCore

QT_END_NAMESPACE


#endif // QT3DCORE_QASPECTENGINE_H
