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

#ifndef QT3DCORE_QASPECTENGINE_P_H
#define QT3DCORE_QASPECTENGINE_P_H

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

#include <Qt3DCore/private/qt3dcore_global_p.h>
#include <Qt3DCore/qnodecreatedchange.h>
#include <QtCore/qsharedpointer.h>

#include <Qt3DCore/private/qaspectfactory_p.h>
#include <Qt3DCore/private/qaspectengine_p.h>
#include <QtCore/private/qobject_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

class QEntity;
class QNode;
class QAspectManager;
class QPostman;
class QScene;

class Q_3DCORE_PRIVATE_EXPORT QAspectEnginePrivate : public QObjectPrivate
{
public:
    QAspectEnginePrivate();
    ~QAspectEnginePrivate();

    Q_DECLARE_PUBLIC(QAspectEngine)

    QAspectFactory m_factory;
    QAspectManager *m_aspectManager;
    QPostman *m_postman;
    QScene *m_scene;
    QSharedPointer<QEntity> m_root;
    QVector<QAbstractAspect*> m_aspects;
    QHash<QString, QAbstractAspect *> m_namedAspects;
    bool m_initialized;
    QAspectEngine::RunMode m_runMode;

    void initialize();
    void shutdown();

    void exitSimulationLoop();

    void initNodeTree(QNode *node);
    void initNode(QNode *node);
    void initEntity(QEntity *entity);
    void addNode(QNode *node);
    void removeNode(QNode *node);

    static QAspectEnginePrivate *get(QAspectEngine *engine);
};

} // Qt3D

QT_END_NAMESPACE

#endif // QT3DCORE_QASPECTENGINE_P_H
