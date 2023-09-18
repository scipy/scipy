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

#ifndef QT3DCORE_QABSTRACTASPECT_P_H
#define QT3DCORE_QABSTRACTASPECT_P_H

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

#include <Qt3DCore/qabstractaspect.h>
#include <Qt3DCore/qnodedestroyedchange.h>

#include <Qt3DCore/private/qaspectjobproviderinterface_p.h>
#include <Qt3DCore/private/qbackendnode_p.h>
#include <Qt3DCore/private/qt3dcore_global_p.h>
#include <Qt3DCore/private/qscenechange_p.h>
#include <QtCore/private/qobject_p.h>

#include <QMutex>
#include <QVector>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

class QAbstractAspect;
class QBackendNode;
class QEntity;
class QAspectManager;
class QAbstractAspectJobManager;
class QChangeArbiter;
class QServiceLocator;

namespace Debug {

class Q_3DCORE_PRIVATE_EXPORT AsynchronousCommandReply : public QObject
{
    Q_OBJECT
public:
    explicit AsynchronousCommandReply(const QString &commandName, QObject *parent = nullptr);

    inline QByteArray data() const { return m_data; }
    inline QString commandName() const { return m_commandName; }
    inline bool isFinished() const { return m_finished; }

    void setFinished(bool finished);
    void setData(const QByteArray &data);

Q_SIGNALS:
    void finished(AsynchronousCommandReply *reply);

private:
    QByteArray m_data;
    QString m_commandName;
    bool m_finished;
};

} // Debug

struct NodeTreeChange
{
    enum NodeTreeChangeType {
        Added = 0,
        Removed = 1
    };
    Qt3DCore::QNodeId id;
    const QMetaObject *metaObj;
    NodeTreeChangeType type;
    Qt3DCore::QNode *node;
};

class Q_3DCORE_PRIVATE_EXPORT QAbstractAspectPrivate
        : public QObjectPrivate
        , public QAspectJobProviderInterface
{
public:
    QAbstractAspectPrivate();
    ~QAbstractAspectPrivate();

    void setRootAndCreateNodes(QEntity *rootObject, const QVector<NodeTreeChange> &nodesTreeChanges);

    QServiceLocator *services() const;
    QAbstractAspectJobManager *jobManager() const;

    QVector<QAspectJobPtr> jobsToExecute(qint64 time) override;
    void jobsDone() override;      // called when all the jobs are completed
    void frameDone() override;     // called when frame is completed (after the jobs), safe to wait until next frame here

    QBackendNode *createBackendNode(const NodeTreeChange &change) const;
    void clearBackendNode(const NodeTreeChange &change) const;
    void syncDirtyFrontEndNodes(const QVector<QNode *> &nodes);
    void syncDirtyFrontEndSubNodes(const QVector<NodeRelationshipChange> &nodes);
    virtual void syncDirtyFrontEndNode(QNode *node, QBackendNode *backend, bool firstTime) const;
    void sendPropertyMessages(QNode *node, QBackendNode *backend) const;

    virtual void onEngineAboutToShutdown();

    // TODO: Make these public in 5.8
    template<class Frontend>
    void unregisterBackendType();
    void unregisterBackendType(const QMetaObject &mo);

    Q_DECLARE_PUBLIC(QAbstractAspect)

    enum NodeMapperInfo {
        DefaultMapper = 0,
        SupportsSyncing = 1 << 0
    };
    using BackendNodeMapperAndInfo = QPair<QBackendNodeMapperPtr, NodeMapperInfo>;
    BackendNodeMapperAndInfo mapperForNode(const QMetaObject *metaObj) const;

    QEntity *m_root;
    QNodeId m_rootId;
    QAspectManager *m_aspectManager;
    QAbstractAspectJobManager *m_jobManager;
    QChangeArbiter *m_arbiter;
    QHash<const QMetaObject*, BackendNodeMapperAndInfo> m_backendCreatorFunctors;
    QMutex m_singleShotMutex;
    QVector<QAspectJobPtr> m_singleShotJobs;

    static QAbstractAspectPrivate *get(QAbstractAspect *aspect);
};

template<class Frontend>
void QAbstractAspectPrivate::unregisterBackendType()
{
    unregisterBackendType(Frontend::staticMetaObject);
}

} // Qt3DCore

QT_END_NAMESPACE

#endif // QT3DCORE_QABSTRACTASPECT_P_H
