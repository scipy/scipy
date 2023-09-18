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

#ifndef QT3DCORE_QCHANGEARBITER_P_H
#define QT3DCORE_QCHANGEARBITER_P_H

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

#include <Qt3DCore/qnodeid.h>
#include <Qt3DCore/qscenechange.h>
#include <QtCore/QFlags>
#include <QtCore/QMutex>
#include <QtCore/QObject>
#include <QtCore/QPair>
#include <QtCore/QReadWriteLock>
#include <QtCore/QThreadStorage>
#include <QtCore/QVariant>
#include <QtCore/QVector>

#include <Qt3DCore/private/qlockableobserverinterface_p.h>
#include <Qt3DCore/private/qt3dcore_global_p.h>
#include <Qt3DCore/private/qscenechange_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

class QNode;
class QObservableInterface;
class QAbstractAspectJobManager;
class QSceneObserverInterface;
class QAbstractPostman;
class QScene;


class Q_3DCORE_PRIVATE_EXPORT QAbstractArbiter : public QLockableObserverInterface
{
public:
    virtual QAbstractPostman *postman() const = 0;
    virtual void addDirtyFrontEndNode(QNode *node) = 0;
    virtual void removeDirtyFrontEndNode(QNode *node) = 0;
    QT_WARNING_PUSH
    QT_WARNING_DISABLE_DEPRECATED
    virtual void addDirtyFrontEndNode(QNode *node, QNode *subNode, const char *property, ChangeFlag change) = 0;
    QT_WARNING_POP
};

class Q_3DCORE_PRIVATE_EXPORT QChangeArbiter final
        : public QObject
        , public QAbstractArbiter
{
    Q_OBJECT
public:
    explicit QChangeArbiter(QObject *parent = nullptr);
    ~QChangeArbiter();

    void initialize(Qt3DCore::QAbstractAspectJobManager *jobManager);

    void syncChanges();

    QT_WARNING_PUSH
    QT_WARNING_DISABLE_DEPRECATED
    void registerObserver(QObserverInterface *observer,
                          QNodeId nodeId,
                          ChangeFlags changeFlags = AllChanges);
    QT_WARNING_POP
    void unregisterObserver(QObserverInterface *observer,
                            QNodeId nodeId);

    void sceneChangeEvent(const QSceneChangePtr &e) override;         // QLockableObserverInterface impl
    void sceneChangeEventWithLock(const QSceneChangePtr &e) override; // QLockableObserverInterface impl
    void sceneChangeEventWithLock(const QSceneChangeList &e) override; // QLockableObserverInterface impl

    void addDirtyFrontEndNode(QNode *node) override;
    QT_WARNING_PUSH
    QT_WARNING_DISABLE_DEPRECATED
    void addDirtyFrontEndNode(QNode *node, QNode *subNode, const char *property, ChangeFlag change) override;
    QT_WARNING_POP
    void removeDirtyFrontEndNode(QNode *node) override;
    QVector<QNode *> takeDirtyFrontEndNodes();
    QVector<NodeRelationshipChange> takeDirtyFrontEndSubNodes();

    void setPostman(Qt3DCore::QAbstractPostman *postman);
    void setScene(Qt3DCore::QScene *scene);

    QAbstractPostman *postman() const final;
    QScene *scene() const;

    static void createUnmanagedThreadLocalChangeQueue(void *changeArbiter);
    static void destroyUnmanagedThreadLocalChangeQueue(void *changeArbiter);
    static void createThreadLocalChangeQueue(void *changeArbiter);
    static void destroyThreadLocalChangeQueue(void *changeArbiter);

Q_SIGNALS:
    void receivedChange();
    void syncedChanges();

protected:
    typedef std::vector<QSceneChangePtr> QChangeQueue;
    typedef QPair<ChangeFlags, QObserverInterface *> QObserverPair;
    typedef QVector<QObserverPair> QObserverList;

    void distributeQueueChanges(QChangeQueue *queue);

    QThreadStorage<QChangeQueue *> *tlsChangeQueue();
    void appendChangeQueue(QChangeQueue *queue);
    void removeChangeQueue(QChangeQueue *queue);
    void appendLockingChangeQueue(QChangeQueue *queue);
    void removeLockingChangeQueue(QChangeQueue *queue);

private:
    mutable QRecursiveMutex m_mutex;
    QAbstractAspectJobManager *m_jobManager;

    // The lists of observers indexed by observable (QNodeId).
    // m_nodeObservations is for observables in the main thread's object tree
    QHash<QNodeId, QObserverList> m_nodeObservations;

    // Each thread has a TLS ChangeQueue so we never need to lock whilst
    // receiving a QSceneChange.
    QThreadStorage<QChangeQueue *> m_tlsChangeQueue;

    // We store a list of the ChangeQueue's from each thread. This will only
    // be accessed from the aspect thread during the syncChanges() phase.
    QList<QChangeQueue *> m_changeQueues;
    QList<QChangeQueue *> m_lockingChangeQueues;
    QAbstractPostman *m_postman;
    QScene *m_scene;

    QVector<QNode *> m_dirtyFrontEndNodes;
    QVector<NodeRelationshipChange> m_dirtySubNodeChanges;
};

} // namespace Qt3DCore

QT_END_NAMESPACE

#endif // QT3DCORE_QCHANGEARBITER_P_H
