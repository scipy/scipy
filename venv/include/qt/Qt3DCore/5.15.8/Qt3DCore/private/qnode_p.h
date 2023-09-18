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

#ifndef QT3DCORE_QNODE_P_H
#define QT3DCORE_QNODE_P_H

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

#include <Qt3DCore/qnode.h>

#include <functional>
#include <vector>

#include <Qt3DCore/private/propertychangehandler_p.h>
#include <Qt3DCore/private/qchangearbiter_p.h>
#include <Qt3DCore/private/qobservableinterface_p.h>
#include <Qt3DCore/private/qt3dcore_global_p.h>
#include <QtCore/private/qobject_p.h>
#include <QQueue>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

class QNode;
class QAspectEngine;

class Q_3DCORE_PRIVATE_EXPORT QNodePrivate : public QObjectPrivate, public QObservableInterface
{
public:
    QNodePrivate();
    ~QNodePrivate();

    void init(QNode *parent);

    virtual void setScene(QScene *scene);
    QScene *scene() const;

    void setArbiter(QLockableObserverInterface *arbiter) override;

    void notifyPropertyChange(const char *name, const QVariant &value);
    void notifyDynamicPropertyChange(const QByteArray &name, const QVariant &value);
    void notifyObservers(const QSceneChangePtr &change) override;

    void insertTree(QNode *treeRoot, int depth = 0);
    void updatePropertyTrackMode();

    void update();
    QT_WARNING_PUSH
    QT_WARNING_DISABLE_DEPRECATED
    void updateNode(QNode *node, const char* property, ChangeFlag change);
    QT_WARNING_POP

    Q_DECLARE_PUBLIC(QNode)

    // For now this just protects access to the m_changeArbiter.
    // Later on we may decide to extend support for multiple observers.
    QAbstractArbiter *m_changeArbiter;
    QMetaObject *m_typeInfo;
    QScene *m_scene;
    mutable QNodeId m_id;
    QNodeId m_parentId; // Store this so we have it even in parent's QObject dtor
    bool m_blockNotifications;
    bool m_hasBackendNode;
    bool m_enabled;
    bool m_notifiedParent;
    QNode::PropertyTrackingMode m_defaultPropertyTrackMode;
    QHash<QString, QNode::PropertyTrackingMode> m_trackedPropertiesOverrides;

    static QNodePrivate *get(QNode *q);
    static const QNodePrivate *get(const QNode *q);
    static void nodePtrDeleter(QNode *q);

    template<typename Caller, typename NodeType>
    using DestructionFunctionPointer = void (Caller::*)(NodeType *);

    template<typename Caller, typename NodeType, typename PropertyType>
    void registerDestructionHelper(NodeType *, DestructionFunctionPointer<Caller, NodeType>, PropertyType);

    template<typename Caller, typename NodeType>
    void registerDestructionHelper(NodeType *node, DestructionFunctionPointer<Caller, NodeType> func, NodeType *&)
    {
        // If the node is destoyed, we make sure not to keep a dangling pointer to it
        Q_Q(QNode);
        auto f = [q, func]() { (static_cast<Caller *>(q)->*func)(nullptr); };
        m_destructionConnections.push_back({node, QObject::connect(node, &QNode::nodeDestroyed, f)});
    }

    template<typename Caller, typename NodeType>
    void registerDestructionHelper(NodeType *node, DestructionFunctionPointer<Caller, NodeType> func, QVector<NodeType*> &)
    {
        // If the node is destoyed, we make sure not to keep a dangling pointer to it
        Q_Q(QNode);
        auto f = [q, func, node]() { (static_cast<Caller *>(q)->*func)(node); };
        m_destructionConnections.push_back({node, QObject::connect(node, &QNode::nodeDestroyed, f)});
    }

    template<typename Caller, typename NodeType>
    void registerDestructionHelper(NodeType *node, DestructionFunctionPointer<Caller, NodeType> func, QList<NodeType*> &)
    {
        // If the node is destoyed, we make sure not to keep a dangling pointer to it
        Q_Q(QNode);
        auto f = [q, func, node]() { (static_cast<Caller *>(q)->*func)(node); };
        m_destructionConnections.push_back({node, QObject::connect(node, &QNode::nodeDestroyed, f)});
    }

    template<typename Caller, typename NodeType>
    void registerDestructionHelper(NodeType *node, DestructionFunctionPointer<Caller, NodeType> func, std::vector<NodeType*> &)
    {
        // If the node is destoyed, we make sure not to keep a dangling pointer to it
        Q_Q(QNode);
        auto f = [q, func, node]() { (static_cast<Caller *>(q)->*func)(node); };
        m_destructionConnections.push_back({node, QObject::connect(node, &QNode::nodeDestroyed, f)});
    }

    template<typename Caller, typename ValueType>
    using DestructionFunctionValue = void (Caller::*)(const ValueType&);

    template<typename Caller, typename NodeType, typename ValueType>
    void registerDestructionHelper(NodeType *node, DestructionFunctionValue<Caller, ValueType> func, NodeType *&,
                                   const ValueType &resetValue)
    {
        // If the node is destoyed, we make sure not to keep a dangling pointer to it
        Q_Q(QNode);
        auto f = [q, func, resetValue]() { (static_cast<Caller *>(q)->*func)(resetValue); };
        m_destructionConnections.push_back({node, QObject::connect(node, &QNode::nodeDestroyed, f)});
    }

    template<typename Caller, typename NodeType>
    void registerPrivateDestructionHelper(NodeType *node, DestructionFunctionPointer<Caller, NodeType> func)
    {
        // If the node is destoyed, we make sure not to keep a dangling pointer to it
        auto f = [this, func, node]() { (static_cast<Caller *>(this)->*func)(node); };
        m_destructionConnections.push_back({node, QObject::connect(node, &QNode::nodeDestroyed, f)});
    }

    void unregisterDestructionHelper(QNode *node)
    {
        m_destructionConnections.erase(std::remove_if(m_destructionConnections.begin(),
                                                      m_destructionConnections.end(),
                                                      [node] (const QPair<QNode *, QMetaObject::Connection> &nodeConnectionPair) {
                                                          if (nodeConnectionPair.first == node) {
                                                              QObject::disconnect(nodeConnectionPair.second);
                                                              return true;
                                                          }
                                                          return false;
                                                      }),
                                       m_destructionConnections.end());
    }

    static const QMetaObject *findStaticMetaObject(const QMetaObject *metaObject);

    void _q_postConstructorInit();
    void _q_ensureBackendNodeCreated();

private:
    void createBackendNode();
    void notifyDestructionChangesAndRemoveFromScene();
    void _q_addChild(QNode *childNode);
    void _q_removeChild(QNode *childNode);
    void _q_setParentHelper(QNode *parent);
    void registerNotifiedProperties();
    void unregisterNotifiedProperties();
    void propertyChanged(int propertyIndex);

    void setSceneHelper(QNode *root);
    void unsetSceneHelper(QNode *root);
    void addEntityComponentToScene(QNode *root);

    friend class PropertyChangeHandler<QNodePrivate>;
    bool m_propertyChangesSetup;
    PropertyChangeHandler<QNodePrivate> m_signals;
    QVector<QPair<QNode *, QMetaObject::Connection>> m_destructionConnections;
};

class NodePostConstructorInit : public QObject
{
    Q_OBJECT
public:
    NodePostConstructorInit(QObject *parent = nullptr);
    virtual ~NodePostConstructorInit();
    void removeNode(QNode *node);
    void addNode(QNode *node);

public Q_SLOTS:
    void processNodes();

private:
    QQueue<QNodePrivate *> m_nodesToConstruct;
    bool m_requestedProcessing;
};

} // namespace Qt3DCore

QT_END_NAMESPACE

#endif // QT3DCORE_NODE_P_H
