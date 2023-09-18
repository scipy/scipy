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

#ifndef QT3DCORE_QNODE_H
#define QT3DCORE_QNODE_H

#include <Qt3DCore/qnodecreatedchange.h>
#include <Qt3DCore/qnodeid.h>
#include <Qt3DCore/qnodecommand.h>
#include <Qt3DCore/qscenechange.h>
#include <Qt3DCore/qt3dcore_global.h>
#include <QtCore/QObject>

#define Q_NODE_NULLPTR static_cast<Qt3DCore::QNode *>(nullptr)

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

class QNode;
class QNodePrivate;
class QEntity;
class QAspectEngine;

#if defined(QT_BUILD_INTERNAL)
class QBackendNodeTester;
#endif

typedef QVector<QNode *> QNodeVector;
typedef QSharedPointer<QNode> QNodePtr;

class Q_3DCORESHARED_EXPORT QNode : public QObject
{
    Q_OBJECT
    Q_PROPERTY(Qt3DCore::QNode *parent READ parentNode WRITE setParent NOTIFY parentChanged)
    Q_PROPERTY(bool enabled READ isEnabled WRITE setEnabled NOTIFY enabledChanged)
    Q_PROPERTY(PropertyTrackingMode defaultPropertyTrackingMode READ defaultPropertyTrackingMode WRITE setDefaultPropertyTrackingMode NOTIFY defaultPropertyTrackingModeChanged REVISION 9)
public:

    enum PropertyTrackingMode : quint16 {
        TrackFinalValues,
        DontTrackValues,
        TrackAllValues
    };
    Q_ENUM(PropertyTrackingMode)

    explicit QNode(QNode *parent = nullptr);
    virtual ~QNode();

    QNodeId id() const;
    QNode *parentNode() const;

    bool notificationsBlocked() const;
    bool blockNotifications(bool block);

    QNodeVector childNodes() const;

    bool isEnabled() const;
    PropertyTrackingMode defaultPropertyTrackingMode() const;

    void setPropertyTracking(const QString &propertyName, PropertyTrackingMode trackMode);
    PropertyTrackingMode propertyTracking(const QString &propertyName) const;
    void clearPropertyTracking(const QString &propertyName);
    void clearPropertyTrackings();

    Q3D_DECL_DEPRECATED QNodeCommand::CommandId sendCommand(const QString &name, const QVariant &data = QVariant(),
                                                          QNodeCommand::CommandId replyTo = QNodeCommand::CommandId());
    Q3D_DECL_DEPRECATED void sendReply(const QNodeCommandPtr &command);

public Q_SLOTS:
    void setParent(QNode *parent);
    void setEnabled(bool isEnabled);
    void setDefaultPropertyTrackingMode(PropertyTrackingMode mode);

Q_SIGNALS:
    void parentChanged(QObject *parent);
    void enabledChanged(bool enabled);
    void defaultPropertyTrackingModeChanged(PropertyTrackingMode mode);
    void nodeDestroyed();

protected:
    explicit QNode(QNodePrivate &dd, QNode *parent = nullptr);
    Q3D_DECL_DEPRECATED void notifyObservers(const QSceneChangePtr &change);
    Q3D_DECL_DEPRECATED virtual void sceneChangeEvent(const QSceneChangePtr &change);

private:
    Q_DECLARE_PRIVATE(QNode)
    Q3D_DECL_DEPRECATED virtual QNodeCreatedChangeBasePtr createNodeCreationChange() const;

    // We only want setParent(QNode *) to be callable
    // when dealing with QNode objects
    void setParent(QObject *) Q_DECL_EQ_DELETE;

    Q_PRIVATE_SLOT(d_func(), void _q_postConstructorInit())
    Q_PRIVATE_SLOT(d_func(), void _q_addChild(Qt3DCore::QNode *))
    Q_PRIVATE_SLOT(d_func(), void _q_removeChild(Qt3DCore::QNode *))
    Q_PRIVATE_SLOT(d_func(), void _q_setParentHelper(Qt3DCore::QNode *))

    friend class QAspectEngine;
    friend class QAspectEnginePrivate;
    friend class QAbstractAspectPrivate;
    friend class QNodeCreatedChangeGenerator;
    friend class QPostman;
    friend class QScene;

#if defined(QT_BUILD_INTERNAL)
    friend class QBackendNodeTester;
#endif
};

inline QNodeId qIdForNode(QNode *node){ return node ? node->id() : QNodeId(); }

template<typename T>
inline QNodeIdVector qIdsForNodes(const T &nodes)
{
    QNodeIdVector ids;
    ids.reserve(nodes.size());
    for (const auto n : nodes)
        ids.push_back(n->id());
    return ids;
}

struct QNodeIdTypePair
{
    QNodeIdTypePair() Q_DECL_NOTHROW
        : id()
        , type(nullptr)
    {}

    explicit QNodeIdTypePair(QNodeId _id, const QMetaObject *_type) Q_DECL_NOTHROW
        : id(_id)
        , type(_type)
    {}

    QNodeId id;
    const QMetaObject *type;
};
QT3D_DECLARE_TYPEINFO(Qt3DCore, QNodeIdTypePair, Q_PRIMITIVE_TYPE)

} // namespace Qt3DCore

QT_END_NAMESPACE

#endif
