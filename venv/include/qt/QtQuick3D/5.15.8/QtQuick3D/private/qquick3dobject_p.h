/****************************************************************************
**
** Copyright (C) 2019 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of Qt Quick 3D.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QSSGOBJECT_P_H
#define QSSGOBJECT_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include "qquick3dobject.h"

#include "qtquick3dglobal_p.h"

#include "qquick3dobjectchangelistener_p.h"

#include <private/qobject_p.h>
#include <private/qquickstate_p.h>
#include <private/qqmlnotifier_p.h>
#include <private/qlazilyallocated_p.h>
#include <private/qssgrendergraphobject_p.h>

QT_BEGIN_NAMESPACE

class QSSGRenderContextInterface;
class QSSGRenderContext;

class Q_QUICK3D_PRIVATE_EXPORT QQuick3DObjectPrivate : public QObjectPrivate
{
    Q_DECLARE_PUBLIC(QQuick3DObject)
public:
    enum class Type {
        Unknown = 0,
        Node, // Node
        Layer, // Node
        Light, // Node
        Camera, // Node
        Model, // Node
        Text, // Node
        Item2D, // Renderable? Node
        SceneEnvironment, // Resource
        DefaultMaterial, // Resource
        PrincipledMaterial, // Resource
        Image, // Resource
        Effect, // Resource
        CustomMaterial, // Resource
        Lightmaps, // Resource
        Geometry, // Resource
        RenderPlugin, // Not used
        LastKnownGraphObjectType,
    };

    static QQuick3DObjectPrivate *get(QQuick3DObject *item) { return item->d_func(); }
    static const QQuick3DObjectPrivate *get(const QQuick3DObject *item) { return item->d_func(); }

    explicit QQuick3DObjectPrivate(Type t);
    ~QQuick3DObjectPrivate() override;
    void init(QQuick3DObject *parent);

    QQmlListProperty<QObject> data();
    QQmlListProperty<QObject> resources();
    QQmlListProperty<QQuick3DObject> children();

    QQmlListProperty<QQuickState> states();
    QQmlListProperty<QQuickTransition> transitions();

    QString state() const;
    void setState(const QString &);

    // data property
    static void data_append(QQmlListProperty<QObject> *, QObject *);
    static int data_count(QQmlListProperty<QObject> *);
    static QObject *data_at(QQmlListProperty<QObject> *, int);
    static void data_clear(QQmlListProperty<QObject> *);

    // resources property
    static QObject *resources_at(QQmlListProperty<QObject> *, int);
    static void resources_append(QQmlListProperty<QObject> *, QObject *);
    static int resources_count(QQmlListProperty<QObject> *);
    static void resources_clear(QQmlListProperty<QObject> *);

    // children property
    static void children_append(QQmlListProperty<QQuick3DObject> *, QQuick3DObject *);
    static int children_count(QQmlListProperty<QQuick3DObject> *);
    static QQuick3DObject *children_at(QQmlListProperty<QQuick3DObject> *, int);
    static void children_clear(QQmlListProperty<QQuick3DObject> *);

    void _q_resourceObjectDeleted(QObject *);
    quint64 _q_createJSWrapper(QV4::ExecutionEngine *engine);

    enum ChangeType {
        Geometry = 0x01,
        SiblingOrder = 0x02,
        Visibility = 0x04,
        Opacity = 0x08,
        Destroyed = 0x10,
        Parent = 0x20,
        Children = 0x40,
        Rotation = 0x80,
        ImplicitWidth = 0x100,
        ImplicitHeight = 0x200,
        Enabled = 0x400,
    };

    Q_DECLARE_FLAGS(ChangeTypes, ChangeType)

    struct ChangeListener
    {
        using ChangeTypes = QQuick3DObjectPrivate::ChangeTypes;

        ChangeListener(QQuick3DObjectChangeListener *l = nullptr, ChangeTypes t = {}) : listener(l), types(t) {}

        ChangeListener(QQuick3DObjectChangeListener *l) : listener(l), types(Geometry) {}

        bool operator==(const ChangeListener &other) const
        {
            return listener == other.listener && types == other.types;
        }

        QQuick3DObjectChangeListener *listener;
        ChangeTypes types;

        QVector<QQuick3DObjectPrivate::ChangeListener> changeListeners;
    };

    struct ExtraData
    {
        ExtraData();

        int hideRefCount;
        QObjectList resourcesList;

    };
    QLazilyAllocated<ExtraData> extra;

    QVector<QQuick3DObjectPrivate::ChangeListener> changeListeners;

    void addItemChangeListener(QQuick3DObjectChangeListener *listener, ChangeTypes types);
    void updateOrAddItemChangeListener(QQuick3DObjectChangeListener *listener, ChangeTypes types);
    void removeItemChangeListener(QQuick3DObjectChangeListener *, ChangeTypes types);

    QQuickStateGroup *_states();
    QQuickStateGroup *_stateGroup;

    enum DirtyType {
        TransformOrigin = 0x00000001,
        Transform = 0x00000002,
        BasicTransform = 0x00000004,
        Position = 0x00000008,
        Size = 0x00000010,

        ZValue = 0x00000020,
        Content = 0x00000040,
        Smooth = 0x00000080,
        OpacityValue = 0x00000100,
        ChildrenChanged = 0x00000200,
        ChildrenStackingChanged = 0x00000400,
        ParentChanged = 0x00000800,

        Clip = 0x00001000,
        Window = 0x00002000,

        EffectReference = 0x00008000,
        Visible = 0x00010000,
        HideReference = 0x00020000,
        Antialiasing = 0x00040000,
        // When you add an attribute here, don't forget to update
        // dirtyToString()

        TransformUpdateMask = TransformOrigin | Transform | BasicTransform | Position | Window,
        ComplexTransformUpdateMask = Transform | Window,
        ContentUpdateMask = Size | Content | Smooth | Window | Antialiasing,
        ChildrenUpdateMask = ChildrenChanged | ChildrenStackingChanged | EffectReference | Window
    };

    quint32 dirtyAttributes;
    QString dirtyToString() const;
    void dirty(DirtyType);
    void addToDirtyList();
    void removeFromDirtyList();
    QQuick3DObject *nextDirtyItem;
    QQuick3DObject **prevDirtyItem;

    bool isResourceNode() const;
    bool isSpatialNode() const;

    void setCulled(bool);

    QSharedPointer<QQuick3DSceneManager> sceneManager;
    int windowRefCount;

    QQuick3DObject *parentItem;

    QList<QQuick3DObject *> childItems;
    mutable QList<QQuick3DObject *> *sortedChildItems;
    QList<QQuick3DObject *> paintOrderChildItems() const;
    void addChild(QQuick3DObject *);
    void removeChild(QQuick3DObject *);
    void siblingOrderChanged();

    void markSortedChildrenDirty(QQuick3DObject *child);

    void refSceneManager(const QSharedPointer<QQuick3DSceneManager> &);
    void derefSceneManager();

    static void refSceneManager(QQuick3DObject *obj,const QSharedPointer<QQuick3DSceneManager> &mgr)
    {
        if (obj)
            QQuick3DObjectPrivate::get(obj)->refSceneManager(mgr);
    }
    static void derefSceneManager(QQuick3DObject *obj)
    {
        if (obj)
            QQuick3DObjectPrivate::get(obj)->derefSceneManager();
    }

    QQuick3DObject *subFocusItem;
    void updateSubFocusItem(QQuick3DObject *scope, bool focus);

    void itemChange(QQuick3DObject::ItemChange, const QQuick3DObject::ItemChangeData &);

    virtual void updatePolish() {}

    QSSGRenderGraphObject *spatialNode = nullptr;

    Type type = Type::Unknown;
    bool componentComplete = true;
    bool culled;
};

Q_DECLARE_OPERATORS_FOR_FLAGS(QQuick3DObjectPrivate::ChangeTypes)
Q_DECLARE_TYPEINFO(QQuick3DObjectPrivate::ChangeListener, Q_PRIMITIVE_TYPE);

QT_END_NAMESPACE

#endif // QSSGOBJECT_P_H
