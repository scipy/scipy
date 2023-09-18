/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtQml module of the Qt Toolkit.
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

#ifndef QQMLDATAMODEL_P_P_H
#define QQMLDATAMODEL_P_P_H

#include "qqmldelegatemodel_p.h"
#include <private/qv4qobjectwrapper_p.h>

#include <QtQml/qqmlcontext.h>
#include <QtQml/qqmlincubator.h>

#include <private/qqmladaptormodel_p.h>
#include <private/qqmlopenmetaobject_p.h>

#include <QtCore/qloggingcategory.h>

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

QT_REQUIRE_CONFIG(qml_delegate_model);

QT_BEGIN_NAMESPACE

Q_DECLARE_LOGGING_CATEGORY(lcItemViewDelegateRecycling)

typedef QQmlListCompositor Compositor;

class QQmlDelegateModelAttachedMetaObject;
class QQmlAbstractDelegateComponent;

class Q_QMLMODELS_PRIVATE_EXPORT QQmlDelegateModelItemMetaType : public QQmlRefCount
{
public:
    QQmlDelegateModelItemMetaType(QV4::ExecutionEngine *engine, QQmlDelegateModel *model, const QStringList &groupNames);
    ~QQmlDelegateModelItemMetaType();

    void initializeMetaObject();
    void initializePrototype();

    int parseGroups(const QStringList &groupNames) const;
    int parseGroups(const QV4::Value &groupNames) const;

    QPointer<QQmlDelegateModel> model;
    const int groupCount;
    QV4::ExecutionEngine * const v4Engine;
    QQmlDelegateModelAttachedMetaObject *metaObject;
    const QStringList groupNames;
    QV4::PersistentValue modelItemProto;
};

class QQmlAdaptorModel;
class QQDMIncubationTask;

class QQmlDelegateModelItem : public QObject
{
    Q_OBJECT
    Q_PROPERTY(int index READ modelIndex NOTIFY modelIndexChanged)
    Q_PROPERTY(int row READ modelRow NOTIFY rowChanged REVISION 12)
    Q_PROPERTY(int column READ modelColumn NOTIFY columnChanged REVISION 12)
    Q_PROPERTY(QObject *model READ modelObject CONSTANT)
public:
    QQmlDelegateModelItem(const QQmlRefPointer<QQmlDelegateModelItemMetaType> &metaType,
                          QQmlAdaptorModel::Accessors *accessor, int modelIndex,
                          int row, int column);
    ~QQmlDelegateModelItem();

    void referenceObject() { ++objectRef; }
    bool releaseObject() { return --objectRef == 0 && !(groups & Compositor::PersistedFlag); }
    bool isObjectReferenced() const { return objectRef != 0 || (groups & Compositor::PersistedFlag); }
    void childContextObjectDestroyed(QObject *childContextObject);

    bool isReferenced() const {
        return scriptRef
                || incubationTask
                || ((groups & Compositor::UnresolvedFlag) && (groups & Compositor::GroupMask));
    }

    void Dispose();

    QObject *modelObject() { return this; }

    void destroyObject();

    static QQmlDelegateModelItem *dataForObject(QObject *object);

    int groupIndex(Compositor::Group group);

    int modelRow() const { return row; }
    int modelColumn() const { return column; }
    int modelIndex() const { return index; }
    virtual void setModelIndex(int idx, int newRow, int newColumn, bool alwaysEmit = false);

    virtual QV4::ReturnedValue get() { return QV4::QObjectWrapper::wrap(v4, this); }

    virtual void setValue(const QString &role, const QVariant &value) { Q_UNUSED(role); Q_UNUSED(value); }
    virtual bool resolveIndex(const QQmlAdaptorModel &, int) { return false; }

    static QV4::ReturnedValue get_model(const QV4::FunctionObject *, const QV4::Value *thisObject, const QV4::Value *argv, int argc);
    static QV4::ReturnedValue get_groups(const QV4::FunctionObject *, const QV4::Value *thisObject, const QV4::Value *argv, int argc);
    static QV4::ReturnedValue set_groups(const QV4::FunctionObject *, const QV4::Value *thisObject, const QV4::Value *argv, int argc);
    static QV4::ReturnedValue get_member(QQmlDelegateModelItem *thisItem, uint flag, const QV4::Value &);
    static QV4::ReturnedValue set_member(QQmlDelegateModelItem *thisItem, uint flag, const QV4::Value &arg);
    static QV4::ReturnedValue get_index(QQmlDelegateModelItem *thisItem, uint flag, const QV4::Value &arg);

    QV4::ExecutionEngine *v4;
    QQmlRefPointer<QQmlDelegateModelItemMetaType> const metaType;
    QQmlContextDataRef contextData;
    QPointer<QObject> object;
    QPointer<QQmlDelegateModelAttached> attached;
    QQDMIncubationTask *incubationTask;
    QQmlComponent *delegate;
    int poolTime;
    int objectRef;
    int scriptRef;
    int groups;
    int index;

Q_SIGNALS:
    void modelIndexChanged();
    Q_REVISION(12) void rowChanged();
    Q_REVISION(12) void columnChanged();

protected:
    void objectDestroyed(QObject *);
    int row;
    int column;
};

namespace QV4 {
namespace Heap {
struct QQmlDelegateModelItemObject : Object {
    inline void init(QQmlDelegateModelItem *item);
    void destroy();
    QQmlDelegateModelItem *item;
};

}
}

struct QQmlDelegateModelItemObject : QV4::Object
{
    V4_OBJECT2(QQmlDelegateModelItemObject, QV4::Object)
    V4_NEEDS_DESTROY
};

void QV4::Heap::QQmlDelegateModelItemObject::init(QQmlDelegateModelItem *item)
{
    Object::init();
    this->item = item;
}

class QQmlReusableDelegateModelItemsPool
{
public:
    void insertItem(QQmlDelegateModelItem *modelItem);
    QQmlDelegateModelItem *takeItem(const QQmlComponent *delegate, int newIndexHint);
    void reuseItem(QQmlDelegateModelItem *item, int newModelIndex);
    void drain(int maxPoolTime, std::function<void(QQmlDelegateModelItem *cacheItem)> releaseItem);
    int size() { return m_reusableItemsPool.size(); }

private:
    QList<QQmlDelegateModelItem *> m_reusableItemsPool;
};

class QQmlDelegateModelPrivate;
class QQDMIncubationTask : public QQmlIncubator
{
public:
    QQDMIncubationTask(QQmlDelegateModelPrivate *l, IncubationMode mode)
        : QQmlIncubator(mode)
        , incubating(nullptr)
        , vdm(l) {}

    void initializeRequiredProperties(QQmlDelegateModelItem *modelItemToIncubate, QObject* object);
    void statusChanged(Status) override;
    void setInitialState(QObject *) override;

    QQmlDelegateModelItem *incubating = nullptr;
    QQmlDelegateModelPrivate *vdm = nullptr;
    QQmlContextData *proxyContext = nullptr;
    QPointer<QObject> proxiedObject  = nullptr; // the proxied object might disapear, so we use a QPointer instead of a raw one
    int index[QQmlListCompositor::MaximumGroupCount];
};


class QQmlDelegateModelGroupEmitter
{
public:
    virtual ~QQmlDelegateModelGroupEmitter() {}
    virtual void emitModelUpdated(const QQmlChangeSet &changeSet, bool reset) = 0;
    virtual void createdPackage(int, QQuickPackage *) {}
    virtual void initPackage(int, QQuickPackage *) {}
    virtual void destroyingPackage(QQuickPackage *) {}

    QIntrusiveListNode emitterNode;
};

typedef QIntrusiveList<QQmlDelegateModelGroupEmitter, &QQmlDelegateModelGroupEmitter::emitterNode> QQmlDelegateModelGroupEmitterList;

class QQmlDelegateModelGroupPrivate : public QObjectPrivate
{
public:
    Q_DECLARE_PUBLIC(QQmlDelegateModelGroup)

    QQmlDelegateModelGroupPrivate() : group(Compositor::Cache), defaultInclude(false) {}

    static QQmlDelegateModelGroupPrivate *get(QQmlDelegateModelGroup *group) {
        return static_cast<QQmlDelegateModelGroupPrivate *>(QObjectPrivate::get(group)); }

    void setModel(QQmlDelegateModel *model, Compositor::Group group);
    bool isChangedConnected();
    void emitChanges(QV4::ExecutionEngine *engine);
    void emitModelUpdated(bool reset);

    void createdPackage(int index, QQuickPackage *package);
    void initPackage(int index, QQuickPackage *package);
    void destroyingPackage(QQuickPackage *package);

    bool parseIndex(const QV4::Value &value, int *index, Compositor::Group *group) const;
    bool parseGroupArgs(
            QQmlV4Function *args, Compositor::Group *group, int *index, int *count, int *groups) const;

    Compositor::Group group;
    QPointer<QQmlDelegateModel> model;
    QQmlDelegateModelGroupEmitterList emitters;
    QQmlChangeSet changeSet;
    QString name;
    bool defaultInclude;
};

class QQmlDelegateModelParts;

class QQmlDelegateModelPrivate : public QObjectPrivate, public QQmlDelegateModelGroupEmitter
{
    Q_DECLARE_PUBLIC(QQmlDelegateModel)
public:
    QQmlDelegateModelPrivate(QQmlContext *);
    ~QQmlDelegateModelPrivate();

    static QQmlDelegateModelPrivate *get(QQmlDelegateModel *m) {
        return static_cast<QQmlDelegateModelPrivate *>(QObjectPrivate::get(m));
    }

    void init();
    void connectModel(QQmlAdaptorModel *model);
    void connectToAbstractItemModel();
    void disconnectFromAbstractItemModel();

    void requestMoreIfNecessary();
    QObject *object(Compositor::Group group, int index, QQmlIncubator::IncubationMode incubationMode);
    QQmlDelegateModel::ReleaseFlags release(QObject *object, QQmlInstanceModel::ReusableFlag reusable = QQmlInstanceModel::NotReusable);
    QVariant variantValue(Compositor::Group group, int index, const QString &name);
    void emitCreatedPackage(QQDMIncubationTask *incubationTask, QQuickPackage *package);
    void emitInitPackage(QQDMIncubationTask *incubationTask, QQuickPackage *package);
    void emitCreatedItem(QQDMIncubationTask *incubationTask, QObject *item) {
        Q_EMIT q_func()->createdItem(incubationTask->index[m_compositorGroup], item); }
    void emitInitItem(QQDMIncubationTask *incubationTask, QObject *item) {
        Q_EMIT q_func()->initItem(incubationTask->index[m_compositorGroup], item); }
    void emitDestroyingPackage(QQuickPackage *package);
    void emitDestroyingItem(QObject *item) { Q_EMIT q_func()->destroyingItem(item); }
    void addCacheItem(QQmlDelegateModelItem *item, Compositor::iterator it);
    void removeCacheItem(QQmlDelegateModelItem *cacheItem);
    void destroyCacheItem(QQmlDelegateModelItem *cacheItem);
    void updateFilterGroup();

    void reuseItem(QQmlDelegateModelItem *item, int newModelIndex, int newGroups);
    void drainReusableItemsPool(int maxPoolTime);
    QQmlComponent *resolveDelegate(int index);

    void addGroups(Compositor::iterator from, int count, Compositor::Group group, int groupFlags);
    void removeGroups(Compositor::iterator from, int count, Compositor::Group group, int groupFlags);
    void setGroups(Compositor::iterator from, int count, Compositor::Group group, int groupFlags);

    void itemsInserted(
            const QVector<Compositor::Insert> &inserts,
            QVarLengthArray<QVector<QQmlChangeSet::Change>, Compositor::MaximumGroupCount> *translatedInserts,
            QHash<int, QList<QQmlDelegateModelItem *> > *movedItems = nullptr);
    void itemsInserted(const QVector<Compositor::Insert> &inserts);
    void itemsRemoved(
            const QVector<Compositor::Remove> &removes,
            QVarLengthArray<QVector<QQmlChangeSet::Change>, Compositor::MaximumGroupCount> *translatedRemoves,
            QHash<int, QList<QQmlDelegateModelItem *> > *movedItems = nullptr);
    void itemsRemoved(const QVector<Compositor::Remove> &removes);
    void itemsMoved(
            const QVector<Compositor::Remove> &removes, const QVector<Compositor::Insert> &inserts);
    void itemsChanged(const QVector<Compositor::Change> &changes);
    void emitChanges();
    void emitModelUpdated(const QQmlChangeSet &changeSet, bool reset) override;
    void delegateChanged(bool add = true, bool remove = true);

    bool insert(Compositor::insert_iterator &before, const QV4::Value &object, int groups);

    int adaptorModelCount() const;

    static void group_append(QQmlListProperty<QQmlDelegateModelGroup> *property, QQmlDelegateModelGroup *group);
    static int group_count(QQmlListProperty<QQmlDelegateModelGroup> *property);
    static QQmlDelegateModelGroup *group_at(QQmlListProperty<QQmlDelegateModelGroup> *property, int index);

    void releaseIncubator(QQDMIncubationTask *incubationTask);
    void incubatorStatusChanged(QQDMIncubationTask *incubationTask, QQmlIncubator::Status status);
    void setInitialState(QQDMIncubationTask *incubationTask, QObject *o);

    QQmlAdaptorModel m_adaptorModel;
    QQmlListCompositor m_compositor;
    QQmlStrongJSQObjectReference<QQmlComponent> m_delegate;
    QQmlAbstractDelegateComponent *m_delegateChooser;
    QMetaObject::Connection m_delegateChooserChanged;
    QQmlDelegateModelItemMetaType *m_cacheMetaType;
    QPointer<QQmlContext> m_context;
    QQmlDelegateModelParts *m_parts;
    QQmlDelegateModelGroupEmitterList m_pendingParts;

    QList<QQmlDelegateModelItem *> m_cache;
    QQmlReusableDelegateModelItemsPool m_reusableItemsPool;
    QList<QQDMIncubationTask *> m_finishedIncubating;
    QList<QByteArray> m_watchedRoles;

    QString m_filterGroup;

    int m_count;
    int m_groupCount;

    QQmlListCompositor::Group m_compositorGroup;
    bool m_complete : 1;
    bool m_delegateValidated : 1;
    bool m_reset : 1;
    bool m_transaction : 1;
    bool m_incubatorCleanupScheduled : 1;
    bool m_waitingToFetchMore : 1;

    union {
        struct {
            QQmlDelegateModelGroup *m_cacheItems;
            QQmlDelegateModelGroup *m_items;
            QQmlDelegateModelGroup *m_persistedItems;
        };
        QQmlDelegateModelGroup *m_groups[Compositor::MaximumGroupCount];
    };
};

class QQmlPartsModel : public QQmlInstanceModel, public QQmlDelegateModelGroupEmitter
{
    Q_OBJECT
    Q_PROPERTY(QString filterOnGroup READ filterGroup WRITE setFilterGroup NOTIFY filterGroupChanged RESET resetFilterGroup)
public:
    QQmlPartsModel(QQmlDelegateModel *model, const QString &part, QObject *parent = nullptr);
    ~QQmlPartsModel();

    QString filterGroup() const;
    void setFilterGroup(const QString &group);
    void resetFilterGroup();
    void updateFilterGroup();
    void updateFilterGroup(Compositor::Group group, const QQmlChangeSet &changeSet);

    int count() const override;
    bool isValid() const override;
    QObject *object(int index, QQmlIncubator::IncubationMode incubationMode = QQmlIncubator::AsynchronousIfNested) override;
    ReleaseFlags release(QObject *item, ReusableFlag reusable = NotReusable) override;
    QVariant variantValue(int index, const QString &role) override;
    QList<QByteArray> watchedRoles() const { return m_watchedRoles; }
    void setWatchedRoles(const QList<QByteArray> &roles) override;
    QQmlIncubator::Status incubationStatus(int index) override;

    int indexOf(QObject *item, QObject *objectContext) const override;

    void emitModelUpdated(const QQmlChangeSet &changeSet, bool reset) override;

    void createdPackage(int index, QQuickPackage *package) override;
    void initPackage(int index, QQuickPackage *package) override;
    void destroyingPackage(QQuickPackage *package) override;

Q_SIGNALS:
    void filterGroupChanged();

private:
    QQmlDelegateModel *m_model;
    QMultiHash<QObject *, QQuickPackage *> m_packaged;
    QString m_part;
    QString m_filterGroup;
    QList<QByteArray> m_watchedRoles;
    QVector<int> m_pendingPackageInitializations; // vector holds model indices
    Compositor::Group m_compositorGroup;
    bool m_inheritGroup;
    bool m_modelUpdatePending = true;
};

class QMetaPropertyBuilder;

class QQmlDelegateModelPartsMetaObject : public QQmlOpenMetaObject
{
public:
    QQmlDelegateModelPartsMetaObject(QObject *parent)
    : QQmlOpenMetaObject(parent) {}

    void propertyCreated(int, QMetaPropertyBuilder &) override;
    QVariant initialValue(int) override;
};

class QQmlDelegateModelParts : public QObject
{
Q_OBJECT
public:
    QQmlDelegateModelParts(QQmlDelegateModel *parent);

    QQmlDelegateModel *model;
    QList<QQmlPartsModel *> models;
};

class QQmlDelegateModelAttachedMetaObject : public QAbstractDynamicMetaObject, public QQmlRefCount
{
public:
    QQmlDelegateModelAttachedMetaObject(
            QQmlDelegateModelItemMetaType *metaType, QMetaObject *metaObject);
    ~QQmlDelegateModelAttachedMetaObject();

    void objectDestroyed(QObject *) override;
    int metaCall(QObject *, QMetaObject::Call, int _id, void **) override;

private:
    QQmlDelegateModelItemMetaType * const metaType;
    QMetaObject * const metaObject;
    const int memberPropertyOffset;
    const int indexPropertyOffset;
};

class PropertyUpdater : public QObject
{
    Q_OBJECT

public:
    PropertyUpdater(QObject *parent);
    QHash<int, QMetaObject::Connection> senderToConnection;
    QHash<int, int> changeSignalIndexToPropertyIndex;
    int updateCount = 0;
public Q_SLOTS:
    void doUpdate();
    void breakBinding();
};

QT_END_NAMESPACE

#endif
