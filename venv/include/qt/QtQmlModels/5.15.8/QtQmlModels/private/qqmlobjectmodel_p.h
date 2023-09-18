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

#ifndef QQMLINSTANCEMODEL_P_H
#define QQMLINSTANCEMODEL_P_H

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

#include <private/qtqmlmodelsglobal_p.h>
#include <private/qqmlincubator_p.h>
#include <QtQml/qqml.h>
#include <QtCore/qobject.h>

QT_REQUIRE_CONFIG(qml_object_model);

QT_BEGIN_NAMESPACE

class QObject;
class QQmlChangeSet;
class QAbstractItemModel;

class Q_QMLMODELS_PRIVATE_EXPORT QQmlInstanceModel : public QObject
{
    Q_OBJECT

    Q_PROPERTY(int count READ count NOTIFY countChanged)
    QML_ANONYMOUS

public:
    enum ReusableFlag {
        NotReusable,
        Reusable
    };

    virtual ~QQmlInstanceModel() {}

    enum ReleaseFlag { Referenced = 0x01, Destroyed = 0x02, Pooled = 0x04 };
    Q_DECLARE_FLAGS(ReleaseFlags, ReleaseFlag)

    virtual int count() const = 0;
    virtual bool isValid() const = 0;
    virtual QObject *object(int index, QQmlIncubator::IncubationMode incubationMode = QQmlIncubator::AsynchronousIfNested) = 0;
    virtual ReleaseFlags release(QObject *object, ReusableFlag reusableFlag = NotReusable) = 0;
    virtual void cancel(int) {}
    QString stringValue(int index, const QString &role) { return variantValue(index, role).toString(); }
    virtual QVariant variantValue(int, const QString &) = 0;
    virtual void setWatchedRoles(const QList<QByteArray> &roles) = 0;
    virtual QQmlIncubator::Status incubationStatus(int index) = 0;

    virtual void drainReusableItemsPool(int maxPoolTime) { Q_UNUSED(maxPoolTime) }
    virtual int poolSize() { return 0; }

    virtual int indexOf(QObject *object, QObject *objectContext) const = 0;
    virtual const QAbstractItemModel *abstractItemModel() const { return nullptr; }

Q_SIGNALS:
    void countChanged();
    void modelUpdated(const QQmlChangeSet &changeSet, bool reset);
    void createdItem(int index, QObject *object);
    void initItem(int index, QObject *object);
    void destroyingItem(QObject *object);
    Q_REVISION(15) void itemPooled(int index, QObject *object);
    Q_REVISION(15) void itemReused(int index, QObject *object);

protected:
    QQmlInstanceModel(QObjectPrivate &dd, QObject *parent = nullptr)
        : QObject(dd, parent) {}

private:
    Q_DISABLE_COPY(QQmlInstanceModel)
};

class QQmlObjectModelAttached;
class QQmlObjectModelPrivate;
class Q_QMLMODELS_PRIVATE_EXPORT QQmlObjectModel : public QQmlInstanceModel
{
    Q_OBJECT
    Q_DECLARE_PRIVATE(QQmlObjectModel)

    Q_PROPERTY(QQmlListProperty<QObject> children READ children NOTIFY childrenChanged DESIGNABLE false)
    Q_CLASSINFO("DefaultProperty", "children")
    QML_NAMED_ELEMENT(ObjectModel)
    QML_ADDED_IN_MINOR_VERSION(1)
    QML_ATTACHED(QQmlObjectModelAttached)

public:
    QQmlObjectModel(QObject *parent=nullptr);
    ~QQmlObjectModel() {}

    int count() const override;
    bool isValid() const override;
    QObject *object(int index, QQmlIncubator::IncubationMode incubationMode = QQmlIncubator::AsynchronousIfNested) override;
    ReleaseFlags release(QObject *object, ReusableFlag reusable = NotReusable) override;
    QVariant variantValue(int index, const QString &role) override;
    void setWatchedRoles(const QList<QByteArray> &) override {}
    QQmlIncubator::Status incubationStatus(int index) override;

    int indexOf(QObject *object, QObject *objectContext) const override;

    QQmlListProperty<QObject> children();

    static QQmlObjectModelAttached *qmlAttachedProperties(QObject *obj);

    Q_REVISION(3) Q_INVOKABLE QObject *get(int index) const;
    Q_REVISION(3) Q_INVOKABLE void append(QObject *object);
    Q_REVISION(3) Q_INVOKABLE void insert(int index, QObject *object);
    Q_REVISION(3) Q_INVOKABLE void move(int from, int to, int n = 1);
    Q_REVISION(3) Q_INVOKABLE void remove(int index, int n = 1);

public Q_SLOTS:
    Q_REVISION(3) void clear();

Q_SIGNALS:
    void childrenChanged();

private:
    Q_DISABLE_COPY(QQmlObjectModel)
};

class QQmlObjectModelAttached : public QObject
{
    Q_OBJECT

public:
    QQmlObjectModelAttached(QObject *parent)
        : QObject(parent), m_index(-1) {}
    ~QQmlObjectModelAttached() {
        attachedProperties.remove(parent());
    }

    Q_PROPERTY(int index READ index NOTIFY indexChanged)
    int index() const { return m_index; }
    void setIndex(int idx) {
        if (m_index != idx) {
            m_index = idx;
            Q_EMIT indexChanged();
        }
    }

    static QQmlObjectModelAttached *properties(QObject *obj) {
        QQmlObjectModelAttached *rv = attachedProperties.value(obj);
        if (!rv) {
            rv = new QQmlObjectModelAttached(obj);
            attachedProperties.insert(obj, rv);
        }
        return rv;
    }

Q_SIGNALS:
    void indexChanged();

public:
    int m_index;

    static QHash<QObject*, QQmlObjectModelAttached*> attachedProperties;
};


QT_END_NAMESPACE

QML_DECLARE_TYPE(QQmlInstanceModel)
QML_DECLARE_TYPE(QQmlObjectModel)

#endif // QQMLINSTANCEMODEL_P_H
