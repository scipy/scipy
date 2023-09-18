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

#ifndef QQMLADAPTORMODEL_P_H
#define QQMLADAPTORMODEL_P_H

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

#include <QtCore/qabstractitemmodel.h>

#include <private/qqmlglobal_p.h>
#include <private/qqmllistaccessor_p.h>
#include <private/qtqmlmodelsglobal_p.h>
#include <private/qqmlguard_p.h>
#include <private/qqmlnullablevalue_p.h>
#include <private/qqmlpropertycache_p.h>

QT_REQUIRE_CONFIG(qml_delegate_model);

QT_BEGIN_NAMESPACE

class QQmlEngine;

class QQmlDelegateModel;
class QQmlDelegateModelItem;
class QQmlDelegateModelItemMetaType;

class Q_QMLMODELS_PRIVATE_EXPORT QQmlAdaptorModel : public QQmlStrongJSQObjectReference<QObject>
{
public:
    class Accessors
    {
    public:
        inline Accessors() {}
        virtual ~Accessors();
        virtual int rowCount(const QQmlAdaptorModel &) const { return 0; }
        virtual int columnCount(const QQmlAdaptorModel &) const { return 0; }
        virtual void cleanup(QQmlAdaptorModel &) const {}

        virtual QVariant value(const QQmlAdaptorModel &, int, const QString &) const {
            return QVariant(); }

        virtual QQmlDelegateModelItem *createItem(
                QQmlAdaptorModel &,
                const QQmlRefPointer<QQmlDelegateModelItemMetaType> &,
                int, int, int) const { return nullptr; }

        virtual bool notify(
                const QQmlAdaptorModel &,
                const QList<QQmlDelegateModelItem *> &,
                int,
                int,
                const QVector<int> &) const { return false; }
        virtual void replaceWatchedRoles(
                QQmlAdaptorModel &,
                const QList<QByteArray> &,
                const QList<QByteArray> &) const {}
        virtual QVariant parentModelIndex(const QQmlAdaptorModel &) const {
            return QVariant(); }
        virtual QVariant modelIndex(const QQmlAdaptorModel &, int) const {
            return QVariant(); }
        virtual bool canFetchMore(const QQmlAdaptorModel &) const { return false; }
        virtual void fetchMore(QQmlAdaptorModel &) const {}

        QScopedPointer<QMetaObject, QScopedPointerPodDeleter> metaObject;
        QQmlRefPointer<QQmlPropertyCache> propertyCache;
    };

    const Accessors *accessors;
    QPersistentModelIndex rootIndex;
    QQmlListAccessor list;

    int modelItemRevision = 0;

    QQmlAdaptorModel();
    ~QQmlAdaptorModel();

    inline QVariant model() const { return list.list(); }
    void setModel(const QVariant &variant, QObject *parent, QQmlEngine *engine);
    void invalidateModel();

    bool isValid() const;
    int count() const;
    int rowCount() const;
    int columnCount() const;
    int rowAt(int index) const;
    int columnAt(int index) const;
    int indexAt(int row, int column) const;

    void useImportVersion(int minorVersion);

    inline bool adaptsAim() const { return qobject_cast<QAbstractItemModel *>(object()); }
    inline QAbstractItemModel *aim() { return static_cast<QAbstractItemModel *>(object()); }
    inline const QAbstractItemModel *aim() const { return static_cast<const QAbstractItemModel *>(object()); }

    inline QVariant value(int index, const QString &role) const {
        return accessors->value(*this, index, role); }
    inline QQmlDelegateModelItem *createItem(
            const QQmlRefPointer<QQmlDelegateModelItemMetaType> &metaType, int index)
    {
        return accessors->createItem(*this, metaType, index, rowAt(index), columnAt(index));
    }
    inline bool hasProxyObject() const {
        return list.type() == QQmlListAccessor::Instance
                || list.type() == QQmlListAccessor::ListProperty
                || list.type() == QQmlListAccessor::ObjectList;
    }

    inline bool notify(
            const QList<QQmlDelegateModelItem *> &items,
            int index,
            int count,
            const QVector<int> &roles) const {
        return accessors->notify(*this, items, index, count, roles); }
    inline void replaceWatchedRoles(
            const QList<QByteArray> &oldRoles, const QList<QByteArray> &newRoles) {
        accessors->replaceWatchedRoles(*this, oldRoles, newRoles); }

    inline QVariant modelIndex(int index) const { return accessors->modelIndex(*this, index); }
    inline QVariant parentModelIndex() const { return accessors->parentModelIndex(*this); }
    inline bool canFetchMore() const { return accessors->canFetchMore(*this); }
    inline void fetchMore() { return accessors->fetchMore(*this); }

protected:
    void objectDestroyed(QObject *) override;
};

class QQmlAdaptorModelProxyInterface
{
public:
    virtual ~QQmlAdaptorModelProxyInterface() {}

    virtual QObject *proxiedObject() = 0;
};

#define QQmlAdaptorModelProxyInterface_iid "org.qt-project.Qt.QQmlAdaptorModelProxyInterface"

Q_DECLARE_INTERFACE(QQmlAdaptorModelProxyInterface, QQmlAdaptorModelProxyInterface_iid)

QT_END_NAMESPACE

#endif
