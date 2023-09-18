/****************************************************************************
**
** Copyright (C) 2017 Ford Motor Company
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtRemoteObjects module of the Qt Toolkit.
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

#ifndef QREMOTEOBJECTS_ABSTRACTITEMMODELREPLICA_H
#define QREMOTEOBJECTS_ABSTRACTITEMMODELREPLICA_H

#include <QtRemoteObjects/qtremoteobjectglobal.h>

#include <QtCore/qabstractitemmodel.h>
#include <QtCore/qitemselectionmodel.h>

QT_BEGIN_NAMESPACE

class QAbstractItemModelReplicaImplementation;

class Q_REMOTEOBJECTS_EXPORT QAbstractItemModelReplica : public QAbstractItemModel
{
    Q_OBJECT
public:
    ~QAbstractItemModelReplica() override;

    QItemSelectionModel* selectionModel() const;

    QVariant data(const QModelIndex & index, int role = Qt::DisplayRole) const override;
    bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole) override;
    QModelIndex parent(const QModelIndex & index) const  override;
    QModelIndex index(int row, int column, const QModelIndex & parent = QModelIndex()) const override;
    bool hasChildren(const QModelIndex & parent = QModelIndex()) const override;
    int rowCount(const QModelIndex & parent = QModelIndex()) const override;
    int columnCount(const QModelIndex & parent = QModelIndex()) const override;
    QVariant headerData(int section, Qt::Orientation orientation, int role) const override;
    Qt::ItemFlags flags(const QModelIndex &index) const override;
    QVector<int> availableRoles() const;
    QHash<int, QByteArray> roleNames() const override;

    bool isInitialized() const;
    bool hasData(const QModelIndex &index, int role) const;

    size_t rootCacheSize() const;
    void setRootCacheSize(size_t rootCacheSize);

Q_SIGNALS:
    void initialized();

private:
    explicit QAbstractItemModelReplica(QAbstractItemModelReplicaImplementation *rep, QtRemoteObjects::InitialAction action, const QVector<int> &rolesHint);
    QScopedPointer<QAbstractItemModelReplicaImplementation> d;
    friend class QAbstractItemModelReplicaImplementation;
    friend class QRemoteObjectNode;
};

QT_END_NAMESPACE

#endif // QREMOTEOBJECTS_ABSTRACTITEMMODELREPLICA_H
