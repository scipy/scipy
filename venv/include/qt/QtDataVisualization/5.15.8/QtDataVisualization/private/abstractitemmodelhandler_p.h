/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Data Visualization module of the Qt Toolkit.
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

//
//  W A R N I N G
//  -------------
//
// This file is not part of the QtDataVisualization API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#ifndef ABSTRACTITEMMODELHANDLER_P_H
#define ABSTRACTITEMMODELHANDLER_P_H

#include "datavisualizationglobal_p.h"
#include <QtCore/QAbstractItemModel>
#include <QtCore/QPointer>
#include <QtCore/QTimer>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class AbstractItemModelHandler : public QObject
{
    Q_OBJECT
public:
    AbstractItemModelHandler(QObject *parent = 0);
    virtual ~AbstractItemModelHandler();

    virtual void setItemModel(QAbstractItemModel *itemModel);
    virtual QAbstractItemModel *itemModel() const;

public Q_SLOTS:
    virtual void handleColumnsInserted(const QModelIndex &parent, int start, int end);
    virtual void handleColumnsMoved(const QModelIndex &sourceParent, int sourceStart,
                                    int sourceEnd, const QModelIndex &destinationParent,
                                    int destinationColumn);
    virtual void handleColumnsRemoved(const QModelIndex &parent, int start, int end);
    virtual void handleDataChanged(const QModelIndex &topLeft, const QModelIndex &bottomRight,
                                   const QVector<int> &roles = QVector<int> ());
    virtual void handleLayoutChanged(const QList<QPersistentModelIndex> &parents = QList<QPersistentModelIndex>(),
                                     QAbstractItemModel::LayoutChangeHint hint = QAbstractItemModel::NoLayoutChangeHint);
    virtual void handleModelReset();
    virtual void handleRowsInserted(const QModelIndex &parent, int start, int end);
    virtual void handleRowsMoved(const QModelIndex &sourceParent, int sourceStart, int sourceEnd,
                                 const QModelIndex &destinationParent, int destinationRow);
    virtual void handleRowsRemoved(const QModelIndex &parent, int start, int end);

    virtual void handleMappingChanged();
    virtual void handlePendingResolve();

Q_SIGNALS:
    void itemModelChanged(const QAbstractItemModel *itemModel);

protected:
    virtual void resolveModel() = 0;

    QPointer<QAbstractItemModel> m_itemModel;  // Not owned
    bool resolvePending;
    QTimer m_resolveTimer;
    bool m_fullReset;

private:
    Q_DISABLE_COPY(AbstractItemModelHandler)
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
