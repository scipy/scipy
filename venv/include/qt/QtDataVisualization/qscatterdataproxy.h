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

#ifndef QSCATTERDATAPROXY_H
#define QSCATTERDATAPROXY_H

#include <QtDataVisualization/qabstractdataproxy.h>
#include <QtDataVisualization/qscatterdataitem.h>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class QScatterDataProxyPrivate;
class QScatter3DSeries;

typedef QVector<QScatterDataItem> QScatterDataArray;

class QT_DATAVISUALIZATION_EXPORT QScatterDataProxy : public QAbstractDataProxy
{
    Q_OBJECT

    Q_PROPERTY(int itemCount READ itemCount NOTIFY itemCountChanged)
    Q_PROPERTY(QScatter3DSeries *series READ series NOTIFY seriesChanged)

public:
    explicit QScatterDataProxy(QObject *parent = nullptr);
    virtual ~QScatterDataProxy();

    QScatter3DSeries *series() const;
    int itemCount() const;
    const QScatterDataArray *array() const;
    const QScatterDataItem *itemAt(int index) const;

    void resetArray(QScatterDataArray *newArray);

    void setItem(int index, const QScatterDataItem &item);
    void setItems(int index, const QScatterDataArray &items);

    int addItem(const QScatterDataItem &item);
    int addItems(const QScatterDataArray &items);

    void insertItem(int index, const QScatterDataItem &item);
    void insertItems(int index, const QScatterDataArray &items);

    void removeItems(int index, int removeCount);

Q_SIGNALS:
    void arrayReset();
    void itemsAdded(int startIndex, int count);
    void itemsChanged(int startIndex, int count);
    void itemsRemoved(int startIndex, int count);
    void itemsInserted(int startIndex, int count);

    void itemCountChanged(int count);
    void seriesChanged(QScatter3DSeries *series);

protected:
    explicit QScatterDataProxy(QScatterDataProxyPrivate *d, QObject *parent = nullptr);
    QScatterDataProxyPrivate *dptr();
    const QScatterDataProxyPrivate *dptrc() const;

private:
    Q_DISABLE_COPY(QScatterDataProxy)

    friend class Scatter3DController;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
