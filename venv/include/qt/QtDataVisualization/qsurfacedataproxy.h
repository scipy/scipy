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

#ifndef QSURFACEDATAPROXY_H
#define QSURFACEDATAPROXY_H

#include <QtDataVisualization/qabstractdataproxy.h>
#include <QtDataVisualization/qsurfacedataitem.h>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class QSurfaceDataProxyPrivate;
class QSurface3DSeries;

typedef QVector<QSurfaceDataItem> QSurfaceDataRow;
typedef QList<QSurfaceDataRow *> QSurfaceDataArray;

class QT_DATAVISUALIZATION_EXPORT QSurfaceDataProxy : public QAbstractDataProxy
{
    Q_OBJECT

    Q_PROPERTY(int rowCount READ rowCount NOTIFY rowCountChanged)
    Q_PROPERTY(int columnCount READ columnCount NOTIFY columnCountChanged)
    Q_PROPERTY(QSurface3DSeries *series READ series NOTIFY seriesChanged)

public:
    explicit QSurfaceDataProxy(QObject *parent = nullptr);
    virtual ~QSurfaceDataProxy();

    QSurface3DSeries *series() const;
    int rowCount() const;
    int columnCount() const;
    const QSurfaceDataArray *array() const;
    const QSurfaceDataItem *itemAt(int rowIndex, int columnIndex) const;
    const QSurfaceDataItem *itemAt(const QPoint &position) const;

    void resetArray(QSurfaceDataArray *newArray);

    void setRow(int rowIndex, QSurfaceDataRow *row);
    void setRows(int rowIndex, const QSurfaceDataArray &rows);

    void setItem(int rowIndex, int columnIndex, const QSurfaceDataItem &item);
    void setItem(const QPoint &position, const QSurfaceDataItem &item);

    int addRow(QSurfaceDataRow *row);
    int addRows(const QSurfaceDataArray &rows);

    void insertRow(int rowIndex, QSurfaceDataRow *row);
    void insertRows(int rowIndex, const QSurfaceDataArray &rows);

    void removeRows(int rowIndex, int removeCount);

Q_SIGNALS:
    void arrayReset();
    void rowsAdded(int startIndex, int count);
    void rowsChanged(int startIndex, int count);
    void rowsRemoved(int startIndex, int count);
    void rowsInserted(int startIndex, int count);
    void itemChanged(int rowIndex, int columnIndex);

    void rowCountChanged(int count);
    void columnCountChanged(int count);
    void seriesChanged(QSurface3DSeries *series);

protected:
    explicit QSurfaceDataProxy(QSurfaceDataProxyPrivate *d, QObject *parent = nullptr);
    QSurfaceDataProxyPrivate *dptr();
    const QSurfaceDataProxyPrivate *dptrc() const;

private:
    Q_DISABLE_COPY(QSurfaceDataProxy)

    friend class Surface3DController;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
