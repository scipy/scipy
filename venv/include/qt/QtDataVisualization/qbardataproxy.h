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

#ifndef QBARDATAPROXY_H
#define QBARDATAPROXY_H

#include <QtDataVisualization/qabstractdataproxy.h>
#include <QtDataVisualization/qbardataitem.h>
#include <QtCore/QVector>
#include <QtCore/QStringList>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class QBarDataProxyPrivate;
class QBar3DSeries;

typedef QVector<QBarDataItem> QBarDataRow;
typedef QList<QBarDataRow *> QBarDataArray;

class QT_DATAVISUALIZATION_EXPORT QBarDataProxy : public QAbstractDataProxy
{
    Q_OBJECT

    Q_PROPERTY(int rowCount READ rowCount NOTIFY rowCountChanged)
    Q_PROPERTY(QStringList rowLabels READ rowLabels WRITE setRowLabels NOTIFY rowLabelsChanged)
    Q_PROPERTY(QStringList columnLabels READ columnLabels WRITE setColumnLabels NOTIFY columnLabelsChanged)
    Q_PROPERTY(QBar3DSeries *series READ series NOTIFY seriesChanged)
public:
    explicit QBarDataProxy(QObject *parent = nullptr);
    virtual ~QBarDataProxy();

    QBar3DSeries *series() const;
    int rowCount() const;

    QStringList rowLabels() const;
    void setRowLabels(const QStringList &labels);
    QStringList columnLabels() const;
    void setColumnLabels(const QStringList &labels);

    const QBarDataArray *array() const;
    const QBarDataRow *rowAt(int rowIndex) const;
    const QBarDataItem *itemAt(int rowIndex, int columnIndex) const;
    const QBarDataItem *itemAt(const QPoint &position) const;

    void resetArray();
    void resetArray(QBarDataArray *newArray);
    void resetArray(QBarDataArray *newArray, const QStringList &rowLabels,
                    const QStringList &columnLabels);

    void setRow(int rowIndex, QBarDataRow *row);
    void setRow(int rowIndex, QBarDataRow *row, const QString &label);
    void setRows(int rowIndex, const QBarDataArray &rows);
    void setRows(int rowIndex, const QBarDataArray &rows, const QStringList &labels);

    void setItem(int rowIndex, int columnIndex, const QBarDataItem &item);
    void setItem(const QPoint &position, const QBarDataItem &item);

    int addRow(QBarDataRow *row);
    int addRow(QBarDataRow *row, const QString &label);
    int addRows(const QBarDataArray &rows);
    int addRows(const QBarDataArray &rows, const QStringList &labels);

    void insertRow(int rowIndex, QBarDataRow *row);
    void insertRow(int rowIndex, QBarDataRow *row, const QString &label);
    void insertRows(int rowIndex, const QBarDataArray &rows);
    void insertRows(int rowIndex, const QBarDataArray &rows, const QStringList &labels);

    void removeRows(int rowIndex, int removeCount, bool removeLabels = true);

Q_SIGNALS:
    void arrayReset();
    void rowsAdded(int startIndex, int count);
    void rowsChanged(int startIndex, int count);
    void rowsRemoved(int startIndex, int count);
    void rowsInserted(int startIndex, int count);
    void itemChanged(int rowIndex, int columnIndex);

    void rowCountChanged(int count);
    void rowLabelsChanged();
    void columnLabelsChanged();
    void seriesChanged(QBar3DSeries *series);

protected:
    explicit QBarDataProxy(QBarDataProxyPrivate *d, QObject *parent = nullptr);
    QBarDataProxyPrivate *dptr();
    const QBarDataProxyPrivate *dptrc() const;

private:
    Q_DISABLE_COPY(QBarDataProxy)

    friend class Bars3DController;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
