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

#ifndef QSCATTER3DSERIES_H
#define QSCATTER3DSERIES_H

#include <QtDataVisualization/qabstract3dseries.h>
#include <QtDataVisualization/qscatterdataproxy.h>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class QScatter3DSeriesPrivate;

class QT_DATAVISUALIZATION_EXPORT QScatter3DSeries : public QAbstract3DSeries
{
    Q_OBJECT
    Q_PROPERTY(QScatterDataProxy *dataProxy READ dataProxy WRITE setDataProxy NOTIFY dataProxyChanged)
    Q_PROPERTY(int selectedItem READ selectedItem WRITE setSelectedItem NOTIFY selectedItemChanged)
    Q_PROPERTY(float itemSize READ itemSize WRITE setItemSize NOTIFY itemSizeChanged)

public:
    explicit QScatter3DSeries(QObject *parent = nullptr);
    explicit QScatter3DSeries(QScatterDataProxy *dataProxy, QObject *parent = nullptr);
    virtual ~QScatter3DSeries();

    void setDataProxy(QScatterDataProxy *proxy);
    QScatterDataProxy *dataProxy() const;

    void setSelectedItem(int index);
    int selectedItem() const;
    static int invalidSelectionIndex();

    void setItemSize(float size);
    float itemSize() const;

Q_SIGNALS:
    void dataProxyChanged(QScatterDataProxy *proxy);
    void selectedItemChanged(int index);
    void itemSizeChanged(float size);

protected:
    explicit QScatter3DSeries(QScatter3DSeriesPrivate *d, QObject *parent = nullptr);
    QScatter3DSeriesPrivate *dptr();
    const QScatter3DSeriesPrivate *dptrc() const;

private:
    Q_DISABLE_COPY(QScatter3DSeries)

    friend class Scatter3DController;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
