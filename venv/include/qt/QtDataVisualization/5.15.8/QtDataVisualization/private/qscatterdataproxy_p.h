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

#ifndef QSCATTERDATAPROXY_P_H
#define QSCATTERDATAPROXY_P_H

#include "qscatterdataproxy.h"
#include "qabstractdataproxy_p.h"
#include "qscatterdataitem.h"

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class QAbstract3DAxis;

class QScatterDataProxyPrivate : public QAbstractDataProxyPrivate
{
    Q_OBJECT
public:
    QScatterDataProxyPrivate(QScatterDataProxy *q);
    virtual ~QScatterDataProxyPrivate();

    void resetArray(QScatterDataArray *newArray);
    void setItem(int index, const QScatterDataItem &item);
    void setItems(int index, const QScatterDataArray &items);
    int addItem(const QScatterDataItem &item);
    int addItems(const QScatterDataArray &items);
    void insertItem(int index, const QScatterDataItem &item);
    void insertItems(int index, const QScatterDataArray &items);
    void removeItems(int index, int removeCount);
    void limitValues(QVector3D &minValues, QVector3D &maxValues, QAbstract3DAxis *axisX,
                     QAbstract3DAxis *axisY, QAbstract3DAxis *axisZ) const;
    bool isValidValue(float axisValue, float value, QAbstract3DAxis *axis) const;

    virtual void setSeries(QAbstract3DSeries *series);
private:
    QScatterDataProxy *qptr();
    QScatterDataArray *m_dataArray;

    friend class QScatterDataProxy;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
