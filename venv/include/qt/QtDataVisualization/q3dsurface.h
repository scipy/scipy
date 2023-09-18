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

#ifndef Q3DSURFACE_H
#define Q3DSURFACE_H

#include <QtDataVisualization/qabstract3dgraph.h>
#include <QtDataVisualization/qvalue3daxis.h>
#include <QtDataVisualization/qsurface3dseries.h>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class Q3DSurfacePrivate;

class QT_DATAVISUALIZATION_EXPORT Q3DSurface : public QAbstract3DGraph
{
    Q_OBJECT
    Q_PROPERTY(QValue3DAxis *axisX READ axisX WRITE setAxisX NOTIFY axisXChanged)
    Q_PROPERTY(QValue3DAxis *axisY READ axisY WRITE setAxisY NOTIFY axisYChanged)
    Q_PROPERTY(QValue3DAxis *axisZ READ axisZ WRITE setAxisZ NOTIFY axisZChanged)
    Q_PROPERTY(QSurface3DSeries *selectedSeries READ selectedSeries NOTIFY selectedSeriesChanged)
    Q_PROPERTY(bool flipHorizontalGrid READ flipHorizontalGrid WRITE setFlipHorizontalGrid NOTIFY flipHorizontalGridChanged)

public:
    explicit Q3DSurface(const QSurfaceFormat *format = nullptr, QWindow *parent = nullptr);
    virtual ~Q3DSurface();

    void addSeries(QSurface3DSeries *series);
    void removeSeries(QSurface3DSeries *series);
    QList<QSurface3DSeries *> seriesList() const;

    // Axes
    void setAxisX(QValue3DAxis *axis);
    QValue3DAxis *axisX() const;
    void setAxisY(QValue3DAxis *axis);
    QValue3DAxis *axisY() const;
    void setAxisZ(QValue3DAxis *axis);
    QValue3DAxis *axisZ() const;
    void addAxis(QValue3DAxis *axis);
    void releaseAxis(QValue3DAxis *axis);
    QList<QValue3DAxis *> axes() const;

    QSurface3DSeries *selectedSeries() const;
    void setFlipHorizontalGrid(bool flip);
    bool flipHorizontalGrid() const;

Q_SIGNALS:
    void axisXChanged(QValue3DAxis *axis);
    void axisYChanged(QValue3DAxis *axis);
    void axisZChanged(QValue3DAxis *axis);
    void selectedSeriesChanged(QSurface3DSeries *series);
    void flipHorizontalGridChanged(bool flip);

private:
    Q3DSurfacePrivate *dptr();
    const Q3DSurfacePrivate *dptrc() const;
    Q_DISABLE_COPY(Q3DSurface)
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
