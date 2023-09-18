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

#ifndef SURFACE3DCONTROLLER_P_H
#define SURFACE3DCONTROLLER_P_H

#include "abstract3dcontroller_p.h"
#include "datavisualizationglobal_p.h"

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class Surface3DRenderer;
class QSurface3DSeries;

struct Surface3DChangeBitField {
    bool selectedPointChanged      : 1;
    bool rowsChanged               : 1;
    bool itemChanged               : 1;
    bool flipHorizontalGridChanged : 1;
    bool surfaceTextureChanged     : 1;

    Surface3DChangeBitField() :
        selectedPointChanged(true),
        rowsChanged(false),
        itemChanged(false),
        flipHorizontalGridChanged(true),
        surfaceTextureChanged(true)
    {
    }
};

class QT_DATAVISUALIZATION_EXPORT Surface3DController : public Abstract3DController
{
    Q_OBJECT

public:
    struct ChangeItem {
        QSurface3DSeries *series;
        QPoint point;
    };
    struct ChangeRow {
        QSurface3DSeries *series;
        int row;
    };

private:
    Surface3DChangeBitField m_changeTracker;
    Surface3DRenderer *m_renderer;
    QPoint m_selectedPoint;
    QSurface3DSeries *m_selectedSeries; // Points to the series for which the point is selected in
                                        // single series selection cases.
    bool m_flatShadingSupported;
    QVector<ChangeItem> m_changedItems;
    QVector<ChangeRow> m_changedRows;
    bool m_flipHorizontalGrid;
    QVector<QSurface3DSeries *> m_changedTextures;

public:
    explicit Surface3DController(QRect rect, Q3DScene *scene = 0);
    ~Surface3DController();

    virtual void initializeOpenGL();
    virtual void synchDataToRenderer();

    void setSelectionMode(QAbstract3DGraph::SelectionFlags mode);
    void setSelectedPoint(const QPoint &position, QSurface3DSeries *series, bool enterSlice);
    virtual void clearSelection();

    inline QSurface3DSeries *selectedSeries() const { return m_selectedSeries; }

    virtual void handleAxisAutoAdjustRangeChangedInOrientation(
            QAbstract3DAxis::AxisOrientation orientation, bool autoAdjust);
    virtual void handleAxisRangeChangedBySender(QObject *sender);
    virtual void handleSeriesVisibilityChangedBySender(QObject *sender);
    virtual void handlePendingClick();
    virtual void adjustAxisRanges();

    static QPoint invalidSelectionPosition();
    bool isFlatShadingSupported();

    virtual void addSeries(QAbstract3DSeries *series);
    virtual void removeSeries(QAbstract3DSeries *series);
    virtual QList<QSurface3DSeries *> surfaceSeriesList();

    void setFlipHorizontalGrid(bool flip);
    bool flipHorizontalGrid() const;

    void updateSurfaceTexture(QSurface3DSeries *series);

public Q_SLOTS:
    void handleArrayReset();
    void handleRowsAdded(int startIndex, int count);
    void handleRowsChanged(int startIndex, int count);
    void handleRowsRemoved(int startIndex, int count);
    void handleRowsInserted(int startIndex, int count);
    void handleItemChanged(int rowIndex, int columnIndex);

    void handleFlatShadingSupportedChange(bool supported);

Q_SIGNALS:
    void selectedSeriesChanged(QSurface3DSeries *series);
    void flipHorizontalGridChanged(bool flip);

private:
    Q_DISABLE_COPY(Surface3DController)
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
