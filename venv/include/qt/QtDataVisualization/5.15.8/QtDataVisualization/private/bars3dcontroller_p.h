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

#ifndef Q3DBARSCONTROLLER_p_H
#define Q3DBARSCONTROLLER_p_H

#include "datavisualizationglobal_p.h"
#include "abstract3dcontroller_p.h"

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class Bars3DRenderer;
class QBar3DSeries;

struct Bars3DChangeBitField {
    bool multiSeriesScalingChanged  : 1;
    bool barSpecsChanged            : 1;
    bool selectedBarChanged         : 1;
    bool rowsChanged                : 1;
    bool itemChanged                : 1;
    bool floorLevelChanged          : 1;

    Bars3DChangeBitField() :
        multiSeriesScalingChanged(true),
        barSpecsChanged(true),
        selectedBarChanged(true),
        rowsChanged(false),
        itemChanged(false),
        floorLevelChanged(false)
    {
    }
};

class QT_DATAVISUALIZATION_EXPORT Bars3DController : public Abstract3DController
{
    Q_OBJECT

public:
    struct ChangeItem {
        QBar3DSeries *series;
        QPoint point;
    };
    struct ChangeRow {
        QBar3DSeries *series;
        int row;
    };

private:
    Bars3DChangeBitField m_changeTracker;
    QVector<ChangeItem> m_changedItems;
    QVector<ChangeRow> m_changedRows;

    // Interaction
    QPoint m_selectedBar;     // Points to row & column in data window.
    QBar3DSeries *m_selectedBarSeries; // Points to the series for which the bar is selected in
                                       // single series selection cases.
    QBar3DSeries *m_primarySeries; // Category axis labels are taken from the primary series

    // Look'n'feel
    bool m_isMultiSeriesUniform;
    bool m_isBarSpecRelative;
    GLfloat m_barThicknessRatio;
    QSizeF m_barSpacing;
    float m_floorLevel;

    // Rendering
    Bars3DRenderer *m_renderer;

public:
    explicit Bars3DController(QRect rect, Q3DScene *scene = 0);
    ~Bars3DController();

    virtual void initializeOpenGL();
    virtual void synchDataToRenderer();

    void setMultiSeriesScaling(bool uniform);
    bool multiSeriesScaling() const;

    // bar thickness, spacing between bars, and is spacing relative to thickness or absolute
    // y -component sets the thickness/spacing of z -direction
    // With relative 0.0f means side-to-side, 1.0f = one thickness in between
    void setBarSpecs(GLfloat thicknessRatio = 1.0f,
                     const QSizeF &spacing = QSizeF(1.0, 1.0),
                     bool relative = true);
    GLfloat barThickness();
    QSizeF barSpacing();
    bool isBarSpecRelative();
    void setFloorLevel(float level);
    float floorLevel() const;

    inline QBar3DSeries *selectedSeries() const { return m_selectedBarSeries; }

    void setSelectionMode(QAbstract3DGraph::SelectionFlags mode);
    void setSelectedBar(const QPoint &position, QBar3DSeries *series, bool enterSlice);
    virtual void clearSelection();

    virtual void handleAxisAutoAdjustRangeChangedInOrientation(
            QAbstract3DAxis::AxisOrientation orientation, bool autoAdjust);
    virtual void handleSeriesVisibilityChangedBySender(QObject *sender);
    virtual void handlePendingClick();

    static QPoint invalidSelectionPosition();

    virtual void setAxisX(QAbstract3DAxis *axis);
    virtual void setAxisZ(QAbstract3DAxis *axis);

    virtual void setPrimarySeries(QBar3DSeries *series);
    virtual QBar3DSeries *primarySeries() const;
    virtual void addSeries(QAbstract3DSeries *series);
    virtual void removeSeries(QAbstract3DSeries *series);
    virtual void insertSeries(int index, QAbstract3DSeries *series);
    virtual QList<QBar3DSeries *> barSeriesList();

    virtual void handleAxisRangeChangedBySender(QObject *sender);
    virtual void adjustAxisRanges();

public Q_SLOTS:
    void handleArrayReset();
    void handleRowsAdded(int startIndex, int count);
    void handleRowsChanged(int startIndex, int count);
    void handleRowsRemoved(int startIndex, int count);
    void handleRowsInserted(int startIndex, int count);
    void handleItemChanged(int rowIndex, int columnIndex);
    void handleDataRowLabelsChanged();
    void handleDataColumnLabelsChanged();

Q_SIGNALS:
    void primarySeriesChanged(QBar3DSeries *series);
    void selectedSeriesChanged(QBar3DSeries *series);

protected:
    virtual QAbstract3DAxis *createDefaultAxis(QAbstract3DAxis::AxisOrientation orientation);

private:
    void adjustSelectionPosition(QPoint &pos, const QBar3DSeries *series);

    Q_DISABLE_COPY(Bars3DController)
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
