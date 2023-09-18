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

#ifndef Q3DSCATTERCONTROLLER_p_H
#define Q3DSCATTERCONTROLLER_p_H

#include "datavisualizationglobal_p.h"
#include "abstract3dcontroller_p.h"

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class Scatter3DRenderer;
class QScatterDataProxy;
class QScatter3DSeries;

struct Scatter3DChangeBitField {
    bool selectedItemChanged : 1;
    bool itemChanged         : 1;

    Scatter3DChangeBitField() :
        selectedItemChanged(true),
        itemChanged(false)
    {
    }
};

class QT_DATAVISUALIZATION_EXPORT Scatter3DController : public Abstract3DController
{
    Q_OBJECT

public:
    struct ChangeItem {
        QScatter3DSeries *series;
        int index;
    };
private:
    Scatter3DChangeBitField m_changeTracker;
    QVector<ChangeItem> m_changedItems;

    // Rendering
    Scatter3DRenderer *m_renderer;
    int m_selectedItem;
    QScatter3DSeries *m_selectedItemSeries; // Points to the series for which the bar is selected
                                            // in single series selection cases.

    struct InsertRemoveRecord {
        bool m_isInsert;
        int m_startIndex;
        int m_count;
        QAbstract3DSeries *m_series;

        InsertRemoveRecord() :
            m_isInsert(false),
            m_startIndex(0),
            m_count(0),
            m_series(0)
        {}

        InsertRemoveRecord(bool isInsert, int startIndex, int count, QAbstract3DSeries *series) :
            m_isInsert(isInsert),
            m_startIndex(startIndex),
            m_count(count),
            m_series(series)
        {}
    };

    QVector<InsertRemoveRecord> m_insertRemoveRecords;
    bool m_recordInsertsAndRemoves;

public:
    explicit Scatter3DController(QRect rect, Q3DScene *scene = 0);
    ~Scatter3DController();

    virtual void initializeOpenGL();

    // Change selection mode
    void setSelectionMode(QAbstract3DGraph::SelectionFlags mode);

    inline QScatter3DSeries *selectedSeries() const { return m_selectedItemSeries; }

    void setSelectedItem(int index, QScatter3DSeries *series);
    static inline int invalidSelectionIndex() { return -1; }
    virtual void clearSelection();

    void synchDataToRenderer();

    virtual void addSeries(QAbstract3DSeries *series);
    virtual void removeSeries(QAbstract3DSeries *series);
    virtual QList<QScatter3DSeries *> scatterSeriesList();

    virtual void handleAxisAutoAdjustRangeChangedInOrientation(
            QAbstract3DAxis::AxisOrientation orientation, bool autoAdjust);
    virtual void handleAxisRangeChangedBySender(QObject *sender);
    virtual void handlePendingClick();
    virtual void adjustAxisRanges();

public Q_SLOTS:
    void handleArrayReset();
    void handleItemsAdded(int startIndex, int count);
    void handleItemsChanged(int startIndex, int count);
    void handleItemsRemoved(int startIndex, int count);
    void handleItemsInserted(int startIndex, int count);

Q_SIGNALS:
    void selectedSeriesChanged(QScatter3DSeries *series);

protected:
    virtual void startRecordingRemovesAndInserts();

private:

    Q_DISABLE_COPY(Scatter3DController)
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
