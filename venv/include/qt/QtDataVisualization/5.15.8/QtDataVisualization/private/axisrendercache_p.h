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

#ifndef AXISRENDERCACHE_P_H
#define AXISRENDERCACHE_P_H

#include "datavisualizationglobal_p.h"
#include "drawer_p.h"
#include <QtCore/QPointer>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class AxisRenderCache : public QObject
{
    Q_OBJECT
public:
    AxisRenderCache();
    virtual ~AxisRenderCache();

    void setDrawer(Drawer *drawer);

    void setType(QAbstract3DAxis::AxisType type);
    inline QAbstract3DAxis::AxisType type() const { return m_type; }
    void setTitle(const QString &title);
    inline const QString &title() const { return m_title; }
    void setLabels(const QStringList &labels);
    inline const QStringList &labels() const { return m_labels; }
    inline void setMin(float min) { m_min = min; }
    inline float min() const { return m_min; }
    inline void setMax(float max) { m_max = max; }
    inline float max() const { return m_max; }
    inline void setSegmentCount(int count) { m_segmentCount = count; m_positionsDirty = true; }
    inline int segmentCount() const { return m_segmentCount; }
    inline void setSubSegmentCount(int count) { m_subSegmentCount = count; m_positionsDirty = true; }
    inline int subSegmentCount() const { return m_subSegmentCount; }
    inline void setLabelFormat(const QString &format) { m_labelFormat = format; }
    inline const QString &labelFormat() const { return m_labelFormat; }
    inline void setReversed(bool enable) { m_reversed = enable; m_positionsDirty = true; }
    inline bool reversed() const { return m_reversed; }
    inline void setFormatter(QValue3DAxisFormatter *formatter)
    {
        m_formatter = formatter; m_positionsDirty = true;
    }
    inline QValue3DAxisFormatter *formatter() const { return m_formatter; }
    inline void setCtrlFormatter(QValue3DAxisFormatter *formatter)
    {
        m_ctrlFormatter = formatter;
    }
    inline QValue3DAxisFormatter *ctrlFormatter() const { return m_ctrlFormatter; }

    inline LabelItem &titleItem() { return m_titleItem; }
    inline QList<LabelItem *> &labelItems() { return m_labelItems; }
    inline float gridLinePosition(int index) { return m_adjustedGridLinePositions.at(index); }
    inline int gridLineCount() { return m_adjustedGridLinePositions.size(); }
    inline float labelPosition(int index) { return m_adjustedLabelPositions.at(index); }
    inline int labelCount() {
        // Some value axis formatters may opt to not show all labels,
        // so use positions array for determining count in that case.
        if (m_type == QAbstract3DAxis::AxisTypeValue)
            return m_adjustedLabelPositions.size();
        else
            return m_labels.size();
    }
    void updateAllPositions();
    inline bool positionsDirty() const { return m_positionsDirty; }
    inline void markPositionsDirty() { m_positionsDirty = true; }
    inline void setTranslate(float translate) { m_translate = translate; m_positionsDirty = true; }
    inline float translate() { return m_translate; }
    inline void setScale(float scale) { m_scale = scale; m_positionsDirty = true; }
    inline float scale() { return m_scale; }
    inline float positionAt(float value)
    {
        if (m_reversed)
            return (1.0f - m_formatter->positionAt(value)) * m_scale + m_translate;
        else
            return m_formatter->positionAt(value) * m_scale + m_translate;
    }
    inline float labelAutoRotation() const { return m_labelAutoRotation; }
    inline void setLabelAutoRotation(float angle) { m_labelAutoRotation = angle; }
    inline bool isTitleVisible() const { return m_titleVisible; }
    inline void setTitleVisible(bool visible) { m_titleVisible = visible; }
    inline bool isTitleFixed() const { return m_titleFixed; }
    inline void setTitleFixed(bool fixed) { m_titleFixed = fixed; }

    void updateTextures();
    void clearLabels();

private:
    int maxLabelWidth(const QStringList &labels) const;

    // Cached axis values
    QAbstract3DAxis::AxisType m_type;
    QString m_title;
    QStringList m_labels;
    float m_min;
    float m_max;
    int m_segmentCount;
    int m_subSegmentCount;
    QString m_labelFormat;
    bool m_reversed;
    QFont m_font;
    QValue3DAxisFormatter *m_formatter;
    QPointer<QValue3DAxisFormatter> m_ctrlFormatter;

    // Renderer items
    Drawer *m_drawer; // Not owned
    LabelItem m_titleItem;
    QList<LabelItem *> m_labelItems;
    QVector<float> m_adjustedGridLinePositions;
    QVector<float> m_adjustedLabelPositions;
    bool m_positionsDirty;
    float m_translate;
    float m_scale;
    float m_labelAutoRotation;
    bool m_titleVisible;
    bool m_titleFixed;

    Q_DISABLE_COPY(AxisRenderCache)
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
