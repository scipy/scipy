/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Charts module of the Qt Toolkit.
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

//  W A R N I N G
//  -------------
//
// This file is not part of the Qt Chart API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#ifndef LEGENDMARKERITEM_P_H
#define LEGENDMARKERITEM_P_H

#include <QtCharts/QChartGlobal>
#include <QtCharts/QLegend>
#include <QGraphicsObject>
#include <QtGui/QFont>
#include <QtGui/QBrush>
#include <QtGui/QPen>
#include <QtWidgets/QGraphicsTextItem>
#include <QtWidgets/QGraphicsLayoutItem>
#include <QtCharts/private/qchartglobal_p.h>

QT_CHARTS_BEGIN_NAMESPACE

class QLegendMarkerPrivate;

class Q_CHARTS_PRIVATE_EXPORT LegendMarkerItem : public QGraphicsObject, public QGraphicsLayoutItem
{
    Q_OBJECT
    Q_INTERFACES(QGraphicsLayoutItem)
public:
    enum ItemType {
        TypeRect,
        TypeLine,
        TypeCircle
    };

    explicit LegendMarkerItem(QLegendMarkerPrivate *marker, QGraphicsObject *parent = nullptr);
    ~LegendMarkerItem();

    void setPen(const QPen &pen);
    QPen pen() const;

    void setBrush(const QBrush &brush);
    QBrush brush() const;

    void setSeriesPen(const QPen &pen);
    void setSeriesBrush(const QBrush &brush);

    void setFont(const QFont &font);
    QFont font() const;

    void setLabel(const QString label);
    QString label() const;

    void setLabelBrush(const QBrush &brush);
    QBrush labelBrush() const;

    void setGeometry(const QRectF &rect);
    QRectF boundingRect() const;
    QRectF markerRect() const;

    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,QWidget *widget = nullptr);
    QSizeF sizeHint (Qt::SizeHint which, const QSizeF &constraint) const;

    void hoverEnterEvent(QGraphicsSceneHoverEvent *event);
    void hoverLeaveEvent(QGraphicsSceneHoverEvent *event);

    QString displayedLabel() const;
    void setToolTip(const QString &tooltip);

    QLegend::MarkerShape markerShape() const;
    void setMarkerShape(QLegend::MarkerShape shape);

    void updateMarkerShapeAndSize();
    QLegend::MarkerShape effectiveMarkerShape() const;
    qreal effectiveMarkerWidth() const;

    ItemType itemType() const { return m_itemType; }

Q_SIGNALS:
    void markerRectChanged();

protected:
    void setItemBrushAndPen();
    void setItemRect();
    bool useMaxWidth() const;

    QLegendMarkerPrivate *m_marker; // Knows
    QRectF m_defaultMarkerRect;
    QRectF m_markerRect;
    QRectF m_boundingRect;
    QGraphicsTextItem *m_textItem;
    QGraphicsItem *m_markerItem;
    qreal m_margin;
    qreal m_space;
    QString m_label;
    QLegend::MarkerShape m_markerShape;

    QBrush m_labelBrush;
    QPen m_pen;
    QBrush m_brush;
    QPen m_seriesPen;
    QBrush m_seriesBrush;
    QFont m_font;
    bool m_hovering;

    ItemType m_itemType;

    friend class QLegendMarker;
    friend class QLegendMarkerPrivate;
    friend class LegendLayout;
};

QT_CHARTS_END_NAMESPACE

#endif // LEGENDMARKERITEM_P_H
