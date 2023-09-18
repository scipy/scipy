/****************************************************************************
**
** Copyright (C) 2015 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the QtLocation module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL3$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see http://www.qt.io/terms-conditions. For further
** information use the contact form at http://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPLv3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or later as published by the Free
** Software Foundation and appearing in the file LICENSE.GPL included in
** the packaging of this file. Please review the following information to
** ensure the GNU General Public License version 2.0 requirements will be
** met: http://www.gnu.org/licenses/gpl-2.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QDECLARATIVEGEOMAPQUICKITEM_H
#define QDECLARATIVEGEOMAPQUICKITEM_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <QtLocation/private/qlocationglobal_p.h>

#include <QtQuick/QQuickItem>
#include <QtQuick/QSGNode>

#include <QtLocation/private/qdeclarativegeomap_p.h>
#include <QtLocation/private/qdeclarativegeomapitembase_p.h>
#include <QtPositioning/qgeoshape.h>

QT_BEGIN_NAMESPACE

class QMapQuickItemMatrix4x4 : public QQuickTransform
{
public:
    QMapQuickItemMatrix4x4(QObject *parent = nullptr);

    void setMatrix(const QMatrix4x4& matrix);
    void applyTo(QMatrix4x4 *matrix) const override;

    QMatrix4x4 m_matrix;
};

class Q_LOCATION_PRIVATE_EXPORT QDeclarativeGeoMapQuickItem : public QDeclarativeGeoMapItemBase
{
    Q_OBJECT
    Q_PROPERTY(QGeoCoordinate coordinate READ coordinate WRITE setCoordinate NOTIFY coordinateChanged)
    Q_PROPERTY(QPointF anchorPoint READ anchorPoint WRITE setAnchorPoint NOTIFY anchorPointChanged)
    Q_PROPERTY(qreal zoomLevel READ zoomLevel WRITE setZoomLevel NOTIFY zoomLevelChanged)
    Q_PROPERTY(QQuickItem *sourceItem READ sourceItem WRITE setSourceItem NOTIFY sourceItemChanged)

public:
    explicit QDeclarativeGeoMapQuickItem(QQuickItem *parent = 0);
    ~QDeclarativeGeoMapQuickItem();

    virtual void setMap(QDeclarativeGeoMap *quickMap, QGeoMap *map) override;

    void setCoordinate(const QGeoCoordinate &coordinate);
    QGeoCoordinate coordinate();

    void setSourceItem(QQuickItem *sourceItem);
    QQuickItem *sourceItem();

    void setAnchorPoint(const QPointF &anchorPoint);
    QPointF anchorPoint() const;

    void setZoomLevel(qreal zoomLevel);
    qreal zoomLevel() const;

    const QGeoShape &geoShape() const override;
    void setGeoShape(const QGeoShape &shape) override;

Q_SIGNALS:
    void coordinateChanged();
    void sourceItemChanged();
    void anchorPointChanged();
    void zoomLevelChanged();

protected:
    void geometryChanged(const QRectF &newGeometry, const QRectF &oldGeometry) override;
    void updatePolish() override;
    bool childMouseEventFilter(QQuickItem *item, QEvent *event) override;

protected Q_SLOTS:
    virtual void afterChildrenChanged() override;
    virtual void afterViewportChanged(const QGeoMapViewportChangeEvent &event) override;

private:
    qreal scaleFactor();
    QGeoCoordinate dragStartCoordinate_;
    QGeoCoordinate coordinate_;
    QGeoRectangle geoshape_;
    QPointer<QQuickItem> sourceItem_;
    QQuickItem *opacityContainer_;
    QPointF anchorPoint_;
    qreal zoomLevel_;
    bool mapAndSourceItemSet_;
    bool updatingGeometry_;
    QMapQuickItemMatrix4x4 *matrix_;

    friend class QDeclarativeGeoMap;
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QDeclarativeGeoMapQuickItem)

#endif
