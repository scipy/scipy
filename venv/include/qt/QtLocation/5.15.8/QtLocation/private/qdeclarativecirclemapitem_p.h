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

#ifndef QDECLARATIVECIRCLEMAPITEM_H
#define QDECLARATIVECIRCLEMAPITEM_H

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
#include <QtLocation/private/qdeclarativegeomapitembase_p.h>
#include <QtLocation/private/qdeclarativepolylinemapitem_p_p.h>
#include <QSGGeometryNode>
#include <QSGFlatColorMaterial>
#include <QtPositioning/QGeoCircle>

QT_BEGIN_NAMESPACE

class QDeclarativeCircleMapItemPrivate;
class Q_LOCATION_PRIVATE_EXPORT QDeclarativeCircleMapItem : public QDeclarativeGeoMapItemBase
{
    Q_OBJECT
    Q_ENUMS(Backend)

    Q_PROPERTY(QGeoCoordinate center READ center WRITE setCenter NOTIFY centerChanged)
    Q_PROPERTY(qreal radius READ radius WRITE setRadius NOTIFY radiusChanged)
    Q_PROPERTY(QColor color READ color WRITE setColor NOTIFY colorChanged)
    Q_PROPERTY(QDeclarativeMapLineProperties *border READ border CONSTANT)
    Q_PROPERTY(Backend backend READ backend WRITE setBackend NOTIFY backendChanged REVISION 15)

public:
    enum Backend {
        Software = 0,
        OpenGL = 1
    };

    explicit QDeclarativeCircleMapItem(QQuickItem *parent = nullptr);
    ~QDeclarativeCircleMapItem() override;

    virtual void setMap(QDeclarativeGeoMap *quickMap, QGeoMap *map) override;
    virtual QSGNode *updateMapItemPaintNode(QSGNode *, UpdatePaintNodeData *) override;

    QGeoCoordinate center();
    void setCenter(const QGeoCoordinate &center);

    qreal radius() const;
    void setRadius(qreal radius);

    QColor color() const;
    void setColor(const QColor &color);

    QDeclarativeMapLineProperties *border();

    bool contains(const QPointF &point) const override;
    const QGeoShape &geoShape() const override;
    void setGeoShape(const QGeoShape &shape) override;

    Backend backend() const;
    void setBackend(Backend b);

Q_SIGNALS:
    void centerChanged(const QGeoCoordinate &center);
    void radiusChanged(qreal radius);
    void colorChanged(const QColor &color);
    void backendChanged();

protected:
    void geometryChanged(const QRectF &newGeometry, const QRectF &oldGeometry) override;
    void updatePolish() override;
    void possiblySwitchBackend(const QGeoCoordinate &oldCenter, qreal oldRadius, const QGeoCoordinate &newCenter, qreal newRadius);

protected Q_SLOTS:
    void markSourceDirtyAndUpdate();
    void onLinePropertiesChanged();
    virtual void afterViewportChanged(const QGeoMapViewportChangeEvent &event) override;

private:
    void updateCirclePath();

private:
    QGeoCircle m_circle;
    QDeclarativeMapLineProperties m_border;
    QColor m_color;
    bool m_dirtyMaterial;
    bool m_updatingGeometry;
    Backend m_backend = Software;

    QScopedPointer<QDeclarativeCircleMapItemPrivate> m_d;

    friend class QDeclarativeCircleMapItemPrivate;
    friend class QDeclarativeCircleMapItemPrivateCPU;
    friend class QDeclarativeCircleMapItemPrivateOpenGL;
};

//////////////////////////////////////////////////////////////////////

QT_END_NAMESPACE

QML_DECLARE_TYPE(QDeclarativeCircleMapItem)

#endif /* QDECLARATIVECIRCLEMAPITEM_H */
