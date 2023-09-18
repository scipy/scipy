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

#ifndef QGEOMAPITEMGEOMETRY_H
#define QGEOMAPITEMGEOMETRY_H

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

#include <QPainterPath>
#include <QPointF>
#include <QRectF>
#include <QVector>
#include <QGeoCoordinate>
#include <QVector2D>
#include <QList>

QT_BEGIN_NAMESPACE

class QSGGeometry;
class QGeoMap;

class Q_LOCATION_PRIVATE_EXPORT QGeoMapItemGeometry
{
public:
    QGeoMapItemGeometry();
    virtual ~QGeoMapItemGeometry();

    inline bool isSourceDirty() const { return sourceDirty_; }
    inline bool isScreenDirty() const { return screenDirty_; }
    inline void markSourceDirty() { sourceDirty_ = true; screenDirty_ = true; }
    inline void markScreenDirty() { screenDirty_ = true; clipToViewport_ = true; }
    inline void markFullScreenDirty() { screenDirty_ = true; clipToViewport_ = false;}
    inline void markClean() { screenDirty_ = (sourceDirty_ = false); clipToViewport_ = true;}
    inline void clearScreen() { screenDirty_ = false; }

    inline void setPreserveGeometry(bool value, const QGeoCoordinate &geoLeftBound = QGeoCoordinate())
    {
        preserveGeometry_ = value;
        if (preserveGeometry_)
            geoLeftBound_ = geoLeftBound;
    }
    inline QGeoCoordinate geoLeftBound() { return geoLeftBound_; }

    inline QRectF sourceBoundingBox() const { return sourceBounds_; }
    inline QRectF screenBoundingBox() const { return screenBounds_; }
    inline void clearBounds() { sourceBounds_ = screenBounds_ = QRectF(); firstPointOffset_ = QPointF(); }

    inline QPointF firstPointOffset() const { return firstPointOffset_; }
    void translate(const QPointF &offset);

    inline const QGeoCoordinate &origin() const { return srcOrigin_; }

    QPainterPath screenOutline() const {
        return screenOutline_;
    }

    virtual bool contains(const QPointF &screenPoint) const {
        return screenOutline_.contains(screenPoint);
    }

    inline QVector2D vertex(quint32 index) const {
        return QVector2D(screenVertices_[index]);
    }

    inline QVector<QPointF> vertices() const { return screenVertices_; }
    inline QVector<quint32> indices() const { return screenIndices_; }

    inline bool isIndexed() const { return (!screenIndices_.isEmpty()); }

    /* Size is # of triangles */
    inline quint32 size() const
    {
        if (isIndexed())
            return screenIndices_.size() / 3;
        else
            return screenVertices_.size() / 3;
    }

    inline void clear() { firstPointOffset_ = QPointF(0,0);
                          screenVertices_.clear(); screenIndices_.clear(); }

    void allocateAndFill(QSGGeometry *geom) const;

    double geoDistanceToScreenWidth(const QGeoMap &map,
                                           const QGeoCoordinate &fromCoord,
                                           const QGeoCoordinate &toCoord);

    static QRectF translateToCommonOrigin(const QList<QGeoMapItemGeometry *> &geoms);

    mutable bool m_dataChanged = false;

private:
    QGeoMapItemGeometry(const QGeoMapItemGeometry &other); // Or else it may crash on copy
    QGeoMapItemGeometry &operator= (const QGeoMapItemGeometry & other); // Or else it may crash on copy

protected:
    bool sourceDirty_;
    bool screenDirty_;
    bool clipToViewport_;
    bool preserveGeometry_;
    QGeoCoordinate geoLeftBound_;

    QPointF firstPointOffset_;

    QPainterPath screenOutline_;

    QRectF sourceBounds_;
    QRectF screenBounds_;

    QGeoCoordinate srcOrigin_;

    QVector<QPointF> screenVertices_;
    QVector<quint32> screenIndices_;
};

QT_END_NAMESPACE

#endif // QGEOMAPITEMGEOMETRY_H
