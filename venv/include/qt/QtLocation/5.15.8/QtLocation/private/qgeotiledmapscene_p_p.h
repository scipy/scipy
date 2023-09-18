/****************************************************************************
**
** Copyright (C) 2018 The Qt Company Ltd.
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
#ifndef QGEOTILEDMAPSCENE_P_P_H
#define QGEOTILEDMAPSCENE_P_P_H

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

#include "qgeotiledmapscene_p.h"
#include <QtCore/private/qobject_p.h>
#include <QtPositioning/private/qdoublevector3d_p.h>
#include <QtQuick/QSGImageNode>
#include <QtQuick/private/qsgdefaultimagenode_p.h>
#include <QtQuick/QQuickWindow>
#include "qgeocameradata_p.h"
#include "qgeotilespec_p.h"

QT_BEGIN_NAMESPACE

class Q_LOCATION_PRIVATE_EXPORT QGeoTiledMapTileContainerNode : public QSGTransformNode
{
public:
    void addChild(const QGeoTileSpec &spec, QSGImageNode *node)
    {
        tiles.insert(spec, node);
        appendChildNode(node);
    }
    QHash<QGeoTileSpec, QSGImageNode *> tiles;
};

class Q_LOCATION_PRIVATE_EXPORT QGeoTiledMapRootNode : public QSGClipNode
{
public:
    QGeoTiledMapRootNode()
        : isTextureLinear(false)
        , geometry(QSGGeometry::defaultAttributes_Point2D(), 4)
        , root(new QSGTransformNode())
        , tiles(new QGeoTiledMapTileContainerNode())
        , wrapLeft(new QGeoTiledMapTileContainerNode())
        , wrapRight(new QGeoTiledMapTileContainerNode())
    {
        setIsRectangular(true);
        setGeometry(&geometry);
        root->appendChildNode(tiles);
        root->appendChildNode(wrapLeft);
        root->appendChildNode(wrapRight);
        appendChildNode(root);
    }

    ~QGeoTiledMapRootNode()
    {
        qDeleteAll(textures);
    }

    void setClipRect(const QRect &rect)
    {
        if (rect != clipRect) {
            QSGGeometry::updateRectGeometry(&geometry, rect);
            QSGClipNode::setClipRect(rect);
            clipRect = rect;
            markDirty(DirtyGeometry);
        }
    }

    void updateTiles(QGeoTiledMapTileContainerNode *root,
                     QGeoTiledMapScenePrivate *d,
                     double camAdjust,
                     QQuickWindow *window,
                     bool ogl);

    bool isTextureLinear;

    QSGGeometry geometry;
    QRect clipRect;

    QSGTransformNode *root;

    QGeoTiledMapTileContainerNode *tiles;        // The majority of the tiles
    QGeoTiledMapTileContainerNode *wrapLeft;     // When zoomed out, the tiles that wrap around on the left.
    QGeoTiledMapTileContainerNode *wrapRight;    // When zoomed out, the tiles that wrap around on the right

    QHash<QGeoTileSpec, QSGTexture *> textures;

#ifdef QT_LOCATION_DEBUG
    double m_sideLengthPixel;
    QMap<double, QList<QGeoTileSpec>> m_droppedTiles;
#endif
};

class Q_LOCATION_PRIVATE_EXPORT QGeoTiledMapScenePrivate : public QObjectPrivate
{
    Q_DECLARE_PUBLIC(QGeoTiledMapScene)
public:
    QGeoTiledMapScenePrivate();
    ~QGeoTiledMapScenePrivate();

    void addTile(const QGeoTileSpec &spec, QSharedPointer<QGeoTileTexture> texture);

    void setVisibleTiles(const QSet<QGeoTileSpec> &visibleTiles);
    void removeTiles(const QSet<QGeoTileSpec> &oldTiles);
    bool buildGeometry(const QGeoTileSpec &spec, QSGImageNode *imageNode, bool &overzooming);
    void updateTileBounds(const QSet<QGeoTileSpec> &tiles);
    void setupCamera();
    inline bool isTiltedOrRotated() { return (m_cameraData.tilt() > 0.0) || (m_cameraData.bearing() > 0.0); }

public:

    QSize m_screenSize; // in pixels
    int m_tileSize; // the pixel resolution for each tile
    QGeoCameraData m_cameraData;
    QRectF m_visibleArea;
    QSet<QGeoTileSpec> m_visibleTiles;

    QDoubleVector3D m_cameraUp;
    QDoubleVector3D m_cameraEye;
    QDoubleVector3D m_cameraCenter;
    QMatrix4x4 m_projectionMatrix;

    // scales up the tile geometry and the camera altitude, resulting in no visible effect
    // other than to control the accuracy of the render by keeping the values in a sensible range
    double m_scaleFactor;

    // rounded down, positive zoom is zooming in, corresponding to reduced altitude
    int m_intZoomLevel;

    // mercatorToGrid transform
    // the number of tiles in each direction for the whole map (earth) at the current zoom level.
    // it is 1<<zoomLevel
    int m_sideLength;
    double m_mapEdgeSize;

    QHash<QGeoTileSpec, QSharedPointer<QGeoTileTexture> > m_textures;
    QVector<QGeoTileSpec> m_updatedTextures;

    // tilesToGrid transform
    int m_minTileX; // the minimum tile index, i.e. 0 to sideLength which is 1<< zoomLevel
    int m_minTileY;
    int m_maxTileX;
    int m_maxTileY;
    int m_tileXWrapsBelow; // the wrap point as a tile index
    bool m_linearScaling;
    bool m_dropTextures;

#ifdef QT_LOCATION_DEBUG
    double m_sideLengthPixel;
    QGeoTiledMapRootNode *m_mapRoot = nullptr;
#endif
};

QT_END_NAMESPACE

#endif // QGEOTILEDMAPSCENE_P_P_H
