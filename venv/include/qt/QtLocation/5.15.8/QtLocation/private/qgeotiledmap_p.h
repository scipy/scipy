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
#ifndef QGEOTILEDMAP_P_H
#define QGEOTILEDMAP_P_H

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

#include <QObject>
#include <QString>
#include <QtLocation/private/qlocationglobal_p.h>
#include <QtLocation/private/qgeomap_p.h>
#include <QtLocation/private/qgeocameradata_p.h>
#include <QtLocation/private/qgeomaptype_p.h>

#include <QtPositioning/private/qdoublevector2d_p.h>

QT_BEGIN_NAMESPACE

class QGeoTileSpec;
class QGeoTileTexture;
class QAbstractGeoTileCache;
class QGeoTiledMapPrivate;
class QGeoTiledMappingManagerEngine;
class QGeoTileRequestManager;

class QQuickWindow;
class QSGNode;

class QPointF;

class Q_LOCATION_PRIVATE_EXPORT QGeoTiledMap : public QGeoMap
{
    Q_OBJECT
    Q_DECLARE_PRIVATE(QGeoTiledMap)
public:
    enum PrefetchStyle { NoPrefetching, PrefetchNeighbourLayer, PrefetchTwoNeighbourLayers };
    QGeoTiledMap(QGeoTiledMappingManagerEngine *engine, QObject *parent);
    virtual ~QGeoTiledMap();

    QAbstractGeoTileCache *tileCache();
    QGeoTileRequestManager *requestManager();
    void updateTile(const QGeoTileSpec &spec);
    void setPrefetchStyle(PrefetchStyle style);

    void prefetchData() override;
    void clearData() override;
    Capabilities capabilities() const override;

    void setCopyrightVisible(bool visible) override;

public Q_SLOTS:
    virtual void clearScene(int mapId);

protected:
    QSGNode *updateSceneGraph(QSGNode *, QQuickWindow *window) override;
    virtual void evaluateCopyrights(const QSet<QGeoTileSpec> &visibleTiles);

    QGeoTiledMap(QGeoTiledMapPrivate &dd, QGeoTiledMappingManagerEngine *engine, QObject *parent);

private Q_SLOTS:
    void handleTileVersionChanged();

private:
    Q_DISABLE_COPY(QGeoTiledMap)

};

QT_END_NAMESPACE

#endif // QGEOMAP_P_H
