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
#ifndef QGEOFILETILECACHE_P_H
#define QGEOFILETILECACHE_P_H

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

#include <QObject>
#include <QCache>
#include "qcache3q_p.h"
#include <QSet>
#include <QMutex>
#include <QTimer>

#include "qgeotilespec_p.h"
#include "qgeotiledmappingmanagerengine_p.h"
#include "qabstractgeotilecache_p.h"

#include <QImage>

QT_BEGIN_NAMESPACE

class QGeoMappingManager;

class QGeoTile;
class QGeoCachedTileMemory;
class QGeoFileTileCache;

class QPixmap;
class QThread;

/* This would be internal to qgeofiletilecache.cpp except that the eviction
 * policy can't be defined without it being concrete here */
class QGeoCachedTileDisk
{
public:
    ~QGeoCachedTileDisk();

    QGeoTileSpec spec;
    QString filename;
    QString format;
    QGeoFileTileCache *cache;
};

/* Custom eviction policy for the disk cache, to avoid deleting all the files
 * when the application closes */
class Q_LOCATION_PRIVATE_EXPORT QCache3QTileEvictionPolicy : public QCache3QDefaultEvictionPolicy<QGeoTileSpec,QGeoCachedTileDisk>
{
protected:
    void aboutToBeRemoved(const QGeoTileSpec &key, QSharedPointer<QGeoCachedTileDisk> obj);
    void aboutToBeEvicted(const QGeoTileSpec &key, QSharedPointer<QGeoCachedTileDisk> obj);
};

class Q_LOCATION_PRIVATE_EXPORT QGeoFileTileCache : public QAbstractGeoTileCache
{
    Q_OBJECT
public:
    QGeoFileTileCache(const QString &directory = QString(), QObject *parent = 0);
    ~QGeoFileTileCache();

    void setMaxDiskUsage(int diskUsage) override;
    int maxDiskUsage() const override;
    int diskUsage() const override;

    void setMaxMemoryUsage(int memoryUsage) override;
    int maxMemoryUsage() const override;
    int memoryUsage() const override;

    void setMinTextureUsage(int textureUsage) override;
    void setExtraTextureUsage(int textureUsage) override;
    int maxTextureUsage() const override;
    int minTextureUsage() const override;
    int textureUsage() const override;
    void clearAll() override;
    void clearMapId(const int mapId);
    void setCostStrategyDisk(CostStrategy costStrategy) override;
    CostStrategy costStrategyDisk() const override;
    void setCostStrategyMemory(CostStrategy costStrategy) override;
    CostStrategy costStrategyMemory() const override;
    void setCostStrategyTexture(CostStrategy costStrategy) override;
    CostStrategy costStrategyTexture() const override;


    QSharedPointer<QGeoTileTexture> get(const QGeoTileSpec &spec) override;

    // can be called without a specific tileCache pointer
    static void evictFromDiskCache(QGeoCachedTileDisk *td);
    static void evictFromMemoryCache(QGeoCachedTileMemory *tm);

    void insert(const QGeoTileSpec &spec,
                const QByteArray &bytes,
                const QString &format,
                QAbstractGeoTileCache::CacheAreas areas = QAbstractGeoTileCache::AllCaches) override;

    static QString tileSpecToFilenameDefault(const QGeoTileSpec &spec, const QString &format, const QString &directory);
    static QGeoTileSpec filenameToTileSpecDefault(const QString &filename);

protected:
    void init() override;
    void printStats() override;
    void loadTiles();

    QString directory() const;

    QSharedPointer<QGeoCachedTileDisk> addToDiskCache(const QGeoTileSpec &spec, const QString &filename);
    bool addToDiskCache(const QGeoTileSpec &spec, const QString &filename, const QByteArray &bytes);
    void addToMemoryCache(const QGeoTileSpec &spec, const QByteArray &bytes, const QString &format);
    QSharedPointer<QGeoTileTexture> addToTextureCache(const QGeoTileSpec &spec, const QImage &image);
    QSharedPointer<QGeoTileTexture> getFromMemory(const QGeoTileSpec &spec);
    QSharedPointer<QGeoTileTexture> getFromDisk(const QGeoTileSpec &spec);

    virtual bool isTileBogus(const QByteArray &bytes) const;
    virtual QString tileSpecToFilename(const QGeoTileSpec &spec, const QString &format, const QString &directory) const;
    virtual QGeoTileSpec filenameToTileSpec(const QString &filename) const;

    QCache3Q<QGeoTileSpec, QGeoCachedTileDisk, QCache3QTileEvictionPolicy > diskCache_;
    QCache3Q<QGeoTileSpec, QGeoCachedTileMemory > memoryCache_;
    QCache3Q<QGeoTileSpec, QGeoTileTexture > textureCache_;

    QString directory_;

    int minTextureUsage_;
    int extraTextureUsage_;
    CostStrategy costStrategyDisk_;
    CostStrategy costStrategyMemory_;
    CostStrategy costStrategyTexture_;
    bool isDiskCostSet_;
    bool isMemoryCostSet_;
    bool isTextureCostSet_;
};

QT_END_NAMESPACE

#endif // QGEOFILETILECACHE_P_H
