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
#ifndef QABSTRACTGEOTILECACHE_P_H
#define QABSTRACTGEOTILECACHE_P_H

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

#include <QImage>

QT_BEGIN_NAMESPACE

class QGeoMappingManager;
class QGeoMappingManagerEngine;

class QGeoTile;
class QAbstractGeoTileCache;

class QThread;

/* This is also used in the mapgeometry */
class Q_LOCATION_PRIVATE_EXPORT QGeoTileTexture
{
public:

    QGeoTileTexture();
    ~QGeoTileTexture();

    QGeoTileSpec spec;
    QImage image;
    bool textureBound;
};

class Q_LOCATION_PRIVATE_EXPORT QAbstractGeoTileCache : public QObject
{
    Q_OBJECT
public:
    enum CostStrategy {
        Unitary,
        ByteSize
    };

    enum CacheArea {
        DiskCache = 0x01,
        MemoryCache = 0x02,
        AllCaches = 0xFF
    };
    Q_DECLARE_FLAGS(CacheAreas, CacheArea)

    virtual ~QAbstractGeoTileCache();

    virtual void setMaxDiskUsage(int diskUsage);
    virtual int maxDiskUsage() const;
    virtual int diskUsage() const;

    virtual void setMaxMemoryUsage(int memoryUsage);
    virtual int maxMemoryUsage() const;
    virtual int memoryUsage() const;

    virtual void setMinTextureUsage(int textureUsage) = 0;
    virtual void setExtraTextureUsage(int textureUsage) = 0;
    virtual int maxTextureUsage() const = 0;
    virtual int minTextureUsage() const = 0;
    virtual int textureUsage() const = 0;
    virtual void clearAll() = 0;
    virtual void setCostStrategyDisk(CostStrategy costStrategy) = 0;
    virtual CostStrategy costStrategyDisk() const = 0;
    virtual void setCostStrategyMemory(CostStrategy costStrategy) = 0;
    virtual CostStrategy costStrategyMemory() const = 0;
    virtual void setCostStrategyTexture(CostStrategy costStrategy) = 0;
    virtual CostStrategy costStrategyTexture() const = 0;

    virtual QSharedPointer<QGeoTileTexture> get(const QGeoTileSpec &spec) = 0;

    virtual void insert(const QGeoTileSpec &spec,
                const QByteArray &bytes,
                const QString &format,
                QAbstractGeoTileCache::CacheAreas areas = QAbstractGeoTileCache::AllCaches) = 0;
    virtual void handleError(const QGeoTileSpec &spec, const QString &errorString);
    virtual void init() = 0;

    static QString baseCacheDirectory();
    static QString baseLocationCacheDirectory();

protected:
    QAbstractGeoTileCache(QObject *parent = 0);
    virtual void printStats() = 0;

    friend class QGeoTiledMappingManagerEngine;
};

Q_DECLARE_OPERATORS_FOR_FLAGS(QAbstractGeoTileCache::CacheAreas)

QT_END_NAMESPACE

#endif // QABSTRACTGEOTILECACHE_P_H
