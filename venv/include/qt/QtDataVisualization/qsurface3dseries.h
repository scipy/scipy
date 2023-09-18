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

#ifndef QSURFACE3DSERIES_H
#define QSURFACE3DSERIES_H

#include <QtDataVisualization/qabstract3dseries.h>
#include <QtDataVisualization/qsurfacedataproxy.h>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class QSurface3DSeriesPrivate;

class QT_DATAVISUALIZATION_EXPORT QSurface3DSeries : public QAbstract3DSeries
{
    Q_OBJECT
    Q_FLAGS(DrawFlag DrawFlags)
    Q_PROPERTY(QSurfaceDataProxy *dataProxy READ dataProxy WRITE setDataProxy NOTIFY dataProxyChanged)
    Q_PROPERTY(QPoint selectedPoint READ selectedPoint WRITE setSelectedPoint NOTIFY selectedPointChanged)
    Q_PROPERTY(bool flatShadingEnabled READ isFlatShadingEnabled WRITE setFlatShadingEnabled NOTIFY flatShadingEnabledChanged)
    Q_PROPERTY(bool flatShadingSupported READ isFlatShadingSupported NOTIFY flatShadingSupportedChanged)
    Q_PROPERTY(DrawFlags drawMode READ drawMode WRITE setDrawMode NOTIFY drawModeChanged)
    Q_PROPERTY(QImage texture READ texture WRITE setTexture NOTIFY textureChanged)
    Q_PROPERTY(QString textureFile READ textureFile WRITE setTextureFile NOTIFY textureFileChanged)

public:
    enum DrawFlag {
        DrawWireframe = 1,
        DrawSurface = 2,
        DrawSurfaceAndWireframe = DrawWireframe | DrawSurface
    };
    Q_DECLARE_FLAGS(DrawFlags, DrawFlag)

    explicit QSurface3DSeries(QObject *parent = nullptr);
    explicit QSurface3DSeries(QSurfaceDataProxy *dataProxy, QObject *parent = nullptr);
    virtual ~QSurface3DSeries();

    void setDataProxy(QSurfaceDataProxy *proxy);
    QSurfaceDataProxy *dataProxy() const;

    void setSelectedPoint(const QPoint &position);
    QPoint selectedPoint() const;
    static QPoint invalidSelectionPosition();

    void setFlatShadingEnabled(bool enabled);
    bool isFlatShadingEnabled() const;

    void setDrawMode(DrawFlags mode);
    QSurface3DSeries::DrawFlags drawMode() const;

    bool isFlatShadingSupported() const;

    void setTexture(const QImage &texture);
    QImage texture() const;
    void setTextureFile(const QString &filename);
    QString textureFile() const;

Q_SIGNALS:
    void dataProxyChanged(QSurfaceDataProxy *proxy);
    void selectedPointChanged(const QPoint &position);
    void flatShadingEnabledChanged(bool enable);
    void flatShadingSupportedChanged(bool enable);
    void drawModeChanged(QSurface3DSeries::DrawFlags mode);
    void textureChanged(const QImage &image);
    void textureFileChanged(const QString &filename);

protected:
    explicit QSurface3DSeries(QSurface3DSeriesPrivate *d, QObject *parent = nullptr);
    QSurface3DSeriesPrivate *dptr();
    const QSurface3DSeriesPrivate *dptrc() const;

private:
    Q_DISABLE_COPY(QSurface3DSeries)

    friend class Surface3DController;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
