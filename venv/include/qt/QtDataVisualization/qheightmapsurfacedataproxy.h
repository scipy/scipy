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

#ifndef QHEIGHTMAPSURFACEDATAPROXY_H
#define QHEIGHTMAPSURFACEDATAPROXY_H

#include <QtDataVisualization/qsurfacedataproxy.h>
#include <QtGui/QImage>
#include <QtCore/QString>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class QHeightMapSurfaceDataProxyPrivate;

class QT_DATAVISUALIZATION_EXPORT QHeightMapSurfaceDataProxy : public QSurfaceDataProxy
{
    Q_OBJECT

    Q_PROPERTY(QImage heightMap READ heightMap WRITE setHeightMap NOTIFY heightMapChanged)
    Q_PROPERTY(QString heightMapFile READ heightMapFile WRITE setHeightMapFile NOTIFY heightMapFileChanged)
    Q_PROPERTY(float minXValue READ minXValue WRITE setMinXValue NOTIFY minXValueChanged)
    Q_PROPERTY(float maxXValue READ maxXValue WRITE setMaxXValue NOTIFY maxXValueChanged)
    Q_PROPERTY(float minZValue READ minZValue WRITE setMinZValue NOTIFY minZValueChanged)
    Q_PROPERTY(float maxZValue READ maxZValue WRITE setMaxZValue NOTIFY maxZValueChanged)

public:
    explicit QHeightMapSurfaceDataProxy(QObject *parent = nullptr);
    explicit QHeightMapSurfaceDataProxy(const QImage &image, QObject *parent = nullptr);
    explicit QHeightMapSurfaceDataProxy(const QString &filename, QObject *parent = nullptr);
    virtual ~QHeightMapSurfaceDataProxy();

    void setHeightMap(const QImage &image);
    QImage heightMap() const;
    void setHeightMapFile(const QString &filename);
    QString heightMapFile() const;

    void setValueRanges(float minX, float maxX, float minZ, float maxZ);
    void setMinXValue(float min);
    float minXValue() const;
    void setMaxXValue(float max);
    float maxXValue() const;
    void setMinZValue(float min);
    float minZValue() const;
    void setMaxZValue(float max);
    float maxZValue() const;

Q_SIGNALS:
    void heightMapChanged(const QImage &image);
    void heightMapFileChanged(const QString &filename);
    void minXValueChanged(float value);
    void maxXValueChanged(float value);
    void minZValueChanged(float value);
    void maxZValueChanged(float value);

protected:
    explicit QHeightMapSurfaceDataProxy(QHeightMapSurfaceDataProxyPrivate *d, QObject *parent = nullptr);
    QHeightMapSurfaceDataProxyPrivate *dptr();
    const QHeightMapSurfaceDataProxyPrivate *dptrc() const;

private:
    Q_DISABLE_COPY(QHeightMapSurfaceDataProxy)

    friend class Surface3DController;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
