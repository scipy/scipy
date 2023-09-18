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

#ifndef QCUSTOM3DVOLUME_H
#define QCUSTOM3DVOLUME_H

#include <QtDataVisualization/qdatavisualizationglobal.h>
#include <QtDataVisualization/QCustom3DItem>
#include <QtGui/QColor>
#include <QtGui/QImage>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class QCustom3DVolumePrivate;

class QT_DATAVISUALIZATION_EXPORT QCustom3DVolume : public QCustom3DItem
{
    Q_OBJECT
    Q_PROPERTY(int textureWidth READ textureWidth WRITE setTextureWidth NOTIFY textureWidthChanged)
    Q_PROPERTY(int textureHeight READ textureHeight WRITE setTextureHeight NOTIFY textureHeightChanged)
    Q_PROPERTY(int textureDepth READ textureDepth WRITE setTextureDepth NOTIFY textureDepthChanged)
    Q_PROPERTY(int sliceIndexX READ sliceIndexX WRITE setSliceIndexX NOTIFY sliceIndexXChanged)
    Q_PROPERTY(int sliceIndexY READ sliceIndexY WRITE setSliceIndexY NOTIFY sliceIndexYChanged)
    Q_PROPERTY(int sliceIndexZ READ sliceIndexZ WRITE setSliceIndexZ NOTIFY sliceIndexZChanged)
    Q_PROPERTY(QVector<QRgb> colorTable READ colorTable WRITE setColorTable NOTIFY colorTableChanged)
    Q_PROPERTY(QVector<uchar> *textureData READ textureData WRITE setTextureData NOTIFY textureDataChanged)
    Q_PROPERTY(float alphaMultiplier READ alphaMultiplier WRITE setAlphaMultiplier NOTIFY alphaMultiplierChanged)
    Q_PROPERTY(bool preserveOpacity READ preserveOpacity WRITE setPreserveOpacity NOTIFY preserveOpacityChanged)
    Q_PROPERTY(bool useHighDefShader READ useHighDefShader WRITE setUseHighDefShader NOTIFY useHighDefShaderChanged)
    Q_PROPERTY(bool drawSlices READ drawSlices WRITE setDrawSlices NOTIFY drawSlicesChanged)
    Q_PROPERTY(bool drawSliceFrames READ drawSliceFrames WRITE setDrawSliceFrames NOTIFY drawSliceFramesChanged)
    Q_PROPERTY(QColor sliceFrameColor READ sliceFrameColor WRITE setSliceFrameColor NOTIFY sliceFrameColorChanged)
    Q_PROPERTY(QVector3D sliceFrameWidths READ sliceFrameWidths WRITE setSliceFrameWidths NOTIFY sliceFrameWidthsChanged)
    Q_PROPERTY(QVector3D sliceFrameGaps READ sliceFrameGaps WRITE setSliceFrameGaps NOTIFY sliceFrameGapsChanged)
    Q_PROPERTY(QVector3D sliceFrameThicknesses READ sliceFrameThicknesses WRITE setSliceFrameThicknesses NOTIFY sliceFrameThicknessesChanged)

public:

    explicit QCustom3DVolume(QObject *parent = nullptr);
    explicit QCustom3DVolume(const QVector3D &position, const QVector3D &scaling,
                             const QQuaternion &rotation, int textureWidth,
                             int textureHeight, int textureDepth,
                             QVector<uchar> *textureData, QImage::Format textureFormat,
                             const QVector<QRgb> &colorTable, QObject *parent = nullptr);
    virtual ~QCustom3DVolume();

    void setTextureWidth(int value);
    int textureWidth() const;
    void setTextureHeight(int value);
    int textureHeight() const;
    void setTextureDepth(int value);
    int textureDepth() const;
    void setTextureDimensions(int width, int height, int depth);
    int textureDataWidth() const;

    void setSliceIndexX(int value);
    int sliceIndexX() const;
    void setSliceIndexY(int value);
    int sliceIndexY() const;
    void setSliceIndexZ(int value);
    int sliceIndexZ() const;
    void setSliceIndices(int x, int y, int z);

    void setColorTable(const QVector<QRgb> &colors);
    QVector<QRgb> colorTable() const;

    void setTextureData(QVector<uchar> *data);
    QVector<uchar> *createTextureData(const QVector<QImage *> &images);
    QVector<uchar> *textureData() const;
    void setSubTextureData(Qt::Axis axis, int index, const uchar *data);
    void setSubTextureData(Qt::Axis axis, int index, const QImage &image);

    void setTextureFormat(QImage::Format format);
    QImage::Format textureFormat() const;

    void setAlphaMultiplier(float mult);
    float alphaMultiplier() const;
    void setPreserveOpacity(bool enable);
    bool preserveOpacity() const;

    void setUseHighDefShader(bool enable);
    bool useHighDefShader() const;

    void setDrawSlices(bool enable);
    bool drawSlices() const;
    void setDrawSliceFrames(bool enable);
    bool drawSliceFrames() const;

    void setSliceFrameColor(const QColor &color);
    QColor sliceFrameColor() const;
    void setSliceFrameWidths(const QVector3D &values);
    QVector3D sliceFrameWidths() const;
    void setSliceFrameGaps(const QVector3D &values);
    QVector3D sliceFrameGaps() const;
    void setSliceFrameThicknesses(const QVector3D &values);
    QVector3D sliceFrameThicknesses() const;

    QImage renderSlice(Qt::Axis axis, int index);

Q_SIGNALS:
    void textureWidthChanged(int value);
    void textureHeightChanged(int value);
    void textureDepthChanged(int value);
    void sliceIndexXChanged(int value);
    void sliceIndexYChanged(int value);
    void sliceIndexZChanged(int value);
    void colorTableChanged();
    void textureDataChanged(QVector<uchar> *data);
    void textureFormatChanged(QImage::Format format);
    void alphaMultiplierChanged(float mult);
    void preserveOpacityChanged(bool enabled);
    void useHighDefShaderChanged(bool enabled);
    void drawSlicesChanged(bool enabled);
    void drawSliceFramesChanged(bool enabled);
    void sliceFrameColorChanged(const QColor &color);
    void sliceFrameWidthsChanged(const QVector3D &values);
    void sliceFrameGapsChanged(const QVector3D &values);
    void sliceFrameThicknessesChanged(const QVector3D &values);

protected:
    QCustom3DVolumePrivate *dptr();
    const QCustom3DVolumePrivate *dptrc() const;

private:
    Q_DISABLE_COPY(QCustom3DVolume)

    friend class Abstract3DRenderer;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
