/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl-3.0.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or (at your option) the GNU General
** Public license version 3 or any later version approved by the KDE Free
** Qt Foundation. The licenses are as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL2 and LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-2.0.html and
** https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QPAINTERVIDEOSURFACE_P_H
#define QPAINTERVIDEOSURFACE_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API. It exists purely as an
// implementation detail. This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <qtmultimediawidgetdefs.h>
#include <QtCore/qsize.h>
#include <QtGui/qimage.h>
#include <QtGui/qmatrix4x4.h>
#include <QtGui/qpaintengine.h>
#include <qabstractvideosurface.h>
#include <qvideoframe.h>

QT_USE_NAMESPACE

QT_BEGIN_NAMESPACE

class Q_AUTOTEST_EXPORT QVideoSurfacePainter
{
public:
    virtual ~QVideoSurfacePainter();

    virtual QList<QVideoFrame::PixelFormat> supportedPixelFormats(
            QAbstractVideoBuffer::HandleType handleType) const = 0;

    virtual bool isFormatSupported(const QVideoSurfaceFormat &format) const = 0;

    virtual QAbstractVideoSurface::Error start(const QVideoSurfaceFormat &format) = 0;
    virtual void stop() = 0;

    virtual QAbstractVideoSurface::Error setCurrentFrame(const QVideoFrame &frame) = 0;

    virtual QAbstractVideoSurface::Error paint(
            const QRectF &target, QPainter *painter, const QRectF &source) = 0;

    virtual void updateColors(int brightness, int contrast, int hue, int saturation) = 0;
    virtual void viewportDestroyed() {}
};

class QOpenGLContext;
class Q_AUTOTEST_EXPORT QPainterVideoSurface : public QAbstractVideoSurface
{
    Q_OBJECT
public:
    explicit QPainterVideoSurface(QObject *parent = 0);
    ~QPainterVideoSurface();

    QList<QVideoFrame::PixelFormat> supportedPixelFormats(
            QAbstractVideoBuffer::HandleType handleType = QAbstractVideoBuffer::NoHandle) const override;

    bool isFormatSupported(const QVideoSurfaceFormat &format) const override;

    bool start(const QVideoSurfaceFormat &format) override;
    void stop() override;

    bool present(const QVideoFrame &frame) override;

    int brightness() const;
    void setBrightness(int brightness);

    int contrast() const;
    void setContrast(int contrast);

    int hue() const;
    void setHue(int hue);

    int saturation() const;
    void setSaturation(int saturation);

    bool isReady() const;
    void setReady(bool ready);

    void paint(QPainter *painter, const QRectF &target, const QRectF &source = QRectF(0, 0, 1, 1));

#if QT_CONFIG(opengl)
    const QOpenGLContext *glContext() const;
    void updateGLContext();

    enum ShaderType
    {
        NoShaders = 0x00,
        FragmentProgramShader = 0x01,
        GlslShader = 0x02
    };

    Q_DECLARE_FLAGS(ShaderTypes, ShaderType)

    ShaderTypes supportedShaderTypes() const;

    ShaderType shaderType() const;
    void setShaderType(ShaderType type);
#endif

public Q_SLOTS:
    void viewportDestroyed();

Q_SIGNALS:
    void frameChanged();

private:
    void createPainter();

    QVideoSurfacePainter *m_painter;
#if QT_CONFIG(opengl)
    QOpenGLContext *m_glContext;
    ShaderTypes m_shaderTypes;
    ShaderType m_shaderType;
#endif
    int m_brightness;
    int m_contrast;
    int m_hue;
    int m_saturation;

    QVideoFrame::PixelFormat m_pixelFormat;
    QSize m_frameSize;
    QRect m_sourceRect;
    bool m_colorsDirty;
    bool m_ready;
};

#if QT_CONFIG(opengl)
Q_DECLARE_OPERATORS_FOR_FLAGS(QPainterVideoSurface::ShaderTypes)
#endif

QT_END_NAMESPACE


#endif
