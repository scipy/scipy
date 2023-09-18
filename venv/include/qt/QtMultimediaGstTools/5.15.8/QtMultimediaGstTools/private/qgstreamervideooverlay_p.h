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

#ifndef QGSTREAMERVIDEOOVERLAY_P_H
#define QGSTREAMERVIDEOOVERLAY_P_H

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

#include <private/qgstreamerbushelper_p.h>
#include <private/qgstreamerbufferprobe_p.h>
#include <QtGui/qwindowdefs.h>
#include <QtCore/qsize.h>

QT_BEGIN_NAMESPACE

class QGstreamerSinkProperties;
class QGstreamerVideoOverlay
        : public QObject
        , public QGstreamerSyncMessageFilter
        , public QGstreamerBusMessageFilter
        , private QGstreamerBufferProbe
{
    Q_OBJECT
    Q_INTERFACES(QGstreamerSyncMessageFilter QGstreamerBusMessageFilter)
public:
    explicit QGstreamerVideoOverlay(QObject *parent = 0, const QByteArray &elementName = QByteArray());
    virtual ~QGstreamerVideoOverlay();

    GstElement *videoSink() const;
    void setVideoSink(GstElement *);
    QSize nativeVideoSize() const;

    void setWindowHandle(WId id);
    void expose();
    void setRenderRectangle(const QRect &rect);

    bool isActive() const;

    Qt::AspectRatioMode aspectRatioMode() const;
    void setAspectRatioMode(Qt::AspectRatioMode mode);

    int brightness() const;
    void setBrightness(int brightness);

    int contrast() const;
    void setContrast(int contrast);

    int hue() const;
    void setHue(int hue);

    int saturation() const;
    void setSaturation(int saturation);

    bool processSyncMessage(const QGstreamerMessage &message) override;
    bool processBusMessage(const QGstreamerMessage &message) override;

Q_SIGNALS:
    void nativeVideoSizeChanged();
    void activeChanged();
    void brightnessChanged(int brightness);
    void contrastChanged(int contrast);
    void hueChanged(int hue);
    void saturationChanged(int saturation);

private:
    void setWindowHandle_helper(WId id);
    void updateIsActive();
    void probeCaps(GstCaps *caps) override;
    static void showPrerollFrameChanged(GObject *, GParamSpec *, QGstreamerVideoOverlay *);

    GstElement *m_videoSink = nullptr;
    QSize m_nativeVideoSize;
    bool m_isActive = false;

    QGstreamerSinkProperties *m_sinkProperties = nullptr;
    WId m_windowId = 0;
};

QT_END_NAMESPACE

#endif // QGSTREAMERVIDEOOVERLAY_P_H

