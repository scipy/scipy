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

#ifndef QGSTREAMERPLAYERSESSION_P_H
#define QGSTREAMERPLAYERSESSION_P_H

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

#include <QtMultimedia/private/qtmultimediaglobal_p.h>
#include <QObject>
#include <QtCore/qmutex.h>
#include <QtNetwork/qnetworkrequest.h>
#include <private/qgstreamerplayercontrol_p.h>
#include <private/qgstreamerbushelper_p.h>
#include <qmediaplayer.h>
#include <qmediastreamscontrol.h>
#include <qaudioformat.h>

#if QT_CONFIG(gstreamer_app)
#include <private/qgstappsrc_p.h>
#endif

#include <gst/gst.h>

QT_BEGIN_NAMESPACE

class QGstreamerBusHelper;
class QGstreamerMessage;

class QGstreamerVideoRendererInterface;
class QGstreamerVideoProbeControl;
class QGstreamerAudioProbeControl;

typedef enum {
  GST_AUTOPLUG_SELECT_TRY,
  GST_AUTOPLUG_SELECT_EXPOSE,
  GST_AUTOPLUG_SELECT_SKIP
} GstAutoplugSelectResult;

class Q_GSTTOOLS_EXPORT QGstreamerPlayerSession
    : public QObject
    , public QGstreamerBusMessageFilter
{
Q_OBJECT
Q_INTERFACES(QGstreamerBusMessageFilter)

public:
    QGstreamerPlayerSession(QObject *parent);
    virtual ~QGstreamerPlayerSession();

    GstElement *playbin() const;
    GstElement *pipeline() const { return m_pipeline; }
    QGstreamerBusHelper *bus() const { return m_busHelper; }

    QNetworkRequest request() const;

    QMediaPlayer::State state() const { return m_state; }
    QMediaPlayer::State pendingState() const { return m_pendingState; }

    qint64 duration() const;
    qint64 position() const;

    int volume() const;
    bool isMuted() const;

    bool isAudioAvailable() const;

    void setVideoRenderer(QObject *renderer);
    QGstreamerVideoRendererInterface *renderer() const { return m_renderer; }
    bool isVideoAvailable() const;

    bool isSeekable() const;

    qreal playbackRate() const;
    void setPlaybackRate(qreal rate);

    QMediaTimeRange availablePlaybackRanges() const;

    QMap<QByteArray ,QVariant> tags() const { return m_tags; }
    QMap<QString,QVariant> streamProperties(int streamNumber) const { return m_streamProperties[streamNumber]; }
    int streamCount() const { return m_streamProperties.count(); }
    QMediaStreamsControl::StreamType streamType(int streamNumber) { return m_streamTypes.value(streamNumber, QMediaStreamsControl::UnknownStream); }

    int activeStream(QMediaStreamsControl::StreamType streamType) const;
    void setActiveStream(QMediaStreamsControl::StreamType streamType, int streamNumber);

    bool processBusMessage(const QGstreamerMessage &message) override;

#if QT_CONFIG(gstreamer_app)
    QGstAppSrc *appsrc() const { return m_appSrc; }
    static void configureAppSrcElement(GObject*, GObject*, GParamSpec*,QGstreamerPlayerSession* _this);
#endif

    bool isLiveSource() const;

    void addProbe(QGstreamerVideoProbeControl* probe);
    void removeProbe(QGstreamerVideoProbeControl* probe);

    void addProbe(QGstreamerAudioProbeControl* probe);
    void removeProbe(QGstreamerAudioProbeControl* probe);

    void endOfMediaReset();

public slots:
    void loadFromUri(const QNetworkRequest &url);
    void loadFromStream(const QNetworkRequest &url, QIODevice *stream);
    bool play();
    bool pause();
    void stop();

    bool seek(qint64 pos);

    void setVolume(int volume);
    void setMuted(bool muted);

    void showPrerollFrames(bool enabled);

signals:
    void durationChanged(qint64 duration);
    void positionChanged(qint64 position);
    void stateChanged(QMediaPlayer::State state);
    void volumeChanged(int volume);
    void mutedStateChanged(bool muted);
    void audioAvailableChanged(bool audioAvailable);
    void videoAvailableChanged(bool videoAvailable);
    void bufferingProgressChanged(int percentFilled);
    void playbackFinished();
    void tagsChanged();
    void streamsChanged();
    void seekableChanged(bool);
    void error(int error, const QString &errorString);
    void invalidMedia();
    void playbackRateChanged(qreal);
    void rendererChanged();
    void pipelineChanged();

private slots:
    void getStreamsInfo();
    void setSeekable(bool);
    void finishVideoOutputChange();
    void updateVideoRenderer();
    void updateVideoResolutionTag();
    void updateVolume();
    void updateMuted();
    void updateDuration();

private:
    static void playbinNotifySource(GObject *o, GParamSpec *p, gpointer d);
    static void handleVolumeChange(GObject *o, GParamSpec *p, gpointer d);
    static void handleMutedChange(GObject *o, GParamSpec *p, gpointer d);
#if !GST_CHECK_VERSION(1,0,0)
    static void insertColorSpaceElement(GstElement *element, gpointer data);
#endif
    static void handleElementAdded(GstBin *bin, GstElement *element, QGstreamerPlayerSession *session);
    static void handleStreamsChange(GstBin *bin, gpointer user_data);
    static GstAutoplugSelectResult handleAutoplugSelect(GstBin *bin, GstPad *pad, GstCaps *caps, GstElementFactory *factory, QGstreamerPlayerSession *session);

    void processInvalidMedia(QMediaPlayer::Error errorCode, const QString& errorString);

    void removeVideoBufferProbe();
    void addVideoBufferProbe();
    void removeAudioBufferProbe();
    void addAudioBufferProbe();
    void flushVideoProbes();
    void resumeVideoProbes();
    bool parsePipeline();
    bool setPipeline(GstElement *pipeline);
    void resetElements();
    void initPlaybin();
    void setBus(GstBus *bus);

    QNetworkRequest m_request;
    QMediaPlayer::State m_state = QMediaPlayer::StoppedState;
    QMediaPlayer::State m_pendingState = QMediaPlayer::StoppedState;
    QGstreamerBusHelper *m_busHelper = nullptr;
    GstElement *m_playbin = nullptr;
    GstElement *m_pipeline = nullptr;

    GstElement *m_videoSink = nullptr;

    GstElement *m_videoOutputBin = nullptr;
    GstElement *m_videoIdentity = nullptr;
#if !GST_CHECK_VERSION(1,0,0)
    GstElement *m_colorSpace = nullptr;
    bool m_usingColorspaceElement = false;
#endif
    GstElement *m_pendingVideoSink = nullptr;
    GstElement *m_nullVideoSink = nullptr;

    GstElement *m_audioSink = nullptr;
    GstElement *m_volumeElement = nullptr;

    GstBus *m_bus = nullptr;
    QObject *m_videoOutput = nullptr;
    QGstreamerVideoRendererInterface *m_renderer = nullptr;

#if QT_CONFIG(gstreamer_app)
    QGstAppSrc *m_appSrc = nullptr;
#endif

    QMap<QByteArray, QVariant> m_tags;
    QList< QMap<QString,QVariant> > m_streamProperties;
    QList<QMediaStreamsControl::StreamType> m_streamTypes;
    QMap<QMediaStreamsControl::StreamType, int> m_playbin2StreamOffset;

    QGstreamerVideoProbeControl *m_videoProbe = nullptr;
    QGstreamerAudioProbeControl *m_audioProbe = nullptr;

    int m_volume = 100;
    qreal m_playbackRate = 1.0;
    bool m_muted = false;
    bool m_audioAvailable = false;
    bool m_videoAvailable = false;
    bool m_seekable = false;

    mutable qint64 m_lastPosition = 0;
    qint64 m_duration = 0;
    int m_durationQueries = 0;

    bool m_displayPrerolledFrame = true;

    enum SourceType
    {
        UnknownSrc,
        SoupHTTPSrc,
        UDPSrc,
        MMSSrc,
        RTSPSrc,
    };
    SourceType m_sourceType = UnknownSrc;
    bool m_everPlayed = false;
    bool m_isLiveSource = false;

    gulong pad_probe_id = 0;
};

QT_END_NAMESPACE

#endif // QGSTREAMERPLAYERSESSION_H
