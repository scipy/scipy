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

#ifndef QMEDIAPLAYER_H
#define QMEDIAPLAYER_H

#include <QtMultimedia/qtmultimediaglobal.h>
#include <QtMultimedia/qmediaobject.h>
#include <QtMultimedia/qmediacontent.h>
#include <QtMultimedia/qmediaenumdebug.h>
#include <QtMultimedia/qaudio.h>

#include <QtNetwork/qnetworkconfiguration.h>

QT_BEGIN_NAMESPACE


class QAbstractVideoSurface;
class QMediaPlaylist;
class QVideoWidget;
class QGraphicsVideoItem;

class QMediaPlayerPrivate;
class Q_MULTIMEDIA_EXPORT QMediaPlayer : public QMediaObject
{
    Q_OBJECT
    Q_PROPERTY(QMediaContent media READ media WRITE setMedia NOTIFY mediaChanged)
    Q_PROPERTY(QMediaContent currentMedia READ currentMedia NOTIFY currentMediaChanged)
    Q_PROPERTY(QMediaPlaylist* playlist READ playlist WRITE setPlaylist)
    Q_PROPERTY(qint64 duration READ duration NOTIFY durationChanged)
    Q_PROPERTY(qint64 position READ position WRITE setPosition NOTIFY positionChanged)
    Q_PROPERTY(int volume READ volume WRITE setVolume NOTIFY volumeChanged)
    Q_PROPERTY(bool muted READ isMuted WRITE setMuted NOTIFY mutedChanged)
    Q_PROPERTY(int bufferStatus READ bufferStatus NOTIFY bufferStatusChanged)
    Q_PROPERTY(bool audioAvailable READ isAudioAvailable NOTIFY audioAvailableChanged)
    Q_PROPERTY(bool videoAvailable READ isVideoAvailable NOTIFY videoAvailableChanged)
    Q_PROPERTY(bool seekable READ isSeekable NOTIFY seekableChanged)
    Q_PROPERTY(qreal playbackRate READ playbackRate WRITE setPlaybackRate NOTIFY playbackRateChanged)
    Q_PROPERTY(State state READ state NOTIFY stateChanged)
    Q_PROPERTY(MediaStatus mediaStatus READ mediaStatus NOTIFY mediaStatusChanged)
    Q_PROPERTY(QAudio::Role audioRole READ audioRole WRITE setAudioRole NOTIFY audioRoleChanged)
    Q_PROPERTY(QString customAudioRole READ customAudioRole WRITE setCustomAudioRole NOTIFY customAudioRoleChanged)
    Q_PROPERTY(QString error READ errorString)
    Q_ENUMS(State)
    Q_ENUMS(MediaStatus)
    Q_ENUMS(Error)

public:
    enum State
    {
        StoppedState,
        PlayingState,
        PausedState
    };

    enum MediaStatus
    {
        UnknownMediaStatus,
        NoMedia,
        LoadingMedia,
        LoadedMedia,
        StalledMedia,
        BufferingMedia,
        BufferedMedia,
        EndOfMedia,
        InvalidMedia
    };

    enum Flag
    {
        LowLatency = 0x01,
        StreamPlayback = 0x02,
        VideoSurface = 0x04
    };
    Q_DECLARE_FLAGS(Flags, Flag)

    enum Error
    {
        NoError,
        ResourceError,
        FormatError,
        NetworkError,
        AccessDeniedError,
        ServiceMissingError,
        MediaIsPlaylist
    };

    explicit QMediaPlayer(QObject *parent = nullptr, Flags flags = Flags());
    ~QMediaPlayer();

    static QMultimedia::SupportEstimate hasSupport(const QString &mimeType,
                                            const QStringList& codecs = QStringList(),
                                                   Flags flags = Flags());
    static QStringList supportedMimeTypes(Flags flags = Flags());

    void setVideoOutput(QVideoWidget *);
    void setVideoOutput(QGraphicsVideoItem *);
    void setVideoOutput(QAbstractVideoSurface *surface);
    void setVideoOutput(const QVector<QAbstractVideoSurface *> &surfaces);

    QMediaContent media() const;
    const QIODevice *mediaStream() const;
    QMediaPlaylist *playlist() const;
    QMediaContent currentMedia() const;

    State state() const;
    MediaStatus mediaStatus() const;

    qint64 duration() const;
    qint64 position() const;

    int volume() const;
    bool isMuted() const;
    bool isAudioAvailable() const;
    bool isVideoAvailable() const;

    int bufferStatus() const;

    bool isSeekable() const;
    qreal playbackRate() const;

    Error error() const;
    QString errorString() const;

#ifndef QT_NO_BEARERMANAGEMENT
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
    QT_DEPRECATED_VERSION_5_15 QNetworkConfiguration currentNetworkConfiguration() const;
QT_WARNING_POP
#endif

    QMultimedia::AvailabilityStatus availability() const override;

    QAudio::Role audioRole() const;
    void setAudioRole(QAudio::Role audioRole);
    QList<QAudio::Role> supportedAudioRoles() const;
    QString customAudioRole() const;
    void setCustomAudioRole(const QString &audioRole);
    QStringList supportedCustomAudioRoles() const;

public Q_SLOTS:
    void play();
    void pause();
    void stop();

    void setPosition(qint64 position);
    void setVolume(int volume);
    void setMuted(bool muted);

    void setPlaybackRate(qreal rate);

    void setMedia(const QMediaContent &media, QIODevice *stream = nullptr);
    void setPlaylist(QMediaPlaylist *playlist);

#ifndef QT_NO_BEARERMANAGEMENT
#ifndef Q_MOC_RUN // moc fails to parse the expanded macro
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
#endif
    QT_DEPRECATED_VERSION_5_15 void setNetworkConfigurations(const QList<QNetworkConfiguration> &configurations);
#ifndef Q_MOC_RUN // moc fails to parse the expanded macro
QT_WARNING_POP
#endif
#endif

Q_SIGNALS:
    void mediaChanged(const QMediaContent &media);
    void currentMediaChanged(const QMediaContent &media);

    void stateChanged(QMediaPlayer::State newState);
    void mediaStatusChanged(QMediaPlayer::MediaStatus status);

    void durationChanged(qint64 duration);
    void positionChanged(qint64 position);

    void volumeChanged(int volume);
    void mutedChanged(bool muted);
    void audioAvailableChanged(bool available);
    void videoAvailableChanged(bool videoAvailable);

    void bufferStatusChanged(int percentFilled);

    void seekableChanged(bool seekable);
    void playbackRateChanged(qreal rate);

    void audioRoleChanged(QAudio::Role role);
    void customAudioRoleChanged(const QString &role);

    void error(QMediaPlayer::Error error);

#ifndef QT_NO_BEARERMANAGEMENT
#ifndef Q_MOC_RUN // moc fails to parse the expanded macro
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
#endif
    QT_DEPRECATED_VERSION_5_15 void networkConfigurationChanged(const QNetworkConfiguration &configuration);
#ifndef Q_MOC_RUN // moc fails to parse the expanded macro
QT_WARNING_POP
#endif
#endif
public:
    bool bind(QObject *) override;
    void unbind(QObject *) override;

private:
    Q_DISABLE_COPY(QMediaPlayer)
    Q_DECLARE_PRIVATE(QMediaPlayer)
    Q_PRIVATE_SLOT(d_func(), void _q_stateChanged(QMediaPlayer::State))
    Q_PRIVATE_SLOT(d_func(), void _q_mediaStatusChanged(QMediaPlayer::MediaStatus))
    Q_PRIVATE_SLOT(d_func(), void _q_error(int, const QString &))
    Q_PRIVATE_SLOT(d_func(), void _q_updateMedia(const QMediaContent&))
    Q_PRIVATE_SLOT(d_func(), void _q_playlistDestroyed())
    Q_PRIVATE_SLOT(d_func(), void _q_handleMediaChanged(const QMediaContent&))
    Q_PRIVATE_SLOT(d_func(), void _q_handlePlaylistLoaded())
    Q_PRIVATE_SLOT(d_func(), void _q_handlePlaylistLoadFailed())
};

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QMediaPlayer::State)
Q_DECLARE_METATYPE(QMediaPlayer::MediaStatus)
Q_DECLARE_METATYPE(QMediaPlayer::Error)

Q_MEDIA_ENUM_DEBUG(QMediaPlayer, State)
Q_MEDIA_ENUM_DEBUG(QMediaPlayer, MediaStatus)
Q_MEDIA_ENUM_DEBUG(QMediaPlayer, Error)

#endif  // QMEDIAPLAYER_H
