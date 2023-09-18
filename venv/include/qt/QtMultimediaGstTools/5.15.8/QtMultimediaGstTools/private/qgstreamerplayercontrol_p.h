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

#ifndef QGSTREAMERPLAYERCONTROL_P_H
#define QGSTREAMERPLAYERCONTROL_P_H

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

#include <QtCore/qstack.h>
#include <qmediaplayercontrol.h>
#include <private/qgsttools_global_p.h>

QT_BEGIN_NAMESPACE

class QMediaPlayerResourceSetInterface;
class QGstreamerPlayerSession;
class Q_GSTTOOLS_EXPORT QGstreamerPlayerControl : public QMediaPlayerControl
{
    Q_OBJECT

public:
    QGstreamerPlayerControl(QGstreamerPlayerSession *session, QObject *parent = 0);
    ~QGstreamerPlayerControl();

    QGstreamerPlayerSession *session() { return m_session; }

    QMediaPlayer::State state() const override;
    QMediaPlayer::MediaStatus mediaStatus() const override;

    qint64 position() const override;
    qint64 duration() const override;

    int bufferStatus() const override;

    int volume() const override;
    bool isMuted() const override;

    bool isAudioAvailable() const override;
    bool isVideoAvailable() const override;
    void setVideoOutput(QObject *output);

    bool isSeekable() const override;
    QMediaTimeRange availablePlaybackRanges() const override;

    qreal playbackRate() const override;
    void setPlaybackRate(qreal rate) override;

    QMediaContent media() const override;
    const QIODevice *mediaStream() const override;
    void setMedia(const QMediaContent&, QIODevice *) override;

    QMediaPlayerResourceSetInterface* resources() const;

public Q_SLOTS:
    void setPosition(qint64 pos) override;

    void play() override;
    void pause() override;
    void stop() override;

    void setVolume(int volume) override;
    void setMuted(bool muted) override;

private Q_SLOTS:
    void updateSessionState(QMediaPlayer::State state);
    void updateMediaStatus();
    void processEOS();
    void setBufferProgress(int progress);

    void handleInvalidMedia();

    void handleResourcesGranted();
    void handleResourcesLost();
    void handleResourcesDenied();

private:
    void playOrPause(QMediaPlayer::State state);

    void pushState();
    void popAndNotifyState();

    QGstreamerPlayerSession *m_session = nullptr;
    QMediaPlayer::State m_userRequestedState = QMediaPlayer::StoppedState;
    QMediaPlayer::State m_currentState = QMediaPlayer::StoppedState;
    QMediaPlayer::MediaStatus m_mediaStatus = QMediaPlayer::NoMedia;
    QStack<QMediaPlayer::State> m_stateStack;
    QStack<QMediaPlayer::MediaStatus> m_mediaStatusStack;

    int m_bufferProgress = -1;
    qint64 m_pendingSeekPosition = -1;
    bool m_setMediaPending = false;
    QMediaContent m_currentResource;
    QIODevice *m_stream = nullptr;

    QMediaPlayerResourceSetInterface *m_resources = nullptr;
};

QT_END_NAMESPACE

#endif
