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

#ifndef QSOUNDEFFECT_QAUDIO_H
#define QSOUNDEFFECT_QAUDIO_H

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

#include <QtCore/qobject.h>
#include <QtCore/qurl.h>
#include "qaudiooutput.h"
#include "qsamplecache_p.h"
#include "qsoundeffect.h"

QT_BEGIN_NAMESPACE

class QSoundEffectPrivate;

class PrivateSoundSource : public QIODevice
{
    friend class QSoundEffectPrivate;
    Q_OBJECT
public:
    PrivateSoundSource(QSoundEffectPrivate *s, const QAudioDeviceInfo &audioDevice = QAudioDeviceInfo());
    ~PrivateSoundSource() {}

    qint64 readData(char *data, qint64 len) override;
    qint64 writeData(const char *data, qint64 len) override;

private Q_SLOTS:
    void sampleReady();
    void decoderError();
    void stateChanged(QAudio::State);

private:
    QUrl m_url;
    int m_loopCount = 1;
    int m_runningCount = 0;
    bool m_playing = false;
    QSoundEffect::Status  m_status = QSoundEffect::Null;
    QAudioOutput *m_audioOutput = nullptr;
    QSample *m_sample = nullptr;
    bool m_muted = false;
    qreal m_volume = 1.0;
    bool m_sampleReady = false;
    qint64 m_offset = 0;
    QString m_category;
    QAudioDeviceInfo m_audioDevice;
    QSoundEffectPrivate *soundeffect = nullptr;
};


class QSoundEffectPrivate : public QObject
{
    friend class PrivateSoundSource;
    Q_OBJECT
public:

    explicit QSoundEffectPrivate(QObject *parent);
    explicit QSoundEffectPrivate(const QAudioDeviceInfo &audioDevice, QObject *parent);
    ~QSoundEffectPrivate();

    static QStringList supportedMimeTypes();

    QUrl source() const;
    void setSource(const QUrl &url);
    int loopCount() const;
    int loopsRemaining() const;
    void setLoopCount(int loopCount);
    qreal volume() const;
    void setVolume(qreal volume);
    bool isMuted() const;
    void setMuted(bool muted);
    bool isLoaded() const;
    bool isPlaying() const;
    QSoundEffect::Status status() const;

    void release();

    QString category() const;
    void setCategory(const QString &);

public Q_SLOTS:
    void play();
    void stop();

Q_SIGNALS:
    void loopsRemainingChanged();
    void volumeChanged();
    void mutedChanged();
    void loadedChanged();
    void playingChanged();
    void statusChanged();
    void categoryChanged();

private:
    void setStatus(QSoundEffect::Status status);
    void setPlaying(bool playing);
    void setLoopsRemaining(int loopsRemaining);

    PrivateSoundSource *d = nullptr;
};

QT_END_NAMESPACE

#endif // QSOUNDEFFECT_QAUDIO_H
