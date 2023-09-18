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

#ifndef QSOUNDEFFECT_H
#define QSOUNDEFFECT_H

#include <QtMultimedia/qtmultimediaglobal.h>
#include <QtCore/qobject.h>
#include <QtCore/qurl.h>
#include <QtCore/qstringlist.h>


QT_BEGIN_NAMESPACE


class QSoundEffectPrivate;
class QAudioDeviceInfo;

class Q_MULTIMEDIA_EXPORT QSoundEffect : public QObject
{
    Q_OBJECT
    Q_CLASSINFO("DefaultMethod", "play()")
    Q_PROPERTY(QUrl source READ source WRITE setSource NOTIFY sourceChanged)
    Q_PROPERTY(int loops READ loopCount WRITE setLoopCount NOTIFY loopCountChanged)
    Q_PROPERTY(int loopsRemaining READ loopsRemaining NOTIFY loopsRemainingChanged)
    Q_PROPERTY(qreal volume READ volume WRITE setVolume NOTIFY volumeChanged)
    Q_PROPERTY(bool muted READ isMuted WRITE setMuted NOTIFY mutedChanged)
    Q_PROPERTY(bool playing READ isPlaying NOTIFY playingChanged)
    Q_PROPERTY(Status status READ status NOTIFY statusChanged)
    Q_PROPERTY(QString category READ category WRITE setCategory NOTIFY categoryChanged)
    Q_ENUMS(Loop)
    Q_ENUMS(Status)

public:
    enum Loop
    {
        Infinite = -2
    };

    enum Status
    {
        Null,
        Loading,
        Ready,
        Error
    };

    explicit QSoundEffect(QObject *parent = nullptr);
    explicit QSoundEffect(const QAudioDeviceInfo &audioDevice, QObject *parent = nullptr);
    ~QSoundEffect();

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
    Status status() const;

    QString category() const;
    void setCategory(const QString &category);

Q_SIGNALS:
    void sourceChanged();
    void loopCountChanged();
    void loopsRemainingChanged();
    void volumeChanged();
    void mutedChanged();
    void loadedChanged();
    void playingChanged();
    void statusChanged();
    void categoryChanged();

public Q_SLOTS:
    void play();
    void stop();

private:
    Q_DISABLE_COPY(QSoundEffect)
    QSoundEffectPrivate *d = nullptr;
};

QT_END_NAMESPACE


#endif // QSOUNDEFFECT_H
