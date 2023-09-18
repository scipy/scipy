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

#ifndef QAUDIOSYSTEM_H
#define QAUDIOSYSTEM_H

#include <QtMultimedia/qtmultimediaglobal.h>
#include <QtMultimedia/qmultimedia.h>

#include <QtMultimedia/qaudio.h>
#include <QtMultimedia/qaudioformat.h>
#include <QtMultimedia/qaudiodeviceinfo.h>

QT_BEGIN_NAMESPACE

class QIODevice;

class Q_MULTIMEDIA_EXPORT QAbstractAudioDeviceInfo : public QObject
{
    Q_OBJECT

public:
    virtual QAudioFormat preferredFormat() const = 0;
    virtual bool isFormatSupported(const QAudioFormat &format) const = 0;
    virtual QString deviceName() const = 0;
    virtual QStringList supportedCodecs() = 0;
    virtual QList<int> supportedSampleRates() = 0;
    virtual QList<int> supportedChannelCounts() = 0;
    virtual QList<int> supportedSampleSizes() = 0;
    virtual QList<QAudioFormat::Endian> supportedByteOrders() = 0;
    virtual QList<QAudioFormat::SampleType> supportedSampleTypes() = 0;
};

class Q_MULTIMEDIA_EXPORT QAbstractAudioOutput : public QObject
{
    Q_OBJECT

public:
    virtual void start(QIODevice *device) = 0;
    virtual QIODevice* start() = 0;
    virtual void stop() = 0;
    virtual void reset() = 0;
    virtual void suspend() = 0;
    virtual void resume() = 0;
    virtual int bytesFree() const = 0;
    virtual int periodSize() const = 0;
    virtual void setBufferSize(int value) = 0;
    virtual int bufferSize() const = 0;
    virtual void setNotifyInterval(int milliSeconds) = 0;
    virtual int notifyInterval() const = 0;
    virtual qint64 processedUSecs() const = 0;
    virtual qint64 elapsedUSecs() const = 0;
    virtual QAudio::Error error() const = 0;
    virtual QAudio::State state() const = 0;
    virtual void setFormat(const QAudioFormat& fmt) = 0;
    virtual QAudioFormat format() const = 0;
    virtual void setVolume(qreal) {}
    virtual qreal volume() const { return 1.0; }
    virtual QString category() const { return QString(); }
    virtual void setCategory(const QString &) { }

Q_SIGNALS:
    void errorChanged(QAudio::Error error);
    void stateChanged(QAudio::State state);
    void notify();
};

class Q_MULTIMEDIA_EXPORT QAbstractAudioInput : public QObject
{
    Q_OBJECT

public:
    virtual void start(QIODevice *device) = 0;
    virtual QIODevice* start() = 0;
    virtual void stop() = 0;
    virtual void reset() = 0;
    virtual void suspend()  = 0;
    virtual void resume() = 0;
    virtual int bytesReady() const = 0;
    virtual int periodSize() const = 0;
    virtual void setBufferSize(int value) = 0;
    virtual int bufferSize() const = 0;
    virtual void setNotifyInterval(int milliSeconds) = 0;
    virtual int notifyInterval() const = 0;
    virtual qint64 processedUSecs() const = 0;
    virtual qint64 elapsedUSecs() const = 0;
    virtual QAudio::Error error() const = 0;
    virtual QAudio::State state() const = 0;
    virtual void setFormat(const QAudioFormat& fmt) = 0;
    virtual QAudioFormat format() const = 0;
    virtual void setVolume(qreal) = 0;
    virtual qreal volume() const = 0;

Q_SIGNALS:
    void errorChanged(QAudio::Error error);
    void stateChanged(QAudio::State state);
    void notify();
};

QT_END_NAMESPACE

#endif // QAUDIOSYSTEM_H
