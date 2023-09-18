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

#ifndef QAUDIODECODERCONTROL_H
#define QAUDIODECODERCONTROL_H

#include <QtMultimedia/qmediacontrol.h>
#include <QtMultimedia/qaudiodecoder.h>

#include <QtCore/qpair.h>

#include <QtMultimedia/qaudiobuffer.h>

QT_BEGIN_NAMESPACE

class QIODevice;
class Q_MULTIMEDIA_EXPORT QAudioDecoderControl : public QMediaControl
{
    Q_OBJECT

public:
    ~QAudioDecoderControl();

    virtual QAudioDecoder::State state() const = 0;

    virtual QString sourceFilename() const = 0;
    virtual void setSourceFilename(const QString &fileName) = 0;

    virtual QIODevice* sourceDevice() const = 0;
    virtual void setSourceDevice(QIODevice *device) = 0;

    virtual void start() = 0;
    virtual void stop() = 0;

    virtual QAudioFormat audioFormat() const = 0;
    virtual void setAudioFormat(const QAudioFormat &format) = 0;

    virtual QAudioBuffer read() = 0;
    virtual bool bufferAvailable() const = 0;

    virtual qint64 position() const = 0;
    virtual qint64 duration() const = 0;

Q_SIGNALS:
    void stateChanged(QAudioDecoder::State newState);
    void formatChanged(const QAudioFormat &format);
    void sourceChanged();

    void error(int error, const QString &errorString);

    void bufferReady();
    void bufferAvailableChanged(bool available);
    void finished();

    void positionChanged(qint64 position);
    void durationChanged(qint64 duration);

protected:
    explicit QAudioDecoderControl(QObject *parent = nullptr);
};

#define QAudioDecoderControl_iid "org.qt-project.qt.audiodecodercontrol/5.0"
Q_MEDIA_DECLARE_CONTROL(QAudioDecoderControl, QAudioDecoderControl_iid)

QT_END_NAMESPACE

#endif  // QAUDIODECODERCONTROL_H
