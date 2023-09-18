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

#ifndef QMEDIARECORDER_H
#define QMEDIARECORDER_H

#include <QtMultimedia/qmultimedia.h>
#include <QtMultimedia/qmediaobject.h>
#include <QtMultimedia/qmediaencodersettings.h>
#include <QtMultimedia/qmediabindableinterface.h>
#include <QtMultimedia/qmediaenumdebug.h>

#include <QtCore/qpair.h>

QT_BEGIN_NAMESPACE

class QUrl;
class QSize;
class QAudioFormat;
QT_END_NAMESPACE

QT_BEGIN_NAMESPACE

class QMediaRecorderService;
class QAudioEncoderSettings;
class QVideoEncoderSettings;

class QMediaRecorderPrivate;
class Q_MULTIMEDIA_EXPORT QMediaRecorder : public QObject, public QMediaBindableInterface
{
    Q_OBJECT
    Q_INTERFACES(QMediaBindableInterface)
    Q_ENUMS(State)
    Q_ENUMS(Status)
    Q_ENUMS(Error)
    Q_PROPERTY(QMediaRecorder::State state READ state NOTIFY stateChanged)
    Q_PROPERTY(QMediaRecorder::Status status READ status NOTIFY statusChanged)
    Q_PROPERTY(qint64 duration READ duration NOTIFY durationChanged)
    Q_PROPERTY(QUrl outputLocation READ outputLocation WRITE setOutputLocation)
    Q_PROPERTY(QUrl actualLocation READ actualLocation NOTIFY actualLocationChanged)
    Q_PROPERTY(bool muted READ isMuted WRITE setMuted NOTIFY mutedChanged)
    Q_PROPERTY(qreal volume READ volume WRITE setVolume NOTIFY volumeChanged)
    Q_PROPERTY(bool metaDataAvailable READ isMetaDataAvailable NOTIFY metaDataAvailableChanged)
    Q_PROPERTY(bool metaDataWritable READ isMetaDataWritable NOTIFY metaDataWritableChanged)
public:

    enum State
    {
        StoppedState,
        RecordingState,
        PausedState
    };

    enum Status {
        UnavailableStatus,
        UnloadedStatus,
        LoadingStatus,
        LoadedStatus,
        StartingStatus,
        RecordingStatus,
        PausedStatus,
        FinalizingStatus
    };

    enum Error
    {
        NoError,
        ResourceError,
        FormatError,
        OutOfSpaceError
    };

    explicit QMediaRecorder(QMediaObject *mediaObject, QObject *parent = nullptr);
    ~QMediaRecorder();

    QMediaObject *mediaObject() const override;

    bool isAvailable() const;
    QMultimedia::AvailabilityStatus availability() const;

    QUrl outputLocation() const;
    bool setOutputLocation(const QUrl &location);

    QUrl actualLocation() const;

    State state() const;
    Status status() const;

    Error error() const;
    QString errorString() const;

    qint64 duration() const;

    bool isMuted() const;
    qreal volume() const;

    QStringList supportedContainers() const;
    QString containerDescription(const QString &format) const;

    QStringList supportedAudioCodecs() const;
    QString audioCodecDescription(const QString &codecName) const;

    QList<int> supportedAudioSampleRates(const QAudioEncoderSettings &settings = QAudioEncoderSettings(),
                                         bool *continuous = nullptr) const;

    QStringList supportedVideoCodecs() const;
    QString videoCodecDescription(const QString &codecName) const;

    QList<QSize> supportedResolutions(const QVideoEncoderSettings &settings = QVideoEncoderSettings(),
                                      bool *continuous = nullptr) const;

    QList<qreal> supportedFrameRates(const QVideoEncoderSettings &settings = QVideoEncoderSettings(),
                                     bool *continuous = nullptr) const;

    QAudioEncoderSettings audioSettings() const;
    QVideoEncoderSettings videoSettings() const;
    QString containerFormat() const;

    void setAudioSettings(const QAudioEncoderSettings &audioSettings);
    void setVideoSettings(const QVideoEncoderSettings &videoSettings);
    void setContainerFormat(const QString &container);

    void setEncodingSettings(const QAudioEncoderSettings &audioSettings,
                             const QVideoEncoderSettings &videoSettings = QVideoEncoderSettings(),
                             const QString &containerMimeType = QString());

    bool isMetaDataAvailable() const;
    bool isMetaDataWritable() const;

    QVariant metaData(const QString &key) const;
    void setMetaData(const QString &key, const QVariant &value);
    QStringList availableMetaData() const;

public Q_SLOTS:
    void record();
    void pause();
    void stop();
    void setMuted(bool muted);
    void setVolume(qreal volume);

Q_SIGNALS:
    void stateChanged(QMediaRecorder::State state);
    void statusChanged(QMediaRecorder::Status status);
    void durationChanged(qint64 duration);
    void mutedChanged(bool muted);
    void volumeChanged(qreal volume);
    void actualLocationChanged(const QUrl &location);

    void error(QMediaRecorder::Error error);

    void metaDataAvailableChanged(bool available);
    void metaDataWritableChanged(bool writable);
    void metaDataChanged();
    void metaDataChanged(const QString &key, const QVariant &value);

    void availabilityChanged(bool available);
    void availabilityChanged(QMultimedia::AvailabilityStatus availability);

protected:
    QMediaRecorder(QMediaRecorderPrivate &dd, QMediaObject *mediaObject, QObject *parent = nullptr);
    bool setMediaObject(QMediaObject *object) override;

    QMediaRecorderPrivate *d_ptr;
private:
    Q_DISABLE_COPY(QMediaRecorder)
    Q_DECLARE_PRIVATE(QMediaRecorder)
    Q_PRIVATE_SLOT(d_func(), void _q_stateChanged(QMediaRecorder::State))
    Q_PRIVATE_SLOT(d_func(), void _q_error(int, const QString &))
    Q_PRIVATE_SLOT(d_func(), void _q_serviceDestroyed())
    Q_PRIVATE_SLOT(d_func(), void _q_notify())
    Q_PRIVATE_SLOT(d_func(), void _q_updateActualLocation(const QUrl &))
    Q_PRIVATE_SLOT(d_func(), void _q_updateNotifyInterval(int))
    Q_PRIVATE_SLOT(d_func(), void _q_applySettings())
    Q_PRIVATE_SLOT(d_func(), void _q_availabilityChanged(QMultimedia::AvailabilityStatus))
};

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QMediaRecorder::State)
Q_DECLARE_METATYPE(QMediaRecorder::Status)
Q_DECLARE_METATYPE(QMediaRecorder::Error)

Q_MEDIA_ENUM_DEBUG(QMediaRecorder, State)
Q_MEDIA_ENUM_DEBUG(QMediaRecorder, Status)
Q_MEDIA_ENUM_DEBUG(QMediaRecorder, Error)

#endif  // QMEDIARECORDER_H
