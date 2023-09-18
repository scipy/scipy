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

#ifndef QCAMERA_H
#define QCAMERA_H

#include <QtCore/qstringlist.h>
#include <QtCore/qpair.h>
#include <QtCore/qsize.h>
#include <QtCore/qpoint.h>
#include <QtCore/qrect.h>

#include <QtMultimedia/qmediacontrol.h>
#include <QtMultimedia/qmediaobject.h>
#include <QtMultimedia/qmediaservice.h>

#include <QtMultimedia/qcameraexposure.h>
#include <QtMultimedia/qcamerafocus.h>
#include <QtMultimedia/qcameraimageprocessing.h>
#include <QtMultimedia/qcameraviewfindersettings.h>

#include <QtMultimedia/qmediaenumdebug.h>

QT_BEGIN_NAMESPACE


class QAbstractVideoSurface;
class QVideoWidget;
class QGraphicsVideoItem;
class QCameraInfo;

class QCameraPrivate;
class Q_MULTIMEDIA_EXPORT QCamera : public QMediaObject
{
    Q_OBJECT
    Q_PROPERTY(QCamera::State state READ state NOTIFY stateChanged)
    Q_PROPERTY(QCamera::Status status READ status NOTIFY statusChanged)
    Q_PROPERTY(QCamera::CaptureModes captureMode READ captureMode WRITE setCaptureMode NOTIFY captureModeChanged)
    Q_PROPERTY(QCamera::LockStatus lockStatus READ lockStatus NOTIFY lockStatusChanged)

    Q_ENUMS(Status)
    Q_ENUMS(State)
    Q_ENUMS(CaptureMode)
    Q_ENUMS(Error)
    Q_ENUMS(LockStatus)
    Q_ENUMS(LockChangeReason)
    Q_ENUMS(LockType)
    Q_ENUMS(Position)
public:
    struct FrameRateRange
    {
        Q_DECL_CONSTEXPR FrameRateRange() Q_DECL_NOTHROW
            : minimumFrameRate(0)
            , maximumFrameRate(0)
        { }

        Q_DECL_CONSTEXPR FrameRateRange(qreal minimum, qreal maximum) Q_DECL_NOTHROW
            : minimumFrameRate(minimum)
            , maximumFrameRate(maximum)
        { }

        qreal minimumFrameRate;
        qreal maximumFrameRate;
    };

    enum Status {
        UnavailableStatus,
        UnloadedStatus,
        LoadingStatus,
        UnloadingStatus,
        LoadedStatus,
        StandbyStatus,
        StartingStatus,
        StoppingStatus,
        ActiveStatus
    };

    enum State {
        UnloadedState,
        LoadedState,
        ActiveState
    };

    enum CaptureMode
    {
        CaptureViewfinder = 0,
        CaptureStillImage = 0x01,
        CaptureVideo = 0x02
    };
    Q_DECLARE_FLAGS(CaptureModes, CaptureMode)

    enum Error
    {
        NoError,
        CameraError,
        InvalidRequestError,
        ServiceMissingError,
        NotSupportedFeatureError
    };

    enum LockStatus
    {
        Unlocked,
        Searching,
        Locked
    };

    enum LockChangeReason {
        UserRequest,
        LockAcquired,
        LockFailed,
        LockLost,
        LockTemporaryLost
    };

    enum LockType
    {
        NoLock = 0,
        LockExposure = 0x01,
        LockWhiteBalance = 0x02,
        LockFocus = 0x04
    };
    Q_DECLARE_FLAGS(LockTypes, LockType)

    enum Position
    {
        UnspecifiedPosition,
        BackFace,
        FrontFace
    };

    explicit QCamera(QObject *parent = nullptr);
    explicit QCamera(const QByteArray& deviceName, QObject *parent = nullptr);
    explicit QCamera(const QCameraInfo& cameraInfo, QObject *parent = nullptr);
    explicit QCamera(QCamera::Position position, QObject *parent = nullptr);
    ~QCamera();

#if QT_DEPRECATED_SINCE(5, 3)
    QT_DEPRECATED static QList<QByteArray> availableDevices();
    QT_DEPRECATED static QString deviceDescription(const QByteArray &device);
#endif

    QMultimedia::AvailabilityStatus availability() const override;

    State state() const;
    Status status() const;

    CaptureModes captureMode() const;
    bool isCaptureModeSupported(CaptureModes mode) const;

    QCameraExposure *exposure() const;
    QCameraFocus *focus() const;
    QCameraImageProcessing *imageProcessing() const;

    void setViewfinder(QVideoWidget *viewfinder);
    void setViewfinder(QGraphicsVideoItem *viewfinder);
    void setViewfinder(QAbstractVideoSurface *surface);

    QCameraViewfinderSettings viewfinderSettings() const;
    void setViewfinderSettings(const QCameraViewfinderSettings &settings);

    QList<QCameraViewfinderSettings> supportedViewfinderSettings(
            const QCameraViewfinderSettings &settings = QCameraViewfinderSettings()) const;

    QList<QSize> supportedViewfinderResolutions(
            const QCameraViewfinderSettings &settings = QCameraViewfinderSettings()) const;

    QList<FrameRateRange> supportedViewfinderFrameRateRanges(
            const QCameraViewfinderSettings &settings = QCameraViewfinderSettings()) const;

    QList<QVideoFrame::PixelFormat> supportedViewfinderPixelFormats(
            const QCameraViewfinderSettings &settings = QCameraViewfinderSettings()) const;

    Error error() const;
    QString errorString() const;

    QCamera::LockTypes supportedLocks() const;
    QCamera::LockTypes requestedLocks() const;

    QCamera::LockStatus lockStatus() const;
    QCamera::LockStatus lockStatus(QCamera::LockType lock) const;

public Q_SLOTS:
    void setCaptureMode(QCamera::CaptureModes mode);

    void load();
    void unload();

    void start();
    void stop();

    void searchAndLock();
    void unlock();

    void searchAndLock(QCamera::LockTypes locks);
    void unlock(QCamera::LockTypes locks);

Q_SIGNALS:
    void stateChanged(QCamera::State state);
    void captureModeChanged(QCamera::CaptureModes);
    void statusChanged(QCamera::Status status);

    void locked();
    void lockFailed();

    void lockStatusChanged(QCamera::LockStatus status, QCamera::LockChangeReason reason);
    void lockStatusChanged(QCamera::LockType lock, QCamera::LockStatus status, QCamera::LockChangeReason reason);

#if QT_DEPRECATED_SINCE(5,15)
    void error(QCamera::Error);
#endif
    void errorOccurred(QCamera::Error);

private:
    Q_DISABLE_COPY(QCamera)
    Q_DECLARE_PRIVATE(QCamera)
    Q_PRIVATE_SLOT(d_func(), void _q_preparePropertyChange(int))
    Q_PRIVATE_SLOT(d_func(), void _q_restartCamera())
    Q_PRIVATE_SLOT(d_func(), void _q_error(int, const QString &))
    Q_PRIVATE_SLOT(d_func(), void _q_updateLockStatus(QCamera::LockType, QCamera::LockStatus, QCamera::LockChangeReason))
    Q_PRIVATE_SLOT(d_func(), void _q_updateState(QCamera::State))
    friend class QCameraInfo;
};

Q_DECLARE_OPERATORS_FOR_FLAGS(QCamera::LockTypes)

QT_WARNING_PUSH
QT_WARNING_DISABLE_CLANG("-Wfloat-equal")
QT_WARNING_DISABLE_GCC("-Wfloat-equal")

Q_DECL_CONSTEXPR Q_INLINE_TEMPLATE bool operator==(const QCamera::FrameRateRange &r1, const QCamera::FrameRateRange &r2) Q_DECL_NOTHROW
{
    return qFuzzyCompare(r1.minimumFrameRate, r2.minimumFrameRate)
        && qFuzzyCompare(r1.maximumFrameRate, r2.maximumFrameRate);
}

QT_WARNING_POP

Q_DECL_CONSTEXPR Q_INLINE_TEMPLATE bool operator!=(const QCamera::FrameRateRange &r1, const QCamera::FrameRateRange &r2) Q_DECL_NOTHROW
{ return !(r1 == r2); }

Q_DECLARE_TYPEINFO(QCamera::FrameRateRange, Q_PRIMITIVE_TYPE);

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QCamera::State)
Q_DECLARE_METATYPE(QCamera::Status)
Q_DECLARE_METATYPE(QCamera::Error)
Q_DECLARE_METATYPE(QCamera::CaptureMode)
Q_DECLARE_METATYPE(QCamera::CaptureModes)
Q_DECLARE_METATYPE(QCamera::LockType)
Q_DECLARE_METATYPE(QCamera::LockStatus)
Q_DECLARE_METATYPE(QCamera::LockChangeReason)
Q_DECLARE_METATYPE(QCamera::Position)

Q_MEDIA_ENUM_DEBUG(QCamera, State)
Q_MEDIA_ENUM_DEBUG(QCamera, Status)
Q_MEDIA_ENUM_DEBUG(QCamera, Error)
Q_MEDIA_ENUM_DEBUG(QCamera, CaptureMode)
Q_MEDIA_ENUM_DEBUG(QCamera, LockType)
Q_MEDIA_ENUM_DEBUG(QCamera, LockStatus)
Q_MEDIA_ENUM_DEBUG(QCamera, LockChangeReason)
Q_MEDIA_ENUM_DEBUG(QCamera, Position)

#endif  // QCAMERA_H
