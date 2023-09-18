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



#ifndef QCAMERAVIEWFINDERSETTINGSCONTROL_H
#define QCAMERAVIEWFINDERSETTINGSCONTROL_H

#include <QtMultimedia/qmediacontrol.h>
#include <QtMultimedia/qcamera.h>

QT_BEGIN_NAMESPACE

// Required for QDoc workaround
class QString;

class Q_MULTIMEDIA_EXPORT QCameraViewfinderSettingsControl : public QMediaControl
{
    Q_OBJECT
public:
    enum ViewfinderParameter {
        Resolution,
        PixelAspectRatio,
        MinimumFrameRate,
        MaximumFrameRate,
        PixelFormat,
        UserParameter = 1000
    };

    ~QCameraViewfinderSettingsControl();

    virtual bool isViewfinderParameterSupported(ViewfinderParameter parameter) const = 0;
    virtual QVariant viewfinderParameter(ViewfinderParameter parameter) const = 0;
    virtual void setViewfinderParameter(ViewfinderParameter parameter, const QVariant &value) = 0;

protected:
    explicit QCameraViewfinderSettingsControl(QObject *parent = nullptr);
};

#define QCameraViewfinderSettingsControl_iid "org.qt-project.qt.cameraviewfindersettingscontrol/5.0"
Q_MEDIA_DECLARE_CONTROL(QCameraViewfinderSettingsControl, QCameraViewfinderSettingsControl_iid)


// Required for QDoc workaround
class QString;

class Q_MULTIMEDIA_EXPORT QCameraViewfinderSettingsControl2 : public QMediaControl
{
    Q_OBJECT
public:
    virtual ~QCameraViewfinderSettingsControl2();

    virtual QList<QCameraViewfinderSettings> supportedViewfinderSettings() const = 0;

    virtual QCameraViewfinderSettings viewfinderSettings() const = 0;
    virtual void setViewfinderSettings(const QCameraViewfinderSettings &settings) = 0;

protected:
    explicit QCameraViewfinderSettingsControl2(QObject *parent = nullptr);
};

#define QCameraViewfinderSettingsControl2_iid "org.qt-project.qt.cameraviewfindersettingscontrol2/5.5"
Q_MEDIA_DECLARE_CONTROL(QCameraViewfinderSettingsControl2, QCameraViewfinderSettingsControl2_iid)

QT_END_NAMESPACE

#endif // QCAMERAVIEWFINDERSETTINGSCONTROL_H
