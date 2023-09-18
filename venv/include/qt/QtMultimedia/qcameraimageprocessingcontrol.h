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

#ifndef QCAMERAIMAGEPROCESSINGCONTROL_H
#define QCAMERAIMAGEPROCESSINGCONTROL_H

#include <QtMultimedia/qmediacontrol.h>
#include <QtMultimedia/qmediaobject.h>

#include <QtMultimedia/qcamera.h>
#include <QtMultimedia/qmediaenumdebug.h>

QT_BEGIN_NAMESPACE

// Required for QDoc workaround
class QString;

class Q_MULTIMEDIA_EXPORT QCameraImageProcessingControl : public QMediaControl
{
    Q_OBJECT
    Q_ENUMS(ProcessingParameter)

public:
    ~QCameraImageProcessingControl();

    enum ProcessingParameter {
        WhiteBalancePreset,
        ColorTemperature,
        Contrast,
        Saturation,
        Brightness,
        Sharpening,
        Denoising,
        ContrastAdjustment,
        SaturationAdjustment,
        BrightnessAdjustment,
        SharpeningAdjustment,
        DenoisingAdjustment,
        ColorFilter,
        ExtendedParameter = 1000
    };

    virtual bool isParameterSupported(ProcessingParameter) const = 0;
    virtual bool isParameterValueSupported(ProcessingParameter parameter, const QVariant &value) const = 0;
    virtual QVariant parameter(ProcessingParameter parameter) const = 0;
    virtual void setParameter(ProcessingParameter parameter, const QVariant &value) = 0;

protected:
    explicit QCameraImageProcessingControl(QObject *parent = nullptr);
};

#define QCameraImageProcessingControl_iid "org.qt-project.qt.cameraimageprocessingcontrol/5.0"
Q_MEDIA_DECLARE_CONTROL(QCameraImageProcessingControl, QCameraImageProcessingControl_iid)

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QCameraImageProcessingControl::ProcessingParameter)

Q_MEDIA_ENUM_DEBUG(QCameraImageProcessingControl, ProcessingParameter)

#endif

