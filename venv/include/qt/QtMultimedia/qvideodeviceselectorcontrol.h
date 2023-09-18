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

#ifndef QVIDEODEVICESELECTORCONTROL_H
#define QVIDEODEVICESELECTORCONTROL_H

#include <QtMultimedia/qmediacontrol.h>

QT_BEGIN_NAMESPACE

// Required for QDoc workaround
class QString;

class Q_MULTIMEDIA_EXPORT QVideoDeviceSelectorControl : public QMediaControl
{
    Q_OBJECT

public:
    virtual ~QVideoDeviceSelectorControl();

    virtual int deviceCount() const = 0;

    virtual QString deviceName(int index) const = 0;
    virtual QString deviceDescription(int index) const = 0;

    virtual int defaultDevice() const = 0;
    virtual int selectedDevice() const = 0;

public Q_SLOTS:
    virtual void setSelectedDevice(int index) = 0;

Q_SIGNALS:
    void selectedDeviceChanged(int index);
    void selectedDeviceChanged(const QString &name);
    void devicesChanged();

protected:
    explicit QVideoDeviceSelectorControl(QObject *parent = nullptr);
};

#define QVideoDeviceSelectorControl_iid "org.qt-project.qt.videodeviceselectorcontrol/5.0"
Q_MEDIA_DECLARE_CONTROL(QVideoDeviceSelectorControl, QVideoDeviceSelectorControl_iid)

QT_END_NAMESPACE

#endif // QVIDEODEVICESELECTORCONTROL_H
