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


#ifndef QMEDIASTREAMSCONTROL_H
#define QMEDIASTREAMSCONTROL_H

#include <QtMultimedia/qmediacontrol.h>
#include <QtMultimedia/qmultimedia.h>
#include <QtMultimedia/qtmultimediaglobal.h>
#include <QtMultimedia/qmediaenumdebug.h>

QT_BEGIN_NAMESPACE

// Required for QDoc workaround
class QString;

class Q_MULTIMEDIA_EXPORT QMediaStreamsControl : public QMediaControl
{
    Q_OBJECT
    Q_ENUMS(SteamType)
public:
    enum StreamType { UnknownStream, VideoStream, AudioStream, SubPictureStream, DataStream };

    virtual ~QMediaStreamsControl();

    virtual int streamCount() = 0;
    virtual StreamType streamType(int streamNumber) = 0;

    virtual QVariant metaData(int streamNumber, const QString &key) = 0;

    virtual bool isActive(int streamNumber) = 0;
    virtual void setActive(int streamNumber, bool state) = 0;

Q_SIGNALS:
    void streamsChanged();
    void activeStreamsChanged();

protected:
    explicit QMediaStreamsControl(QObject *parent = nullptr);
};

#define QMediaStreamsControl_iid "org.qt-project.qt.mediastreamscontrol/5.0"
Q_MEDIA_DECLARE_CONTROL(QMediaStreamsControl, QMediaStreamsControl_iid)

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QMediaStreamsControl::StreamType)

Q_MEDIA_ENUM_DEBUG(QMediaStreamsControl, StreamType)

#endif // QMEDIASTREAMSCONTROL_H

