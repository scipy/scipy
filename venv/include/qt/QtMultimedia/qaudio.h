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


#ifndef QAUDIO_H
#define QAUDIO_H

#include <QtMultimedia/qtmultimediaglobal.h>
#include <QtMultimedia/qmultimedia.h>

#include <QtCore/qmetatype.h>

QT_BEGIN_NAMESPACE

//QTM_SYNC_HEADER_EXPORT QAudio

// Class forward declaration required for QDoc bug
class QString;
namespace QAudio
{
    enum Error { NoError, OpenError, IOError, UnderrunError, FatalError };
    enum State { ActiveState, SuspendedState, StoppedState, IdleState, InterruptedState };
    enum Mode { AudioInput, AudioOutput };

    enum Role {
        UnknownRole,
        MusicRole,
        VideoRole,
        VoiceCommunicationRole,
        AlarmRole,
        NotificationRole,
        RingtoneRole,
        AccessibilityRole,
        SonificationRole,
        GameRole,
        CustomRole
    };

    enum VolumeScale {
        LinearVolumeScale,
        CubicVolumeScale,
        LogarithmicVolumeScale,
        DecibelVolumeScale
    };

    Q_MULTIMEDIA_EXPORT qreal convertVolume(qreal volume, VolumeScale from, VolumeScale to);
}

#ifndef QT_NO_DEBUG_STREAM
Q_MULTIMEDIA_EXPORT QDebug operator<<(QDebug dbg, QAudio::Error error);
Q_MULTIMEDIA_EXPORT QDebug operator<<(QDebug dbg, QAudio::State state);
Q_MULTIMEDIA_EXPORT QDebug operator<<(QDebug dbg, QAudio::Mode mode);
Q_MULTIMEDIA_EXPORT QDebug operator<<(QDebug dbg, QAudio::Role role);
Q_MULTIMEDIA_EXPORT QDebug operator<<(QDebug dbg, QAudio::VolumeScale role);
#endif

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QAudio::Error)
Q_DECLARE_METATYPE(QAudio::State)
Q_DECLARE_METATYPE(QAudio::Mode)
Q_DECLARE_METATYPE(QAudio::Role)
Q_DECLARE_METATYPE(QAudio::VolumeScale)

#endif // QAUDIO_H
