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


#ifndef QMEDIAPLAYLISTCONTROL_P_H
#define QMEDIAPLAYLISTCONTROL_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API. It exists purely as an
// implementation detail. This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <QtCore/qobject.h>
#include "qmediacontrol.h"
#include <private/qmediaplaylistnavigator_p.h>


QT_BEGIN_NAMESPACE


class QMediaPlaylistProvider;

class Q_MULTIMEDIA_EXPORT QMediaPlaylistControl : public QMediaControl
{
    Q_OBJECT

public:
    virtual ~QMediaPlaylistControl();

    virtual QMediaPlaylistProvider* playlistProvider() const = 0;
    virtual bool setPlaylistProvider(QMediaPlaylistProvider *playlist) = 0;

    virtual int currentIndex() const = 0;
    virtual void setCurrentIndex(int position) = 0;
    virtual int nextIndex(int steps) const = 0;
    virtual int previousIndex(int steps) const = 0;

    virtual void next() = 0;
    virtual void previous() = 0;

    virtual QMediaPlaylist::PlaybackMode playbackMode() const = 0;
    virtual void setPlaybackMode(QMediaPlaylist::PlaybackMode mode) = 0;

Q_SIGNALS:
    void playlistProviderChanged();
    void currentIndexChanged(int position);
    void currentMediaChanged(const QMediaContent&);
    void playbackModeChanged(QMediaPlaylist::PlaybackMode mode);

protected:
    QMediaPlaylistControl(QObject *parent = nullptr);
};

#define QMediaPlaylistControl_iid "org.qt-project.qt.mediaplaylistcontrol/5.0"
Q_MEDIA_DECLARE_CONTROL(QMediaPlaylistControl, QMediaPlaylistControl_iid)

QT_END_NAMESPACE


#endif // QMEDIAPLAYLISTCONTROL_P_H
