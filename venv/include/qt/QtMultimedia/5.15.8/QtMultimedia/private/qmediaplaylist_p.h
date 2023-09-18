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

#ifndef QMEDIAPLAYLIST_P_H
#define QMEDIAPLAYLIST_P_H

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

#include "qmediaplaylist.h"
#include "qmediaplaylistcontrol_p.h"
#include "qmediaplayer.h"
#include "qmediaplayercontrol.h"
#include "qmedianetworkplaylistprovider_p.h"
#include "qmediaobject_p.h"

#include <QtCore/qdebug.h>

#ifdef Q_MOC_RUN
# pragma Q_MOC_EXPAND_MACROS
#endif

QT_BEGIN_NAMESPACE


class QMediaPlaylistControl;
class QMediaPlaylistProvider;
class QMediaPlaylistReader;
class QMediaPlaylistWriter;
class QMediaPlayerControl;

class QMediaPlaylistPrivate
{
    Q_DECLARE_PUBLIC(QMediaPlaylist)
public:
    QMediaPlaylistPrivate()
        :mediaObject(nullptr),
        control(nullptr),
        networkPlaylistControl(nullptr),
        error(QMediaPlaylist::NoError)
    {
    }

    virtual ~QMediaPlaylistPrivate() {}

    void _q_loadFailed(QMediaPlaylist::Error error, const QString &errorString)
    {
        this->error = error;
        this->errorString = errorString;

        emit q_ptr->loadFailed();
    }

    void _q_mediaObjectDeleted()
    {
        Q_Q(QMediaPlaylist);
        mediaObject = nullptr;
        if (control != networkPlaylistControl)
            control = nullptr;
        q->setMediaObject(nullptr);
    }

    QMediaObject *mediaObject;

    QMediaPlaylistControl *control;
    QMediaPlaylistProvider *playlist() const { return control->playlistProvider(); }

    QMediaPlaylistControl *networkPlaylistControl;

    bool readItems(QMediaPlaylistReader *reader);
    bool writeItems(QMediaPlaylistWriter *writer);

    void syncControls(QMediaPlaylistControl *oldControl, QMediaPlaylistControl *newControl,
                      int *removedStart, int *removedEnd,
                      int *insertedStart, int *insertedEnd);

    QMediaPlaylist::Error error;
    QString errorString;

    QMediaPlaylist *q_ptr;
};


class QMediaNetworkPlaylistControl : public QMediaPlaylistControl
{
    Q_OBJECT
public:
    QMediaNetworkPlaylistControl(QObject *parent)
        :QMediaPlaylistControl(parent)
    {
        QMediaPlaylistProvider *playlist = new QMediaNetworkPlaylistProvider(this);
        m_navigator = new QMediaPlaylistNavigator(playlist,this);
        m_navigator->setPlaybackMode(QMediaPlaylist::Sequential);

        connect(m_navigator, SIGNAL(currentIndexChanged(int)), SIGNAL(currentIndexChanged(int)));
        connect(m_navigator, SIGNAL(activated(QMediaContent)), SIGNAL(currentMediaChanged(QMediaContent)));
        connect(m_navigator, SIGNAL(playbackModeChanged(QMediaPlaylist::PlaybackMode)), SIGNAL(playbackModeChanged(QMediaPlaylist::PlaybackMode)));
    }

    ~QMediaNetworkPlaylistControl() {}

    QMediaPlaylistProvider* playlistProvider() const override { return m_navigator->playlist(); }
    bool setPlaylistProvider(QMediaPlaylistProvider *mediaPlaylist) override
    {
        m_navigator->setPlaylist(mediaPlaylist);
        emit playlistProviderChanged();
        return true;
    }

    int currentIndex() const override { return m_navigator->currentIndex(); }
    void setCurrentIndex(int position) override { m_navigator->jump(position); }
    int nextIndex(int steps) const override { return m_navigator->nextIndex(steps); }
    int previousIndex(int steps) const override { return m_navigator->previousIndex(steps); }

    void next() override { m_navigator->next(); }
    void previous() override { m_navigator->previous(); }

    QMediaPlaylist::PlaybackMode playbackMode() const override { return m_navigator->playbackMode(); }
    void setPlaybackMode(QMediaPlaylist::PlaybackMode mode) override { m_navigator->setPlaybackMode(mode); }

private:
    QMediaPlaylistNavigator *m_navigator;
};


QT_END_NAMESPACE


#endif // QMEDIAPLAYLIST_P_H
