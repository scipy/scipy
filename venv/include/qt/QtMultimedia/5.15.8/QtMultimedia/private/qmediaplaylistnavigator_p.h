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

#ifndef QMEDIAPLAYLISTNAVIGATOR_P_H
#define QMEDIAPLAYLISTNAVIGATOR_P_H

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

#include "qmediaplaylistprovider_p.h"
#include "qmediaplaylist.h"
#include <QtCore/qobject.h>

QT_BEGIN_NAMESPACE


class QMediaPlaylistNavigatorPrivate;
class Q_MULTIMEDIA_EXPORT QMediaPlaylistNavigator : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QMediaPlaylist::PlaybackMode playbackMode READ playbackMode WRITE setPlaybackMode NOTIFY playbackModeChanged)
    Q_PROPERTY(int currentIndex READ currentIndex WRITE jump NOTIFY currentIndexChanged)
    Q_PROPERTY(QMediaContent currentItem READ currentItem)

public:
    QMediaPlaylistNavigator(QMediaPlaylistProvider *playlist, QObject *parent = nullptr);
    virtual ~QMediaPlaylistNavigator();

    QMediaPlaylistProvider *playlist() const;
    void setPlaylist(QMediaPlaylistProvider *playlist);

    QMediaPlaylist::PlaybackMode playbackMode() const;

    QMediaContent currentItem() const;
    QMediaContent nextItem(int steps = 1) const;
    QMediaContent previousItem(int steps = 1) const;

    QMediaContent itemAt(int position) const;

    int currentIndex() const;
    int nextIndex(int steps = 1) const;
    int previousIndex(int steps = 1) const;

public Q_SLOTS:
    void next();
    void previous();

    void jump(int);

    void setPlaybackMode(QMediaPlaylist::PlaybackMode mode);

Q_SIGNALS:
    void activated(const QMediaContent &content);
    void currentIndexChanged(int);
    void playbackModeChanged(QMediaPlaylist::PlaybackMode mode);

    void surroundingItemsChanged();

protected:
    QMediaPlaylistNavigatorPrivate *d_ptr;

private:
    Q_DISABLE_COPY(QMediaPlaylistNavigator)
    Q_DECLARE_PRIVATE(QMediaPlaylistNavigator)

    Q_PRIVATE_SLOT(d_func(), void _q_mediaInserted(int start, int end))
    Q_PRIVATE_SLOT(d_func(), void _q_mediaRemoved(int start, int end))
    Q_PRIVATE_SLOT(d_func(), void _q_mediaChanged(int start, int end))
};

QT_END_NAMESPACE


#endif // QMEDIAPLAYLISTNAVIGATOR_P_H
