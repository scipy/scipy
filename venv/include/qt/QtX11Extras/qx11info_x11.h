/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtGui module of the Qt Toolkit.
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

#ifndef QX11INFO_X11_H
#define QX11INFO_X11_H

#include <QtCore/qnamespace.h>
#include "QtX11Extras/qtx11extrasglobal.h"

#include <xcb/xcb.h>

typedef struct _XDisplay Display;

QT_BEGIN_NAMESPACE

class Q_X11EXTRAS_EXPORT QX11Info
{
public:
    enum PeekOption {
        PeekDefault = 0,
        PeekFromCachedIndex = 1
    };
    Q_DECLARE_FLAGS(PeekOptions, PeekOption)

    static bool isPlatformX11();

    static int appDpiX(int screen=-1);
    static int appDpiY(int screen=-1);

#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    static unsigned long appRootWindow(int screen=-1);
#else
    static quint32 appRootWindow(int screen=-1);
#endif
    static int appScreen();

    static unsigned long appTime();
    static unsigned long appUserTime();

    static void setAppTime(unsigned long time);
    static void setAppUserTime(unsigned long time);

    static unsigned long getTimestamp();

    static QByteArray nextStartupId();
    static void setNextStartupId(const QByteArray &id);

    static Display *display();
    static xcb_connection_t *connection();

    static bool isCompositingManagerRunning(int screen = -1);

    static qint32 generatePeekerId();
    static bool removePeekerId(qint32 peekerId);
    typedef bool (*PeekerCallback)(xcb_generic_event_t *event, void *peekerData);
    static bool peekEventQueue(PeekerCallback peeker, void *peekerData = nullptr,
                               PeekOptions option = PeekDefault, qint32 peekerId = -1);

private:
    QX11Info();
};

Q_DECLARE_OPERATORS_FOR_FLAGS(QX11Info::PeekOptions)

QT_END_NAMESPACE

#endif // QX11INFO_X11_H

