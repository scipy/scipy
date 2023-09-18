/****************************************************************************
**
** Copyright (C) 2016 Kurt Pattyn <pattyn.kurt@gmail.com>.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtWebSockets module of the Qt Toolkit.
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
#ifndef QWEBSOCKETCORSAUTHENTICATOR_H
#define QWEBSOCKETCORSAUTHENTICATOR_H

#include "QtWebSockets/qwebsockets_global.h"
#include <QtCore/QScopedPointer>

QT_BEGIN_NAMESPACE

class QWebSocketCorsAuthenticatorPrivate;

class Q_WEBSOCKETS_EXPORT QWebSocketCorsAuthenticator
{
    Q_DECLARE_PRIVATE(QWebSocketCorsAuthenticator)

public:
    explicit QWebSocketCorsAuthenticator(const QString &origin);
    ~QWebSocketCorsAuthenticator();
    explicit QWebSocketCorsAuthenticator(const QWebSocketCorsAuthenticator &other);

#ifdef Q_COMPILER_RVALUE_REFS
    QWebSocketCorsAuthenticator(QWebSocketCorsAuthenticator &&other);
    QWebSocketCorsAuthenticator &operator =(QWebSocketCorsAuthenticator &&other);
#endif

    void swap(QWebSocketCorsAuthenticator &other);

    QWebSocketCorsAuthenticator &operator =(const QWebSocketCorsAuthenticator &other);

    QString origin() const;

    void setAllowed(bool allowed);
    bool allowed() const;

private:
    QScopedPointer<QWebSocketCorsAuthenticatorPrivate> d_ptr;
};

QT_END_NAMESPACE

#endif // QWEBSOCKETCORSAUTHENTICATOR_H
