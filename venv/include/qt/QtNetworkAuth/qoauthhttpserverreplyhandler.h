/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Network Auth module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QOAUTHHTTPSERVERREPLYHANDLER_H
#define QOAUTHHTTPSERVERREPLYHANDLER_H

#ifndef QT_NO_HTTP

#include <QtNetworkAuth/qoauthglobal.h>
#include <QtNetworkAuth/qoauthoobreplyhandler.h>

#include <QtNetwork/qhostaddress.h>

QT_BEGIN_NAMESPACE

class QUrlQuery;

class QOAuthHttpServerReplyHandlerPrivate;
class Q_OAUTH_EXPORT QOAuthHttpServerReplyHandler : public QOAuthOobReplyHandler
{
    Q_OBJECT

public:
    explicit QOAuthHttpServerReplyHandler(QObject *parent = nullptr);
    explicit QOAuthHttpServerReplyHandler(quint16 port, QObject *parent = nullptr);
    explicit QOAuthHttpServerReplyHandler(const QHostAddress &address, quint16 port,
                                          QObject *parent = nullptr);
    ~QOAuthHttpServerReplyHandler();

    QString callback() const override;

    QString callbackPath() const;
    void setCallbackPath(const QString &path);

    QString callbackText() const;
    void setCallbackText(const QString &text);

    quint16 port() const;

    bool listen(const QHostAddress &address = QHostAddress::Any, quint16 port = 0);
    void close();
    bool isListening() const;

private:
    Q_DECLARE_PRIVATE(QOAuthHttpServerReplyHandler)
    QScopedPointer<QOAuthHttpServerReplyHandlerPrivate> d_ptr;
};

QT_END_NAMESPACE

#endif // QT_NO_HTTP

#endif // QOAUTHHTTPSERVERREPLYHANDLER_H
