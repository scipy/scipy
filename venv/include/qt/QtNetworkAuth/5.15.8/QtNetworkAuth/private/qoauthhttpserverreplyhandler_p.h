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

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of the Network Access API.  This header file may change from
// version to version without notice, or even be removed.
//
// We mean it.
//

#ifndef QOAUTHHTTPSERVERREPLYHANDLER_P_H
#define QOAUTHHTTPSERVERREPLYHANDLER_P_H

#ifndef QT_NO_HTTP

#include <QtNetworkAuth/qoauthglobal.h>
#include <QtNetworkAuth/qoauthhttpserverreplyhandler.h>

#include <private/qobject_p.h>

#include <QtNetwork/qtcpserver.h>

QT_BEGIN_NAMESPACE

class QOAuthHttpServerReplyHandlerPrivate
{
    Q_DECLARE_PUBLIC(QOAuthHttpServerReplyHandler)

public:
    explicit QOAuthHttpServerReplyHandlerPrivate(QOAuthHttpServerReplyHandler *p);
    ~QOAuthHttpServerReplyHandlerPrivate();

    QTcpServer httpServer;
    QString text;
    QHostAddress listenAddress = QHostAddress::LocalHost;
    QString path;

private:
    void _q_clientConnected();
    void _q_readData(QTcpSocket *socket);
    void _q_answerClient(QTcpSocket *socket, const QUrl &url);

    struct QHttpRequest {
        quint16 port = 0;

        bool readMethod(QTcpSocket *socket);
        bool readUrl(QTcpSocket *socket);
        bool readStatus(QTcpSocket *socket);
        bool readHeader(QTcpSocket *socket);

        enum class State {
            ReadingMethod,
            ReadingUrl,
            ReadingStatus,
            ReadingHeader,
            ReadingBody,
            AllDone
        } state = State::ReadingMethod;
        QByteArray fragment;

        enum class Method {
            Unknown,
            Head,
            Get,
            Put,
            Post,
            Delete,
        } method = Method::Unknown;
        QUrl url;
        QPair<quint8, quint8> version;
        QMap<QByteArray, QByteArray> headers;
    };

    QMap<QTcpSocket *, QHttpRequest> clients;

    QOAuthHttpServerReplyHandler *q_ptr;
};

QT_END_NAMESPACE

#endif // QT_NO_HTTP

#endif // QOAUTHHTTPSERVERREPLYHANDLER_P_H
