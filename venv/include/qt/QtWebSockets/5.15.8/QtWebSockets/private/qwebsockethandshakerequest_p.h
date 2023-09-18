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

#ifndef QWEBSOCKETHANDSHAKEREQUEST_P_H
#define QWEBSOCKETHANDSHAKEREQUEST_P_H
//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <QtCore/QObject>
#include <QtCore/QMap>
#include <QtCore/QString>
#include <QtCore/QUrl>

#include "qwebsocketprotocol.h"

QT_BEGIN_NAMESPACE

class QTextStream;

class Q_AUTOTEST_EXPORT QWebSocketHandshakeRequest
{
    Q_DISABLE_COPY(QWebSocketHandshakeRequest)

public:
    QWebSocketHandshakeRequest(int port, bool isSecure);
    virtual ~QWebSocketHandshakeRequest();

    void clear();

    int port() const;
    bool isSecure() const;
    bool isValid() const;
    QMap<QString, QString> headers() const;
    QList<QWebSocketProtocol::Version> versions() const;
    QString key() const;
    QString origin() const;
    QList<QString> protocols() const;
    QList<QString> extensions() const;
    QUrl requestUrl() const;
    QString resourceName() const;
    QString host() const;

    void readHandshake(QTextStream &textStream, int maxHeaderLineLength, int maxHeaders);

private:

    int m_port;
    bool m_isSecure;
    bool m_isValid;
    QMultiMap<QString, QString> m_headers;
    QList<QWebSocketProtocol::Version> m_versions;
    QString m_key;
    QString m_origin;
    QList<QString> m_protocols;
    QList<QString> m_extensions;
    QUrl m_requestUrl;
};

QT_END_NAMESPACE

#endif // QWEBSOCKETHANDSHAKEREQUEST_P_H
