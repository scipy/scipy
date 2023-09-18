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

#ifndef QABSTRACTOAUTH2_P_H
#define QABSTRACTOAUTH2_P_H

#ifndef QT_NO_HTTP

#include <private/qabstractoauth_p.h>

#include <QtNetworkAuth/qoauthglobal.h>
#include <QtNetworkAuth/qabstractoauth2.h>

#include <QtCore/qurl.h>
#include <QtCore/qstring.h>
#include <QtCore/qpointer.h>

#include <QtNetwork/qnetworkreply.h>

QT_BEGIN_NAMESPACE

class QNetworkAccessManager;

class QAbstractOAuth2Private : public QAbstractOAuthPrivate
{
    Q_DECLARE_PUBLIC(QAbstractOAuth2)

public:
    QAbstractOAuth2Private(const QPair<QString, QString> &clientCredentials,
                           const QUrl &authorizationUrl, QNetworkAccessManager *manager = nullptr);
    ~QAbstractOAuth2Private();

    static QString generateRandomState();
    QNetworkRequest createRequest(QUrl url, const QVariantMap *parameters = nullptr);

    void prepareRequestImpl(QNetworkRequest *request,
                            const QByteArray &verb,
                            const QByteArray &body) override;

    QString clientIdentifierSharedKey;
    QString scope;
    QString state = generateRandomState();
    QString userAgent = QStringLiteral("QtOAuth/1.0 (+https://www.qt.io)");
    QString responseType;
    const QString bearerFormat = QStringLiteral("Bearer %1"); // Case sensitive
    QDateTime expiresAt;
    QString refreshToken;

    struct OAuth2KeyString
    {
        static const QString accessToken;
        static const QString apiKey;
        static const QString clientIdentifier;
        static const QString clientSharedSecret;
        static const QString code;
        static const QString error;
        static const QString errorDescription;
        static const QString errorUri;
        static const QString expiresIn;
        static const QString grantType;
        static const QString redirectUri;
        static const QString refreshToken;
        static const QString responseType;
        static const QString scope;
        static const QString state;
        static const QString tokenType;
    };
};

QT_END_NAMESPACE

#endif // QT_NO_HTTP

#endif // QABSTRACTOAUTH2_P_H
