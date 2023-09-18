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

#ifndef QOAUTH1_P_H
#define QOAUTH1_P_H

#ifndef QT_NO_HTTP

#include <private/qabstractoauth_p.h>

#include <QtNetworkAuth/qoauth1.h>
#include <QtNetworkAuth/qoauthglobal.h>

#include <QtCore/qurl.h>
#include <QtCore/qpair.h>
#include <QtCore/qobject.h>
#include <QtCore/qstring.h>

#include <QtNetwork/qnetworkreply.h>
#include <QtNetwork/qnetworkaccessmanager.h>

QT_BEGIN_NAMESPACE

class QOAuth1Signature;

class QOAuth1Private : public QAbstractOAuthPrivate
{
    Q_DECLARE_PUBLIC(QOAuth1)

public:
    QOAuth1Private(const QPair<QString, QString> &clientCredentials,
                   QNetworkAccessManager *networkAccessManager = nullptr);

    void appendCommonHeaders(QVariantMap *headers);
    void appendSignature(QAbstractOAuth::Stage stage,
                         QVariantMap *headers,
                         const QUrl &url,
                         QNetworkAccessManager::Operation operation,
                         const QVariantMap parameters);

    QNetworkReply *requestToken(QNetworkAccessManager::Operation operation,
                                const QUrl &url,
                                const QPair<QString, QString> &token,
                                const QVariantMap &additionalParameters);

    QString signatureMethodString() const;
    QByteArray generateSignature(const QVariantMap &parameters,
                                 const QUrl &url,
                                 QNetworkAccessManager::Operation operation) const;
    QByteArray generateSignature(const QVariantMap &parameters,
                                 const QUrl &url,
                                 const QByteArray &verb) const;
    QByteArray formatSignature(const QOAuth1Signature &signature) const;

    QVariantMap createOAuthBaseParams() const;

    void prepareRequestImpl(QNetworkRequest *request,
                            const QByteArray &verb,
                            const QByteArray &body) override;

    void _q_onTokenRequestError(QNetworkReply::NetworkError error);
    void _q_tokensReceived(const QVariantMap &tokens);

    QString clientIdentifierSharedKey;
    QString tokenSecret;
    QString verifier;
    QUrl temporaryCredentialsUrl;
    QUrl tokenCredentialsUrl;
    QOAuth1::SignatureMethod signatureMethod = QOAuth1::SignatureMethod::Hmac_Sha1;
    const QString oauthVersion = QStringLiteral("1.0");
    bool tokenRequested = false;

    struct OAuth1KeyString
    {
        static const QString oauthCallback;
        static const QString oauthCallbackConfirmed;
        static const QString oauthConsumerKey;
        static const QString oauthNonce;
        static const QString oauthSignature;
        static const QString oauthSignatureMethod;
        static const QString oauthTimestamp;
        static const QString oauthToken;
        static const QString oauthTokenSecret;
        static const QString oauthVerifier;
        static const QString oauthVersion;
    };
};

QT_END_NAMESPACE

#endif // QT_NO_HTTP

#endif // QOAUTH1_P_H
