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

#ifndef QOAUTH1SIGNATURE_P_H
#define QOAUTH1SIGNATURE_P_H

#include <QtNetworkAuth/qoauth1signature.h>

#include <QtCore/qurl.h>
#include <QtCore/qstring.h>
#include <QtCore/qvariant.h>
#include <QtCore/qshareddata.h>

QT_BEGIN_NAMESPACE

class QOAuth1SignaturePrivate : public QSharedData
{
public:
    QOAuth1SignaturePrivate() = default;
    QOAuth1SignaturePrivate(const QUrl &url, QOAuth1Signature::HttpRequestMethod method,
                            const QVariantMap &parameters,
                            const QString &clientSharedKey = QString(),
                            const QString &tokenSecret = QString());

    QByteArray signatureBaseString() const;
    QByteArray secret() const;
    static QByteArray parameterString(const QVariantMap &parameters);
    static QByteArray encodeHeaders(const QVariantMap &headers);


    QOAuth1Signature::HttpRequestMethod method = QOAuth1Signature::HttpRequestMethod::Post;
    QByteArray customVerb;
    QUrl url;
    QString clientSharedKey;
    QString tokenSecret;
    QMultiMap<QString, QVariant> parameters;

    static QOAuth1SignaturePrivate shared_null;
};

QT_END_NAMESPACE

#endif // QOAUTH1SIGNATURE_P_H
