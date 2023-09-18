/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtXmlPatterns module of the Qt Toolkit.
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

#ifndef QXMLSCHEMAVALIDATOR_H
#define QXMLSCHEMAVALIDATOR_H

#include <QtCore/QUrl>
#include <QtXmlPatterns/QXmlNamePool>

QT_BEGIN_NAMESPACE


class QAbstractMessageHandler;
class QAbstractUriResolver;
class QIODevice;
class QNetworkAccessManager;
class QUrl;
class QXmlNamePool;
class QXmlSchema;
class QXmlSchemaValidatorPrivate;

class Q_XMLPATTERNS_EXPORT QXmlSchemaValidator
{
    public:
        QXmlSchemaValidator();
        QXmlSchemaValidator(const QXmlSchema &schema);
        ~QXmlSchemaValidator();

        void setSchema(const QXmlSchema &schema);

        bool validate(const QUrl &source) const;
        bool validate(QIODevice *source, const QUrl &documentUri = QUrl()) const;
        bool validate(const QByteArray &data, const QUrl &documentUri = QUrl()) const;

        QXmlNamePool namePool() const;
        QXmlSchema schema() const;

        void setMessageHandler(QAbstractMessageHandler *handler);
        QAbstractMessageHandler *messageHandler() const;

        void setUriResolver(const QAbstractUriResolver *resolver);
        const QAbstractUriResolver *uriResolver() const;

        void setNetworkAccessManager(QNetworkAccessManager *networkmanager);
        QNetworkAccessManager *networkAccessManager() const;

    private:
        QXmlSchemaValidatorPrivate* const d;

        Q_DISABLE_COPY(QXmlSchemaValidator)
};

QT_END_NAMESPACE

#endif
