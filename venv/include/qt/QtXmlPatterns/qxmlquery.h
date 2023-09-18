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

#ifndef QXMLQUERY_H
#define QXMLQUERY_H

#include <QtCore/QUrl>
#include <QtXmlPatterns/QAbstractXmlNodeModel>
#include <QtXmlPatterns/QAbstractXmlReceiver>
#include <QtXmlPatterns/QXmlNamePool>

QT_BEGIN_NAMESPACE


class QAbstractMessageHandler;
class QAbstractUriResolver;
class QIODevice;
class QNetworkAccessManager;
class QXmlName;
class QXmlNodeIndex;
class QXmlQueryPrivate;
class QXmlResultItems;
class QXmlSerializer;

/* The members in the namespace QPatternistSDK are internal, not part of the public API, and
 * unsupported. Using them leads to undefined behavior. */
namespace QPatternistSDK
{
    class TestCase;
}

namespace QPatternist
{
    class XsdSchemaParser;
    class XsdValidatingInstanceReader;
    class VariableLoader;
}

class Q_XMLPATTERNS_EXPORT QXmlQuery
{
public:
    enum QueryLanguage
    {
        XQuery10                                = 1,
        XSLT20                                  = 2,
        XmlSchema11IdentityConstraintSelector   = 1024,
        XmlSchema11IdentityConstraintField      = 2048,
        XPath20                                 = 4096
    };

    QXmlQuery();
    QXmlQuery(const QXmlQuery &other);
    QXmlQuery(const QXmlNamePool &np);
    QXmlQuery(QueryLanguage queryLanguage,
              const QXmlNamePool &np = QXmlNamePool());
    ~QXmlQuery();
    QXmlQuery &operator=(const QXmlQuery &other);

    void setMessageHandler(QAbstractMessageHandler *messageHandler);
    QAbstractMessageHandler *messageHandler() const;

    void setQuery(const QString &sourceCode, const QUrl &documentURI = QUrl());
    void setQuery(QIODevice *sourceCode, const QUrl &documentURI = QUrl());
    void setQuery(const QUrl &queryURI, const QUrl &baseURI = QUrl());

    QXmlNamePool namePool() const;

    void bindVariable(const QXmlName &name, const QXmlItem &value);
    void bindVariable(const QString &localName, const QXmlItem &value);

    void bindVariable(const QXmlName &name, QIODevice *);
    void bindVariable(const QString &localName, QIODevice *);
    void bindVariable(const QXmlName &name, const QXmlQuery &query);
    void bindVariable(const QString &localName, const QXmlQuery &query);

    bool isValid() const;

    void evaluateTo(QXmlResultItems *result) const;
    bool evaluateTo(QAbstractXmlReceiver *callback) const;
    bool evaluateTo(QStringList *target) const;
    bool evaluateTo(QIODevice *target) const;
    bool evaluateTo(QString *output) const;

    void setUriResolver(const QAbstractUriResolver *resolver);
    const QAbstractUriResolver *uriResolver() const;

    void setFocus(const QXmlItem &item);
    bool setFocus(const QUrl &documentURI);
    bool setFocus(QIODevice *document);
    bool setFocus(const QString &focus);

    void setInitialTemplateName(const QXmlName &name);
    void setInitialTemplateName(const QString &name);
    QXmlName initialTemplateName() const;

    void setNetworkAccessManager(QNetworkAccessManager *newManager);
    QNetworkAccessManager *networkAccessManager() const;

    QueryLanguage queryLanguage() const;
private:
    friend class QXmlName;
    friend class QXmlSerializer;
    friend class QPatternistSDK::TestCase;
    friend class QPatternist::XsdSchemaParser;
    friend class QPatternist::XsdValidatingInstanceReader;
    friend class QPatternist::VariableLoader;
    template<typename TInputType> friend bool setFocusHelper(QXmlQuery *const queryInstance,
                                                             const TInputType &focusValue);
    QXmlQueryPrivate *d;
};

QT_END_NAMESPACE

#endif
