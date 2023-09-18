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

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#ifndef Patternist_XsdXPathExpression_H
#define Patternist_XsdXPathExpression_H

#include <private/qanyuri_p.h>
#include <private/qnamedschemacomponent_p.h>
#include <private/qxsdannotated_p.h>

#include <QtCore/QList>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Represents a XSD assertion object.
     *
     * @ingroup Patternist_schema
     * @author Tobias Koenig <tobias.koenig@nokia.com>
     * @see <a href="http://www.w3.org/TR/xmlschema11-1/#x">XPathExpression Definition</a>
     */
    class XsdXPathExpression : public NamedSchemaComponent, public XsdAnnotated
    {
        public:
            typedef QExplicitlySharedDataPointer<XsdXPathExpression> Ptr;
            typedef QList<XsdXPathExpression::Ptr> List;

            /**
             * Sets the list of namespace @p bindings of the XPath expression.
             *
             * @see <a href="http://www.w3.org/TR/xmlschema11-1/#x-namespace_bindings">Namespace Bindings Definition</a>
             *
             * @note We can't use a QSet<QXmlName> here, as the hash method does not take the prefix
             *       in account, so we loose entries.
             */
            void setNamespaceBindings(const QList<QXmlName> &bindings);

            /**
             * Returns the list of namespace bindings of the XPath expression.
             */
            QList<QXmlName> namespaceBindings() const;

            /**
             * Sets the default namespace of the XPath expression.
             *
             * @see <a href="http://www.w3.org/TR/xmlschema11-1/#x-default_namespace">Default Namespace Definition</a>
             */
            void setDefaultNamespace(const AnyURI::Ptr &defaultNamespace);

            /**
             * Returns the default namespace of the XPath expression.
             */
            AnyURI::Ptr defaultNamespace() const;

            /**
             * Sets the base @p uri of the XPath expression.
             *
             * @see <a href="http://www.w3.org/TR/xmlschema11-1/#x-base_URI">Base URI Definition</a>
             */
            void setBaseURI(const AnyURI::Ptr &uri);

            /**
             * Returns the base uri of the XPath expression.
             */
            AnyURI::Ptr baseURI() const;

            /**
             * Sets the @p expression string of the XPath expression.
             *
             * @see <a href="http://www.w3.org/TR/xmlschema11-1/#x-expression">Expression Definition</a>
             */
            void setExpression(const QString &expression);

            /**
             * Returns the expression string of the XPath expression.
             */
            QString expression() const;

        private:
            QList<QXmlName> m_namespaceBindings;
            AnyURI::Ptr     m_defaultNamespace;
            AnyURI::Ptr     m_baseURI;
            QString         m_expression;
    };
}

QT_END_NAMESPACE

#endif
