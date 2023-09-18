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

#ifndef Patternist_XsdSchemaContext_H
#define Patternist_XsdSchemaContext_H

#include <private/qnamedschemacomponent_p.h>
#include <private/qreportcontext_p.h>
#include <private/qschematypefactory_p.h>
#include <private/qxsdschematoken_p.h>
#include <private/qxsdschema_p.h>
#include <private/qxsdschemachecker_p.h>
#include <private/qxsdschemaresolver_p.h>

#include <QtCore/QUrl>
#include <QtNetwork/QNetworkAccessManager>
#include <QtXmlPatterns/QAbstractMessageHandler>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short A context for schema parsing and validation.
     *
     * This class provides the infrastructure for error reporting and
     * network access. Additionally it stores objects that are used by
     * both, the parser and the validator.
     *
     * @ingroup Patternist_schema
     * @author Tobias Koenig <tobias.koenig@nokia.com>
     */
    class XsdSchemaContext : public ReportContext
    {
        public:
            /**
             * A smart pointer wrapping XsdSchemaContext instances.
             */
            typedef QExplicitlySharedDataPointer<XsdSchemaContext> Ptr;

            /**
             * Creates a new schema context object.
             *
             * @param namePool The name pool all names belong to.
             */
            XsdSchemaContext(const NamePool::Ptr &namePool);

            /**
             * Returns the name pool of the schema context.
             */
            virtual NamePool::Ptr namePool() const;

            /**
             * Sets the base URI for the main schema.
             *
             * The main schema is the one that includes resp. imports
             * all the other schema files.
             */
            virtual void setBaseURI(const QUrl &uri);

            /**
             * Returns the base URI of the main schema.
             */
            virtual QUrl baseURI() const;

            /**
             * Sets the network access manager that should be used
             * to access referenced schema definitions.
             */
            void setNetworkAccessManager(QNetworkAccessManager *accessManager);

            /**
             * Returns the network access manager that is used to
             * access referenced schema definitions.
             */
            virtual QNetworkAccessManager* networkAccessManager() const;

            /**
             * Sets the message @p handler used by the context for error reporting.
             */
            void setMessageHandler(QAbstractMessageHandler *handler);

            /**
             * Returns the message handler used by the context for
             * error reporting.
             */
            virtual QAbstractMessageHandler* messageHandler() const;

            /**
             * Always returns an empty source location.
             */
            virtual QSourceLocation locationFor(const SourceLocationReflection *const reflection) const;

            /**
             * Sets the uri @p resolver that is used for resolving URIs in the
             * schema parser.
             */
            void setUriResolver(const QAbstractUriResolver *resolver);

            /**
             * Returns the uri resolver that is used for resolving URIs in the
             * schema parser.
             */
            virtual const QAbstractUriResolver* uriResolver() const;

            /**
             * Returns the list of facets for the given simple @p type.
             */
            XsdFacet::Hash facetsForType(const AnySimpleType::Ptr &type) const;

            /**
             * Returns a schema type factory that contains some predefined schema types.
             */
            SchemaTypeFactory::Ptr schemaTypeFactory() const;

            /**
             * The following variables should not be accessed directly.
             */
            mutable SchemaTypeFactory::Ptr                 m_schemaTypeFactory;
            mutable QHash<SchemaType::Ptr, XsdFacet::Hash> m_builtinTypesFacetList;

        private:
            QHash<SchemaType::Ptr, XsdFacet::Hash> setupBuiltinTypesFacetList() const;

            NamePool::Ptr                                 m_namePool;
            QNetworkAccessManager*                        m_networkAccessManager;
            QUrl                                          m_baseURI;
            const QAbstractUriResolver*                   m_uriResolver;
            QAbstractMessageHandler*                      m_messageHandler;
    };
}

QT_END_NAMESPACE

#endif
