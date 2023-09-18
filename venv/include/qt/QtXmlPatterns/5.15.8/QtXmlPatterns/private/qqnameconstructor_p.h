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

#ifndef Patternist_QNameConstructor_H
#define Patternist_QNameConstructor_H

#include <private/qsinglecontainer_p.h>
#include <private/qbuiltintypes_p.h>
#include <private/qpatternistlocale_p.h>
#include <private/qxpathhelper_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Creates an @c xs:QName value from a lexical QName using
     * statically known namespace bindings.
     *
     * @see QQNameValue
     * @see QXmlUtils
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_expressions
     */
    class QNameConstructor : public SingleContainer
    {
    public:

        QNameConstructor(const Expression::Ptr &source,
                         const NamespaceResolver::Ptr &nsResolver);

        virtual Item evaluateSingleton(const DynamicContext::Ptr &) const;

        virtual SequenceType::List expectedOperandTypes() const;

        virtual SequenceType::Ptr staticType() const;

        virtual ExpressionVisitorResult::Ptr accept(const ExpressionVisitor::Ptr &visitor) const;

        /**
         * Expands @p lexicalQName, which is a lexical representation of a QName such as "x:body", into
         * a QName using @p nsResolver to supply the namespace bindings.
         *
         * If @p lexicalQName is lexically invalid @p InvalidQName is raised via @p context, or if
         * no namespace binding does not exists for a prefix(if any) in @p lexicalQName, @p NoBinding
         * is raised via @p context.
         *
         * If @p asForAttribute is @c true, the name is considered to be for an
         * attribute in some way, and @p lexicalQName will not pick up the
         * default namespace if it doesn't have a prefix.
         *
         * @p nsResolver is parameterized meaning the function can be instantiated with either
         * DynamicContext or StaticContext.
         *
         * @see QQNameValue
         * @see QXmlUtils
         */
        template<typename TReportContext,
                 const ReportContext::ErrorCode InvalidQName,
                 const ReportContext::ErrorCode NoBinding>
        static
        QXmlName expandQName(const QString &lexicalQName,
                             const TReportContext &context,
                             const NamespaceResolver::Ptr &nsResolver,
                             const SourceLocationReflection *const r,
                             const bool asForAttribute = false);

        /**
         * Resolves the namespace prefix @p prefix to its namespace if it exists, or
         * raised ReportContext::XPST0081 otherwise.
         *
         * @returns the namespace URI corresponding to @p prefix
         */
        static QXmlName::NamespaceCode namespaceForPrefix(const QXmlName::PrefixCode prefix,
                                                          const StaticContext::Ptr &context,
                                                          const SourceLocationReflection *const r);

        virtual const SourceLocationReflection *actualReflection() const;

    private:
        const NamespaceResolver::Ptr m_nsResolver;
    };

    template<typename TReportContext,
             const ReportContext::ErrorCode InvalidQName,
             const ReportContext::ErrorCode NoBinding>
    QXmlName QNameConstructor::expandQName(const QString &lexicalQName,
                                           const TReportContext &context,
                                           const NamespaceResolver::Ptr &nsResolver,
                                           const SourceLocationReflection *const r,
                                           const bool asForAttribute)
    {
        Q_ASSERT(nsResolver);
        Q_ASSERT(context);

        if(XPathHelper::isQName(lexicalQName))
        {
            QString prefix;
            QString local;
            XPathHelper::splitQName(lexicalQName, prefix, local);
            const QXmlName::NamespaceCode nsCode = asForAttribute && prefix.isEmpty() ? QXmlName::NamespaceCode(StandardNamespaces::empty)
                                                                                   : (nsResolver->lookupNamespaceURI(context->namePool()->allocatePrefix(prefix)));

            if(nsCode == NamespaceResolver::NoBinding)
            {
                context->error(QtXmlPatterns::tr("No namespace binding exists for "
                                  "the prefix %1 in %2").arg(formatKeyword(prefix),
                                                             formatKeyword(lexicalQName)),
                               NoBinding,
                               r);
                return QXmlName(); /* Silence compiler warning. */
            }
            else
                return context->namePool()->allocateQName(context->namePool()->stringForNamespace(nsCode), local, prefix);
        }
        else
        {
            context->error(QtXmlPatterns::tr("%1 is an invalid %2")
                              .arg(formatData(lexicalQName))
                              .arg(formatType(context->namePool(), BuiltinTypes::xsQName)),
                           InvalidQName,
                           r);
            return QXmlName(); /* Silence compiler warning. */
        }
    }
}

QT_END_NAMESPACE

#endif
