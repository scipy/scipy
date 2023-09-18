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

#ifndef Patternist_ContextItem_H
#define Patternist_ContextItem_H

#include <private/qemptycontainer_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Implements the context item, the dot: <tt>.</tt>.
     *
     * @see <a href="http://www.w3.org/TR/xpath20/#id-context-item-expression">XML Path Language
     * (XPath) 2.0, 3.1.4 Context Item Expression</a>
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_expressions
     */
    class ContextItem : public EmptyContainer
    {
    public:
        /**
         * @p expr is possibly used for error reporting. If this context item has been
         * created implicitly, such as for the expression <tt>fn:string()</tt>, @p expr
         * should be passed a valid pointer to the Expression that this context
         * item is generated for.
         */
        inline ContextItem(const Expression::Ptr &expr = Expression::Ptr()) : m_expr(expr)
        {
        }

        virtual Item evaluateSingleton(const DynamicContext::Ptr &context) const;
        virtual SequenceType::Ptr staticType() const;

        /**
         * @returns always DisableElimination and RequiresContextItem
         */
        virtual Expression::Properties properties() const;

        virtual ExpressionVisitorResult::Ptr accept(const ExpressionVisitor::Ptr &visitor) const;

        /**
         * Overridden to store a pointer to StaticContext::contextItemType().
         */
        virtual Expression::Ptr compress(const StaticContext::Ptr &context);

        /**
         * Overridden to store a pointer to StaticContext::contextItemType().
         */
        virtual Expression::Ptr typeCheck(const StaticContext::Ptr &context,
                                          const SequenceType::Ptr &reqType);

        /**
         * @returns always IDContextItem
         */
        virtual ID id() const;

        /**
         * @returns always BuiltinTypes::item;
         */
        virtual ItemType::Ptr expectedContextItemType() const;

        virtual const SourceLocationReflection *actualReflection() const;
        virtual void announceFocusType(const ItemType::Ptr &type);

    private:
        ItemType::Ptr           m_itemType;
        const Expression::Ptr   m_expr;
    };
}

QT_END_NAMESPACE

#endif
