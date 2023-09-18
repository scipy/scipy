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

#ifndef Patternist_GenericPredicate_H
#define Patternist_GenericPredicate_H

#include <private/qpaircontainer_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short A predicate that can handle all kinds of predicates and
     * is therefore not very efficient, but can cope with all the tricky scenarios.
     *
     * @see FirstItemPredicate
     * @see TruthPredicate
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_expressions
     */
    class GenericPredicate : public PairContainer
    {
    public:

        /**
         * Creates a predicate expression that filters the items gained
         * from evaluating @p sourceExpression through the filter @p predicateExpression.
         *
         * This function performs type analyzis on the passed expressions, and may
         * return more specialized expressions depending on the analyzis.
         *
         * If @p predicateExpression is an invalid predicate, an error is issued
         * via the @p context.
         */
        static Expression::Ptr create(const Expression::Ptr &sourceExpression,
                                      const Expression::Ptr &predicateExpression,
                                      const StaticContext::Ptr &context,
                                      const QSourceLocation &location);

        static Expression::Ptr createFirstItem(const Expression::Ptr &sourceExpression);

        /**
         * Creates a source iterator which is passed to the ItemMappingIterator
         * and the Focus. The ItemMappingIterator modifies it with
         * its QAbstractXmlForwardIterator::next() calls, and since the Focus references the same QAbstractXmlForwardIterator,
         * the focus is automatically moved.
         */
        virtual Item::Iterator::Ptr evaluateSequence(const DynamicContext::Ptr &context) const;

        /**
         * Doesn't return the first item from calling evaluateSequence(), but does the mapping
         * manually. This avoid allocating an ItemMappingIterator.
         */
        virtual Item evaluateSingleton(const DynamicContext::Ptr &context) const;

        inline Item mapToItem(const Item &subject,
                                   const DynamicContext::Ptr &) const;

        virtual SequenceType::List expectedOperandTypes() const;
        virtual SequenceType::Ptr staticType() const;
        virtual ExpressionVisitorResult::Ptr accept(const ExpressionVisitor::Ptr &visitor) const;
        virtual ID id() const;

        /**
         * @returns always CreatesFocusForLast.
         */
        virtual Properties properties() const;

        virtual QString description() const;

    protected:

        /**
         * Creates a GenericPredicate which filters the items from the @p sourceExpression
         * through @p predicate.
         *
         * This constructor is protected. The proper way to create predicates is via the static
         * create() function.
         */
        GenericPredicate(const Expression::Ptr &sourceExpression,
                         const Expression::Ptr &predicate);

        /**
         * @returns the ItemType of the first operand's staticType().
         */
        virtual ItemType::Ptr newFocusType() const;

    private:
        typedef QExplicitlySharedDataPointer<const GenericPredicate> ConstPtr;
    };
}

QT_END_NAMESPACE

#endif
