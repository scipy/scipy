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

#ifndef Patternist_IndexOfIterator_H
#define Patternist_IndexOfIterator_H

#include <private/qitem_p.h>
#include <private/qatomiccomparator_p.h>
#include <private/qcomparisonplatform_p.h>
#include <private/qdynamiccontext_p.h>
#include <private/qexpression_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Forms an QAbstractXmlForwardIterator over a sequence of integers, which each is the position
     * of where a search parameter appeared in another QAbstractXmlForwardIterator.
     *
     * @see <a href="http://www.w3.org/TR/xpath-functions/#func-index-of">XQuery 1.0
     * and XPath 2.0 Functions and Operators, 15.1.3 fn:index-of</a>
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_iterators
     */
    class IndexOfIterator : public Item::Iterator
                          , public ComparisonPlatform<IndexOfIterator, false>
                          , public SourceLocationReflection
    {
    public:

        /**
         * Creates an IndexOfIterator, whose next() function returns integers being
         * the index positions of where @p searchParam was found in @p inputSequence.
         *
         * @param comp the AtomicComparator to be used for comparing values. This may be @c null,
         * meaning the IndexOfIterator iterator will dynamically determine what comparator to use
         * on an item per item basis, which is slower.
         * @param searchParam the item which should be compared to the items in @p inputSequence.
         * @param inputSequence the input sequence which indexes of the @p searchParam should
         * be returned for.
         * @param context the usual DynamicContext
         * @param expr the Expression that this IndexOfIterator is evaluating
         * for. It is used for error reporting, via actualReflection().
         */
        IndexOfIterator(const Item::Iterator::Ptr &inputSequence,
                        const Item &searchParam,
                        const AtomicComparator::Ptr &comp,
                        const DynamicContext::Ptr &context,
                        const Expression::ConstPtr &expr);

        virtual Item next();
        virtual Item current() const;
        virtual xsInteger position() const;
        virtual Item::Iterator::Ptr copy() const;

        inline AtomicComparator::Operator operatorID() const
        {
            return AtomicComparator::OperatorEqual;
        }

        virtual const SourceLocationReflection *actualReflection() const;

    private:
        const Item::Iterator::Ptr   m_seq;
        const Item                  m_searchParam;
        const DynamicContext::Ptr   m_context;
        const Expression::ConstPtr  m_expr;
        Item                        m_current;
        xsInteger                   m_position;
        xsInteger                   m_seqPos;
    };
}

QT_END_NAMESPACE

#endif
