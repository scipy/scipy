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

#ifndef Patternist_LiteralSequence_H
#define Patternist_LiteralSequence_H

#include <private/qemptycontainer_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Houses a sequence of atomic values, making it available as an Expression.
     *
     * This is not only literals that can be created via the XQuery syntax(strings and numbers), but
     * all other atomic values, such as <tt>xs:date</tt> or <tt>xs:time</tt>. It is not guaranteed
     * that consecutive atomic values are represented in a LiteralSequence.
     *
     * @see <a href="http://www.w3.org/TR/xquery/#id-literals">XQuery 1.0: An XML Query Language,
     * 3.1.1 Literals</a>
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_expressions
     */
    class LiteralSequence : public EmptyContainer
    {
    public:
        /**
         * Creates a LiteralSequence that represents @p item.
         *
         * @param list the list of item. No entry may be @c null. The list
         * must at least be two entries large.
         */
        LiteralSequence(const Item::List &list);

        virtual Item::Iterator::Ptr evaluateSequence(const DynamicContext::Ptr &context) const;
        virtual SequenceType::Ptr staticType() const;

        virtual ExpressionVisitorResult::Ptr accept(const ExpressionVisitor::Ptr &visitor) const;

        virtual ID id() const;

        virtual Properties properties() const;

    private:
        const Item::List m_list;
    };
}

QT_END_NAMESPACE

#endif
