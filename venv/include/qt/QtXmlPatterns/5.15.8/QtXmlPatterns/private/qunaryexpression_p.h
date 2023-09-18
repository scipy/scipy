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

#ifndef Patternist_UnaryExpression_H
#define Patternist_UnaryExpression_H

QT_BEGIN_NAMESPACE

#include <private/qarithmeticexpression_p.h>

namespace QPatternist
{
    /**
     * @short Implements XPath 2.0 unary expression, <tt>(-|+)ValueExpr</tt>.
     *
     * UnaryExpression is implemented by rewriting the expression <tt>operator [expr]</tt>
     * to the ArithmeticExpression <tt>0 operator [expr]</tt>. For example, the expression
     * <tt>+3</tt> becomes <tt>0 + 3</tt>, and <tt>-nodetest</tt> becomes <tt>0 - nodetest</tt>.
     *
     * On top of that expression ArithmeticExpression does the usual type
     * checking conversion. The only thing this class do, is to overide
     * evaluateSingleton() and calls Numeric::toNegated(). The reason this
     * UnaryExpression is needed at all and that <tt>0 - [expr]</tt> is
     * insufficent is that <tt>0 - xs:double(0)</tt> needs to return -0,
     * instead of 0. I know no other difference.
     *
     * In most cases the constant propagation optimization rewrites UnaryExpression into
     * a value, an instance of a sub-class of the Numeric class, wrapped with
     * Literal.
     *
     * Beyond the mathematical implication the unary expression have, it also
     * have the significant effect that it may invoke type promotion or that an expression
     * may contain a type error. For example, the expression "+'a string'" contains a type error, since
     * no unary operator is defined for @c xs:string. This is the reason why the '+' unary
     * operator isn't ignored.
     *
     * @see <a href="http://www.w3.org/TR/xpath20/#id-arithmetic">XML Path Language
     * (XPath) 2.0, 3.4 Arithmetic Expressions</a>
     * @see <a href="http://www.w3.org/TR/xpath-functions/#func-numeric-unary-plus">XQuery 1.0 and XPath
     * 2.0 Functions and Operators, 6.2.7 op:numeric-unary-plus</a>
     * @see <a href="http://www.w3.org/TR/xpath-functions/#func-numeric-unary-minus">XQuery 1.0 and XPath
     * 2.0 Functions and Operators, 6.2.8 op:numeric-unary-minus</a>
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_expressions
     */
    class UnaryExpression : public ArithmeticExpression
    {
    public:
        UnaryExpression(const AtomicMathematician::Operator op,
                        const Expression::Ptr &operand,
                        const StaticContext::Ptr &context);

        virtual Item evaluateSingleton(const DynamicContext::Ptr &context) const;

    private:
        Q_DISABLE_COPY(UnaryExpression)
    };
}

QT_END_NAMESPACE

#endif
