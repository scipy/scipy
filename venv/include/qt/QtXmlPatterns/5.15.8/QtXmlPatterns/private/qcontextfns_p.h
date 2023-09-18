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

#ifndef Patternist_ContextFNs_H
#define Patternist_ContextFNs_H

#include <private/qfunctioncall_p.h>

/**
 * @file
 * @short Contains classes implementing the functions found in
 * <a href="http://www.w3.org/TR/xpath-functions/#context">XQuery 1.0 and
 * XPath 2.0 Functions and Operators, 16 Context Functions</a>.
 *
 * @ingroup Patternist_functions
 */

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Implements the function <tt>fn:position()</tt>.
     *
     * @see <a href="http://www.w3.org/TR/xpath-functions/#func-position">XQuery 1.0
     * and XPath 2.0 Functions and Operators, 16.1 fn:position</a>
     * @ingroup Patternist_functions
     * @author Frans Englich <frans.englich@nokia.com>
     */
    class PositionFN : public FunctionCall
    {
    public:
        virtual Item evaluateSingleton(const DynamicContext::Ptr &context) const;
    };

    /**
     * @short Implements the function <tt>fn:last()</tt>.
     *
     * @see <a href="http://www.w3.org/TR/xpath-functions/#func-last">XQuery 1.0
     * and XPath 2.0 Functions and Operators, 16.2 fn:last</a>
     * @ingroup Patternist_functions
     * @author Frans Englich <frans.englich@nokia.com>
     */
    class LastFN : public FunctionCall
    {
    public:
        virtual Item evaluateSingleton(const DynamicContext::Ptr &context) const;
    };

    /**
     * @short Implements the function <tt>fn:implicit-timezone()</tt>.
     *
     * @see <a href="http://www.w3.org/TR/xpath-functions/#func-implicit-timezone">XQuery 1.0
     * and XPath 2.0 Functions and Operators, 16.6 fn:implicit-timezone</a>
     * @ingroup Patternist_functions
     * @author Frans Englich <frans.englich@nokia.com>
     */
    class ImplicitTimezoneFN : public FunctionCall
    {
    public:
        virtual Item evaluateSingleton(const DynamicContext::Ptr &context) const;
    };

    /**
     * @short Implements the function <tt>fn:current-dateTime()</tt>.
     *
     * @see <a href="http://www.w3.org/TR/xpath-functions/#func-current-dateTime">XQuery 1.0
     * and XPath 2.0 Functions and Operators, 16.3 fn:current-dateTime</a>
     * @ingroup Patternist_functions
     * @author Frans Englich <frans.englich@nokia.com>
     */
    class CurrentDateTimeFN : public FunctionCall
    {
    public:
        virtual Item evaluateSingleton(const DynamicContext::Ptr &context) const;
    };

    /**
     * @short Implements the function <tt>fn:current-date()</tt>.
     *
     * @see <a href="http://www.w3.org/TR/xpath-functions/#func-current-date">XQuery 1.0
     * and XPath 2.0 Functions and Operators, 16.4 fn:current-date</a>
     * @ingroup Patternist_functions
     * @author Frans Englich <frans.englich@nokia.com>
     */
    class CurrentDateFN : public FunctionCall
    {
    public:
        virtual Item evaluateSingleton(const DynamicContext::Ptr &context) const;
    };

    /**
     * @short Implements the function <tt>fn:current-time()</tt>.
     *
     * @see <a href="http://www.w3.org/TR/xpath-functions/#func-current-time">XQuery 1.0
     * and XPath 2.0 Functions and Operators, 16.5 fn:current-date</a>
     * @ingroup Patternist_functions
     * @author Frans Englich <frans.englich@nokia.com>
     */
    class CurrentTimeFN : public FunctionCall
    {
    public:
        virtual Item evaluateSingleton(const DynamicContext::Ptr &context) const;
    };

    /**
     * @short Implements the function <tt>fn:default-collation()</tt>.
     *
     * This is done by rewriting to StaticContext::defaultCollation() at the typeCheck() stage.
     *
     * @see <a href="http://www.w3.org/TR/xpath-functions/#func-default-collation">XQuery 1.0
     * and XPath 2.0 Functions and Operators, 16.7 fn:default-collation</a>
     * @ingroup Patternist_functions
     * @author Frans Englich <frans.englich@nokia.com>
     */
    class DefaultCollationFN : public FunctionCall
    {
    public:
        virtual Expression::Ptr typeCheck(const StaticContext::Ptr &context,
                                          const SequenceType::Ptr &reqType);
    };

    /**
     * @short Implements the function <tt>fn:static-base-uri()</tt>.
     *
     * This is done by rewriting to StaticContext::baseURI() at the typeCheck() stage.
     *
     * @see <a href="http://www.w3.org/TR/xpath-functions/#func-static-base-uri">XQuery 1.0
     * and XPath 2.0 Functions and Operators, 16.8 fn:static-base-uri</a>
     * @ingroup Patternist_functions
     * @author Frans Englich <frans.englich@nokia.com>
     */
    class StaticBaseURIFN : public FunctionCall
    {
    public:
        virtual Expression::Ptr typeCheck(const StaticContext::Ptr &context,
                                          const SequenceType::Ptr &reqType);
    };
}

QT_END_NAMESPACE

#endif
