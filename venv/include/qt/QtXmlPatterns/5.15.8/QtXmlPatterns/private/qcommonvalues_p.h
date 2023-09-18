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

#ifndef Patternist_CommonValues_H
#define Patternist_CommonValues_H

#include <private/qdaytimeduration_p.h>
#include <private/qyearmonthduration_p.h>
#include <private/qemptyiterator_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{

    /**
     * @short A collection of common values.
     *
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_xdm
     * @todo Documentation is missing/incomplete
     */
    class CommonValues
    {
    public:
        /**
         * An empty, zero-length string.
         *
         * @note It is not @c null, but empty.
         */
        static const AtomicValue::Ptr EmptyString;

        /**
         * The string "true", the lexical representation of
         * @c xs:boolean's value @c true.
         */
        static const AtomicValue::Ptr TrueString;

        /**
         * The string "false", the lexical representation of
         * @c xs:boolean's value @c false.
         */
        static const AtomicValue::Ptr UntypedAtomicFalse;

        /**
         * The string "true", the lexical representation of
         * @c xs:boolean's value @c true.
         */
        static const AtomicValue::Ptr UntypedAtomicTrue;

        /**
         * The string "false", the lexical representation of
         * @c xs:boolean's value @c false.
         */
        static const AtomicValue::Ptr FalseString;

        /**
         * @returns a Boolean instance carrying the boolean value @c true.
         * Use this value instead of Boolean::fromValue() if you
         * know what boolean value you need.
         */
        static const AtomicValue::Ptr BooleanTrue;

        /**
         * @returns a Boolean instance carrying the boolean value @c true.
         * Use this value instead of Boolean::fromValue() if you
         * know what boolean value you need.
         */
        static const AtomicValue::Ptr BooleanFalse;

        /**
         * Not-a-Numeric typed as @c xs:double.
         */
        static const AtomicValue::Ptr DoubleNaN;

        /**
         * Not-a-Number typed as @c xs:float, <tt>xs:float("NaN")</tt>.
         */
        static const AtomicValue::Ptr FloatNaN;

        /**
         * Zero(0) typed as @c xs:integer, <tt>xs:integer("0")</tt>.
         */
        static const Item IntegerZero;

        /**
         * An empty, "", @c xs:anyURI.
         */
        static const AtomicValue::Ptr EmptyAnyURI;

        /**
         * The empty sequence.
         */
        static const EmptyIterator<Item>::Ptr emptyIterator;

        /**
         * <tt>xs:float("-INF")</tt>
         */
        static const AtomicValue::Ptr NegativeInfFloat;

        /**
         * <tt>xs:float("INF")</tt>
         */
        static const AtomicValue::Ptr InfFloat;

        /**
         * <tt>xs:double("-INF")</tt>
         */
        static const AtomicValue::Ptr NegativeInfDouble;

        /**
         * <tt>xs:double("INF")</tt>
         */
        static const AtomicValue::Ptr InfDouble;

        /**
         * <tt>xs:float("1")</tt>
         */
        static const AtomicValue::Ptr FloatOne;
        /**
         * <tt>xs:double("1")</tt>
         */
        static const AtomicValue::Ptr DoubleOne;
        /**
         * <tt>xs:decimal("1")</tt>
         */
        static const AtomicValue::Ptr DecimalOne;

        /**
         * <tt>xs:integer("1")</tt>
         */
        static const Item IntegerOne;

        /**
         * <tt>xs:integer("-1")</tt>
         */
        static const Item IntegerOneNegative;

        /**
         * <tt>xs:double("0")</tt>
         */
        static const AtomicValue::Ptr DoubleZero;

        /**
         * <tt>xs:float("0")</tt>
         */
        static const AtomicValue::Ptr FloatZero;
        /**
         * <tt>xs:integer("0")</tt>
         */
        static const AtomicValue::Ptr DecimalZero;

        /**
         * The @c xs:dayTimeDuration value PT0S
         */
        static const DayTimeDuration::Ptr DayTimeDurationZero;

        /**
         * The @c xs:yearMonthDuration value P0M
         */
        static const YearMonthDuration::Ptr YearMonthDurationZero;

    private:
        /**
         * The constructor is private because this class is not meant to be instantiated,
         * but should only be used via its static const members.
         */
        inline CommonValues();

        Q_DISABLE_COPY(CommonValues)
    };
}

QT_END_NAMESPACE

#endif
