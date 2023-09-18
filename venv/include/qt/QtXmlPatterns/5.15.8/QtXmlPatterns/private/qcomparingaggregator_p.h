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

#ifndef Patternist_ComparingAggregator_H
#define Patternist_ComparingAggregator_H

/**
 * @file qcomparingaggregator_p.h
 * @short Contains the implementations for the functions <tt>fn:max()</tt>, MaxFN,
 * and <tt>fn:min()</tt>, MinFN, and the class ComparingAggregator.
 */

#include <private/qabstractfloat_p.h>
#include <private/qdecimal_p.h>
#include <private/qcastingplatform_p.h>
#include <private/qcomparisonplatform_p.h>
#include <private/qliteral_p.h>
#include <private/qaggregator_p.h>
#include <private/quntypedatomicconverter_p.h>
#include <private/qpatternistlocale_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Base class for the implementations of the <tt>fn:min()</tt> and <tt>fn:max()</tt> function.
     *
     * What function that more specifically is
     * followed, depends on how the constructor is called.
     *
     * @see MaxFN
     * @see MinFN
     * @ingroup Patternist_functions
     * @author Frans Englich <frans.englich@nokia.com>
     */
    template <AtomicComparator::Operator oper, AtomicComparator::ComparisonResult result>
    class ComparingAggregator : public Aggregator,
                                public ComparisonPlatform<ComparingAggregator<oper, result>,
                                                          true, AtomicComparator::AsValueComparison,
                                                          ReportContext::FORG0006>,
                                public CastingPlatform<ComparingAggregator<oper, result>, true>
    {
    public:
        virtual Item evaluateSingleton(const DynamicContext::Ptr &context) const;
        virtual Expression::Ptr typeCheck(const StaticContext::Ptr &context,
                                          const SequenceType::Ptr &reqType);

        inline AtomicComparator::Operator operatorID() const
        {
            return oper;
        }

        inline ItemType::Ptr targetType() const
        {
            return BuiltinTypes::xsDouble;
        }

    private:
        inline Item applyNumericPromotion(const Item &old,
                                               const Item &nev,
                                               const Item &newVal) const;

        using ComparisonPlatform<ComparingAggregator<oper, result>,
                                 true,
                                 AtomicComparator::AsValueComparison,
                                 ReportContext::FORG0006>::comparator;
        using ComparisonPlatform<ComparingAggregator<oper, result>,
                                 true,
                                 AtomicComparator::AsValueComparison,
                                 ReportContext::FORG0006>::fetchComparator;
        using CastingPlatform<ComparingAggregator<oper, result>, true>::cast;
    };

#include "qcomparingaggregator_tpl_p.h"

    /**
     * @short An instantiation of ComparingAggregator suitable for <tt>fn:max()</tt>.
     *
     * @ingroup Patternist_functions
     */
    typedef ComparingAggregator<AtomicComparator::OperatorGreaterThan, AtomicComparator::GreaterThan> MaxFN;

    /**
     * @short An instantiation of ComparingAggregator suitable for <tt>fn:max()</tt>.
     *
     * @ingroup Patternist_functions
     */
    typedef ComparingAggregator<AtomicComparator::OperatorLessThan, AtomicComparator::LessThan> MinFN;
}

QT_END_NAMESPACE

#endif
