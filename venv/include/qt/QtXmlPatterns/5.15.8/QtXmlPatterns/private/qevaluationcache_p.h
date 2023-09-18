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

#ifndef Patternist_EvaluationCache_H
#define Patternist_EvaluationCache_H

#include <private/qcachingiterator_p.h>
#include <private/qcommonsequencetypes_p.h>
#include <private/qnodebuilder_p.h>
#include <private/qoperandsiterator_p.h>
#include <private/qsinglecontainer_p.h>
#include <private/qvariabledeclaration_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Evaluates to the same result as its operand, but ensures the
     * operand is evaluated once even if this Expression is evaluated several
     * times.
     *
     * EvaluationCache does this in a pipelined way, by delivering items from
     * its cache, which is stored in the DynamicContext. If the cache has less
     * items than what the caller requests, EvaluationCache continues to
     * deliver but this time from the source, which it also populates into the
     * cache.
     *
     * EvaluationCache is used as an optimization in order to avoid running
     * expensive code paths multiple times, but also is sometimes a necessity:
     * for instance, when objects must be unique, such as potentially in the
     * case of node identity.
     *
     * EvaluationCache is in particular used for variables, whose sole purpose
     * is to store it once(at least conceptually) and then use it in multiple
     * places.
     *
     * In some cases an EvaluationCache isn't necessary. For instance, when a
     * variable is only referenced once. In those cases EvaluationCache removes
     * itself as an optimization; implemented in compress().
     *
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_expressions
     */
    template<bool IsForGlobal>
    class EvaluationCache : public SingleContainer
    {
    public:
        EvaluationCache(const Expression::Ptr &operand,
                        const VariableDeclaration::Ptr &varDecl,
                        const VariableSlotID slot);

        virtual Item evaluateSingleton(const DynamicContext::Ptr &context) const;
        virtual Item::Iterator::Ptr evaluateSequence(const DynamicContext::Ptr &context) const;
        virtual Expression::Ptr compress(const StaticContext::Ptr &context);

        virtual SequenceType::Ptr staticType() const;

        /**
         * The first operand must be exactly one @c xs:string.
         */
        virtual SequenceType::List expectedOperandTypes() const;

        virtual ExpressionVisitorResult::Ptr accept(const ExpressionVisitor::Ptr &visitor) const;
        virtual Properties properties() const;
        virtual Expression::Ptr typeCheck(const StaticContext::Ptr &context,
                                          const SequenceType::Ptr &reqType);
        virtual const SourceLocationReflection *actualReflection() const;

        inline VariableSlotID slot() const
        {
            return m_varSlot;
        }

    private:
        static DynamicContext::Ptr topFocusContext(const DynamicContext::Ptr &context);
        const VariableDeclaration::Ptr m_declaration;
        /**
         * This variable must not be called m_slot. If it so, a compiler bug on
         * HP-UX-aCC-64 is triggered in the constructor initializor. See the
         * preprocessor output.
         *
         * Note that this is the cache slot, and is disjoint to any variable's
         * regular slot.
         */
        const VariableSlotID            m_varSlot;
    };

#include "qevaluationcache_tpl_p.h"
}

QT_END_NAMESPACE

#endif
