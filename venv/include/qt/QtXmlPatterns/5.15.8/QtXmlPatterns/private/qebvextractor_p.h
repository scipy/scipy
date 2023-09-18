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

#ifndef Patternist_EBVExtractor_H
#define Patternist_EBVExtractor_H

#include <private/qsinglecontainer_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Calculates the Effective Boolean Value of its operand.
     *
     * EBVExtractor performs functionality wise the same as @c fn:boolean(),
     * but does it without the dependencies which FunctionCall requires.
     *
     * There is code-duplication going on with BooleanFN.
     *
     * @see BooleanFN
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_expressions
     */
    class EBVExtractor : public SingleContainer
    {
    public:
        EBVExtractor(const Expression::Ptr &operand);
        virtual bool evaluateEBV(const DynamicContext::Ptr &context) const;
        virtual SequenceType::List expectedOperandTypes() const;
        virtual const SourceLocationReflection *actualReflection() const;
        virtual Expression::Ptr typeCheck(const StaticContext::Ptr &context,
                                          const SequenceType::Ptr &reqType);
        virtual ExpressionVisitorResult::Ptr accept(const ExpressionVisitor::Ptr &visitor) const;

        /**
         * @short Returns always @c xs:boolean.
         */
        virtual SequenceType::Ptr staticType() const;

        template<typename TSubClass, typename ThisType>
        static Expression::Ptr typeCheck(const StaticContext::Ptr &context,
                                         const SequenceType::Ptr &reqType,
                                         ThisType *const caller)
        {
            if(*CommonSequenceTypes::EBV->itemType() == *reqType->itemType())
                return caller->operands().first()->typeCheck(context, reqType);
            else
                return caller->TSubClass::typeCheck(context, reqType);
        }

    };
}

QT_END_NAMESPACE

#endif
