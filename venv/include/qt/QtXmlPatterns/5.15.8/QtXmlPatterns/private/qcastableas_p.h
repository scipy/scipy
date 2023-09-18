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

#ifndef Patternist_CastableAs_H
#define Patternist_CastableAs_H

#include <private/qsinglecontainer_p.h>
#include <private/qcastingplatform_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{

    /**
     * @short Implements XPath 2.0's <tt>castable as</tt> expression.
     *
     * @see <a href="http://www.w3.org/TR/xpath20/#id-castable">XML Path Language
     * (XPath) 2.0, 3.10.3 Castable</a>
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_expressions
     */
    class CastableAs : public SingleContainer,
                       public CastingPlatform<CastableAs, false>
    {
    public:
        CastableAs(const Expression::Ptr &operand,
                             const SequenceType::Ptr &targetType);

        virtual bool evaluateEBV(const DynamicContext::Ptr &) const;

        /**
         * Overridden to const fold to @c true when the target type
         * is a type which casting to always succeeds. This is
         * the type identical to the target type, <tt>xs:string</tt>,
         * and <tt>xs:untypedAtomic</tt>.
         */
        virtual Expression::Ptr compress(const StaticContext::Ptr &context);
        virtual Expression::Ptr typeCheck(const StaticContext::Ptr &context,
                                          const SequenceType::Ptr &reqType);

        virtual SequenceType::List expectedOperandTypes() const;
        virtual SequenceType::Ptr staticType() const;

        virtual ExpressionVisitorResult::Ptr accept(const ExpressionVisitor::Ptr &visitor) const;

        inline ItemType::Ptr targetType() const
        {
            return m_targetType->itemType();
        }

    private:
        const SequenceType::Ptr m_targetType;
    };
}

QT_END_NAMESPACE

#endif
