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

#ifndef Patternist_EmptySequence_H
#define Patternist_EmptySequence_H

#include <private/qemptycontainer_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Implements the value instance of empty sequence: <tt>()</tt>.
     *
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_xdm
     */
    class EmptySequence : public EmptyContainer
    {
    public:
        /**
         * @short Creates an EmptySequence that is a replacement for @p
         * replacementFor.
         *
         * @see EmptySequence()
         */
        static Expression::Ptr create(const Expression *const replacementFor,
                                      const StaticContext::Ptr &context);


        /**
         * @short Creates an instance of EmptySequence.
         *
         * @note In most cases create() should be used, since it takes care of
         * adjusting source location annotations.
         *
         * @see create()
         */
        inline EmptySequence()
        {
        }

        virtual QString stringValue() const;

        /**
         * @returns always an empty iterator, an instance of EmptyIterator.
         */
        virtual Item::Iterator::Ptr evaluateSequence(const DynamicContext::Ptr &) const;

        /**
         * @returns always @c null.
         */
        virtual Item evaluateSingleton(const DynamicContext::Ptr &) const;

        /**
         * Does nothing.
         */
        virtual void evaluateToSequenceReceiver(const DynamicContext::Ptr &) const;

        /**
         * @returns always @c false.
         */
        virtual bool evaluateEBV(const DynamicContext::Ptr &context) const;

        /**
         * @returns always CommonSequenceTypes::Empty
         */
        virtual ItemType::Ptr type() const;

        /**
         * @returns always CommonSequenceTypes::Empty
         */
        virtual SequenceType::Ptr staticType() const;

        virtual ExpressionVisitorResult::Ptr accept(const ExpressionVisitor::Ptr &visitor) const;
        virtual ID id() const;
        virtual Properties properties() const;
    };
}

QT_END_NAMESPACE

#endif
