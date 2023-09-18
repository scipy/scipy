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

#ifndef Patternist_VariableReference_H
#define Patternist_VariableReference_H

#include <private/qemptycontainer_p.h>

QT_BEGIN_NAMESPACE

template<typename T> class QList;

namespace QPatternist
{
    /**
     * @short Baseclass for classes being references to variables.
     *
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_expressions
     */
    class Q_AUTOTEST_EXPORT VariableReference : public EmptyContainer
    {
    public:
        typedef QExplicitlySharedDataPointer<VariableReference> Ptr;
        typedef QList<VariableReference::Ptr> List;

        /**
         * Creates a VariableReference.
         *
         * @param slot must be a valid slot. That is, zero or larger.
         */
        VariableReference(const VariableSlotID slot);
        ~VariableReference();

        /**
         * @returns the slot that this reference communicates through.
         *
         * This is a slot in the DynamicContext. Which one, depends on the
         * type, which this VariableReference does not have information about.
         * For instance, it could DynamicContext::expressionVariable() or
         * DynamicContext::rangeVariable().
         */
        inline VariableSlotID slot() const;

        /**
         * @returns DisableElimination
         */
        virtual Properties properties() const;

    private:
        /**
         * The slot. Same as returned by slot().
         *
         *  This variable is not called m_slot, because that creates a weird
         *  compiler error on hpuxi-acc. See the preprocessor output. EvaluationCache::m_varSlot
         *  is a similar workaround.
         */
        const VariableSlotID m_varSlot;
    };

    inline VariableSlotID VariableReference::slot() const
    {
        return m_varSlot;
    }

}

QT_END_NAMESPACE

#endif
