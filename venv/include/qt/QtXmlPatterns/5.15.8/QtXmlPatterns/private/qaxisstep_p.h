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

#ifndef Patternist_AxisStep_H
#define Patternist_AxisStep_H

#include <private/qemptycontainer_p.h>
#include <private/qitem_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short A step in a path expression that with an axis and a node test evaluates
     * to a sequence of nodes from the context item.
     *
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_expressions
     */
    class Q_AUTOTEST_EXPORT AxisStep : public EmptyContainer
    {
    public:
        AxisStep(const QXmlNodeModelIndex::Axis axis,
                 const ItemType::Ptr &nodeTest);

        virtual Item::Iterator::Ptr evaluateSequence(const DynamicContext::Ptr &) const;
        virtual Item evaluateSingleton(const DynamicContext::Ptr &) const;

        /**
         * Returns @p node if it matches the node test this step is using, otherwise @c null.
         */
        inline Item mapToItem(const QXmlNodeModelIndex &node,
                              const DynamicContext::Ptr &context) const;

        virtual SequenceType::List expectedOperandTypes() const;
        virtual SequenceType::Ptr staticType() const;

        /**
         * Rewrites to ParentNodeAxis, if possible.
         */
        virtual Expression::Ptr typeCheck(const StaticContext::Ptr &context,
                                          const SequenceType::Ptr &reqType);

        /**
         * @returns always BuiltinTypes::node;
         */
        virtual ItemType::Ptr expectedContextItemType() const;

        virtual ExpressionVisitorResult::Ptr accept(const ExpressionVisitor::Ptr &visitor) const;

        virtual Properties properties() const;

        /**
         * @returns the axis this step is using.
         */
        QXmlNodeModelIndex::Axis axis() const;

        /**
         * @returns the node test this step is using.
         */
        inline ItemType::Ptr nodeTest() const
        {
            return m_nodeTest;
        }

        void setNodeTest(const ItemType::Ptr &nev)
        {
            m_nodeTest = nev;
        }

        /**
         * @short Prints the EBNF name corresponding to @p axis.
         *
         * For instance, for QXmlNodeModelIndex::Child, "child" is returned.
         *
         * Apart from being used in this class, it is used in the SDK.
         */
        static QString axisName(const QXmlNodeModelIndex::Axis axis);

        virtual ID id() const;
        virtual PatternPriority patternPriority() const;

        inline void setAxis(const QXmlNodeModelIndex::Axis newAxis);

    private:
        typedef QExplicitlySharedDataPointer<const AxisStep> ConstPtr;

        static const QXmlNodeModelIndex::NodeKind s_whenAxisNodeKindEmpty[];

        /**
         * @returns @c true when the axis @p axis and a node test testing node of
         * type @p nodeKind always produces an empty sequence. One such example
         * is <tt>attribute::comment()</tt>.
         */
        static bool isAlwaysEmpty(const QXmlNodeModelIndex::Axis axis,
                                  const QXmlNodeModelIndex::NodeKind nodeKind);

        /**
         * The reason this variable is mutable, is that in the case of XSL-T patterns,
         * we do quite some reordering.
         */
        QXmlNodeModelIndex::Axis    m_axis;
        ItemType::Ptr               m_nodeTest;
    };

    void AxisStep::setAxis(const QXmlNodeModelIndex::Axis newAxis)
    {
        m_axis = newAxis;
    }

}

QT_END_NAMESPACE

#endif
