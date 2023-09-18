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

#ifndef Patternist_AtomicType_H
#define Patternist_AtomicType_H

#include <private/qanysimpletype_p.h>
#include <private/qatomiccasterlocator_p.h>
#include <private/qatomiccomparatorlocator_p.h>
#include <private/qatomicmathematicianlocator_p.h>
#include <private/qatomictypedispatch_p.h>
#include <private/qitemtype_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    class Item;
    class SourceLocationReflection;

    /**
     * @short Base class for all classes that implements atomic types.
     *
     * AtomicType does not implement @c xs:anyAtomicType, it is the C++
     * base class for classes that implement atomic types, such as @c xs:anyAtomicType.
     *
     * @ingroup Patternist_types
     * @author Frans Englich <frans.englich@nokia.com>
     */
    class AtomicType : public ItemType,
                       public AnySimpleType
    {
    public:

        typedef QExplicitlySharedDataPointer<AtomicType> Ptr;

        virtual ~AtomicType();

        /**
         * Implements a generic algorithm which relies on wxsTypeMatches().
         *
         * @returns @c true depending on if @p item is an atomic type, and that
         * AtomicValue::itemType()'s SequenceType::itemType() matches this type.
         */
        bool itemMatches(const Item &item) const override;

        /**
         * @returns the result of SharedQXmlName::displayName(), of the SharedQName
         * object returned from the name() function.
         */
        QString displayName(const NamePool::Ptr &np) const override;

        /**
         * returns always @c false
         */
        bool isNodeType() const override;

        /**
         * returns always @c true
         */
        bool isAtomicType() const override;

        /**
         * Determines whether @p other is equal to this type, or is a
         * sub-type of this type.
         *
         * The implementation is generic, relying on operator==()
         * and xdtSuperType().
         */
        bool xdtTypeMatches(const ItemType::Ptr &other) const override;

        /**
         * @returns always 'this'
         */
        ItemType::Ptr atomizedType() const override;

        /**
         * @returns always SchemaType::SimpleTypeAtomic
         */
        TypeCategory category() const override;

        /**
         * @returns DerivationRestriction
         */
        DerivationMethod derivationMethod() const override;

        virtual AtomicTypeVisitorResult::Ptr
        accept(const QExplicitlySharedDataPointer<AtomicTypeVisitor> &visitor,
               const SourceLocationReflection *const) const = 0;

        virtual AtomicTypeVisitorResult::Ptr
        accept(const ParameterizedAtomicTypeVisitor::Ptr &visitor,
               const qint16 param,
               const SourceLocationReflection *const) const = 0;

        virtual AtomicComparatorLocator::Ptr comparatorLocator() const = 0;
        virtual AtomicMathematicianLocator::Ptr mathematicianLocator() const = 0;
        virtual AtomicCasterLocator::Ptr casterLocator() const = 0;

    protected:
        AtomicType();

    };
}

QT_END_NAMESPACE

#endif
