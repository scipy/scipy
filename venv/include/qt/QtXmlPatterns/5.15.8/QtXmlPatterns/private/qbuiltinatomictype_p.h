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

#ifndef Patternist_BuiltinAtomicType_H
#define Patternist_BuiltinAtomicType_H

#include <private/qatomictype_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{

    /**
     * @short Instances of this class represents types that are sub-classes
     * of @c xs:anyAtomicType.
     *
     * Retrieving instances of builtin types is done
     * via BuiltinTypesFactory::createSchemaType(), not by instantiating this
     * class manually.
     *
     * @ingroup Patternist_types
     * @author Frans Englich <frans.englich@nokia.com>
     */
    class BuiltinAtomicType : public AtomicType
    {
    public:

        typedef QExplicitlySharedDataPointer<BuiltinAtomicType> Ptr;

        /**
         * @returns always @c false
         */
        bool isAbstract() const override;

        /**
         * @returns the base type as specified in the constructors baseType argument.
         */
        SchemaType::Ptr wxsSuperType() const override;

        /**
         * @returns the same type as wxsSuperType(), except for the type @c xs:anyAtomicType, which
         * returns item()
         */
        ItemType::Ptr xdtSuperType() const override;

        AtomicComparatorLocator::Ptr comparatorLocator() const override;
        AtomicMathematicianLocator::Ptr mathematicianLocator() const override;
        AtomicCasterLocator::Ptr casterLocator() const override;

    protected:
        friend class BuiltinTypes;

        /**
         * @param baseType the type that is the super type of the constructed
         * atomic type. In the case of AnyAtomicType, @c null is passed.
         * @param comp the AtomicComparatorLocator this type should return. May be @c null.
         * @param mather similar to @p comp, this is the AtomicMathematicianLocator
         * that's appropriate for this type May be @c null.
         * @param casterLocator the CasterLocator that locates classes performing
         * casting with this type. May be @c null.
         */
        BuiltinAtomicType(const AtomicType::Ptr &baseType,
                          const AtomicComparatorLocator::Ptr &comp,
                          const AtomicMathematicianLocator::Ptr &mather,
                          const AtomicCasterLocator::Ptr &casterLocator);

    private:
        const AtomicType::Ptr                   m_superType;
        const AtomicComparatorLocator::Ptr      m_comparatorLocator;
        const AtomicMathematicianLocator::Ptr   m_mathematicianLocator;
        const AtomicCasterLocator::Ptr          m_casterLocator;
    };
}

QT_END_NAMESPACE

#endif
