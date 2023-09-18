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

#ifndef Patternist_AnySimpleType_H
#define Patternist_AnySimpleType_H

#include <private/qanytype_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    class AtomicType;

    /**
     * @short Represents the @c xs:anySimpleType item type.
     *
     * @ingroup Patternist_types
     * @see <a href="http://www.w3.org/TR/xmlschema-2/#dt-anySimpleType">XML Schema Part 2:
     * Datatypes Second Edition, The simple ur-type definition</a>
     * @author Frans Englich <frans.englich@nokia.com>
     */
    class AnySimpleType : public AnyType
    {
    public:
        typedef QExplicitlySharedDataPointer<AnySimpleType> Ptr;
        typedef QList<AnySimpleType::Ptr> List;
        friend class BuiltinTypes;

        ~AnySimpleType();

        QXmlName name(const NamePool::Ptr &np) const override;

        /**
         * @returns always @c xs:anySimpleType
         */
        QString displayName(const NamePool::Ptr &np) const override;

        /**
         * @returns always BuiltinTypes::xsAnyType
         */
        SchemaType::Ptr wxsSuperType() const override;

        /**
         * xs:anySimpleType is the special "simple ur-type". Therefore this function
         * returns SchemaType::None
         *
         * @returns SchemaType::None
         */
        TypeCategory category() const override;

        /**
         * The simple ur-type is a "special restriction of the ur-type definition",
         * according to XML Schema Part 2: Datatypes Second Edition about xs:anySimpleType
         *
         * @returns DERIVATION_RESTRICTION
         */
        SchemaType::DerivationMethod derivationMethod() const override;

        /**
         * Always returns @c true.
         */
        bool isSimpleType() const override;

        /**
         * Always returns @c false.
         */
        bool isComplexType() const override;

    protected:
        AnySimpleType();

    };
}

QT_END_NAMESPACE

#endif
