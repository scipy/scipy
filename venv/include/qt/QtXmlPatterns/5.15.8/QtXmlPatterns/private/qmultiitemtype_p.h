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

#ifndef Patternist_MultiItemType_H
#define Patternist_MultiItemType_H

#include <QList>

#include <private/qitemtype_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Represents multiple types such as <tt>document()</tt> @em or <tt>xs:integer</tt>.
     *
     * In some situations two or more different types are allowed. For example, XQuery's
     * @c validate expression accepts document or element nodes(but not attribute
     * nodes, for example). MultiItemType is useful in such situations, its constructor
     * takes a list of ItemType instances which its member functions treats as a wholeness.
     *
     * For example, xdtTypeMatches() returns @c true if any of the represented types matches.
     *
     * @ingroup Patternist_types
     * @author Frans Englich <frans.englich@nokia.com>
     */
    class MultiItemType : public ItemType
    {
    public:
        /**
         * Creates a MultiItemType representing the types in @p typeList. @p typeList must
         * contain two or more types.
         */
        MultiItemType(const ItemType::List &typeList);

        /**
         * The display name are the names concatenated with "|" as separator. For example,
         * if this MultiItemType represents the types <tt>document()</tt>, <tt>xs:integer</tt>,
         * and <tt>xs:anyAtomicType</tt>, the display name is
         * "document() | xs:integer | xs:anyAtomicType".
         */
        QString displayName(const NamePool::Ptr &np) const override;

        /**
         * If any of the types this MultiItemType represents matches @p item, it is
         * considered a match.
         *
         * @returns @c true if any of the housed ItemType instances matches @p item, otherwise @c false
         */
        bool itemMatches(const Item &item) const override;

        /**
         * If any of the types this MultiItemType represents matches @p other, it is
         * considered a match.
         *
         * @returns @c true if any of the housed ItemType instances matches @p other, otherwise @c false
         */
        bool xdtTypeMatches(const ItemType::Ptr &other) const override;

        /**
         * @returns @c true if any of the represented types is a node type.
         */
        bool isNodeType() const override;

        /**
         * @returns @c true if any of the represented types is an atomic type.
         */
        bool isAtomicType() const override;

        /**
         * Determines the union type of all the represented types super types. For example,
         * if the represented types are <tt>xs:integer</tt>, <tt>document()</tt>
         * and <tt>xs:string</tt>, <tt>item()</tt> is returned.
         */
        ItemType::Ptr xdtSuperType() const override;

        /**
         * Determines the union type of all the represented types atomized types. For example,
         * if the represented types are <tt>xs:integer</tt> and <tt>document()</tt>,
         * <tt>xs:anyAtomicType</tt> is returned, because that's the super type of <tt>xs:integer</tt>
         * and <tt>xs:untypedAtomic</tt>.
         */
        ItemType::Ptr atomizedType() const override;

    private:
        const ItemType::List m_types;
        const ItemType::List::const_iterator m_end;
    };
}

QT_END_NAMESPACE

#endif
