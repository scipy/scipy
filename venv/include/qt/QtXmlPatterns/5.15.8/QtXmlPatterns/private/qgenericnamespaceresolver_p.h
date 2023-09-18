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

#ifndef Patternist_GenericNamespaceResolver_H
#define Patternist_GenericNamespaceResolver_H

#include <QHash>

#include <private/qnamespaceresolver_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Generic namespace resolver which resolves lookups against entries in a QHash.
     *
     * @ingroup Patternist
     * @author Frans Englich <frans.englich@nokia.com>
     */
    class GenericNamespaceResolver : public NamespaceResolver
    {
    public:
        GenericNamespaceResolver(const Bindings &list);
        virtual void addBinding(const QXmlName nb);

        virtual QXmlName::NamespaceCode lookupNamespaceURI(const QXmlName::PrefixCode prefix) const;

        /**
         * Returns a GenericNamespaceResolver containing the following bindings:
         *
         * - <tt>xml</tt> = <tt>http://www.w3.org/XML/1998/namespace</tt>
         * - <tt>xs</tt> = <tt>http://www.w3.org/2001/XMLSchema</tt>
         * - <tt>xsi</tt> = <tt>http://www.w3.org/2001/XMLSchema-instance</tt>
         * - <tt>fn</tt> = <tt>http://www.w3.org/2005/xpath-functions</tt>
         * - <tt>xdt</tt> = <tt>http://www.w3.org/2005/xpath-datatypes</tt>
         * - no prefix = empty namespace
         */
        static NamespaceResolver::Ptr defaultXQueryBindings();

        /**
         * Returns a GenericNamespaceResolver containing the following bindings:
         *
         * - <tt>xml</tt> = <tt>http://www.w3.org/XML/1998/namespace</tt>
         * - no prefix = empty namespace
         */
        static NamespaceResolver::Ptr defaultXSLTBindings();

        virtual Bindings bindings() const;

    private:
        /**
         * The key is the prefix, the value the namespace URI.
         */
        Bindings m_bindings;
    };
}

QT_END_NAMESPACE

#endif
