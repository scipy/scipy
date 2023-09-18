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

#ifndef Patternist_NamespaceBinding_H
#define Patternist_NamespaceBinding_H

template<typename T> class QVector;

#include <QXmlName>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Represents a namespace binding: a prefix, and a namespace URI.
     *
     * @author Frans Englich <frans.englich@nokia.com>
     */
    class NamespaceBinding
    {
    public:
        enum
        {
            InvalidCode = -1
        };

        typedef QVector<NamespaceBinding> Vector;

        inline NamespaceBinding() : m_prefix(InvalidCode),
                                    m_namespace(InvalidCode)
        {
        }

        inline NamespaceBinding(const QXmlName::PrefixCode p,
                                const QXmlName::NamespaceCode n) : m_prefix(p),
                                                                m_namespace(n)
        {
        }

        inline bool operator==(const NamespaceBinding &other) const
        {
            return m_prefix == other.m_prefix &&
                   m_namespace == other.m_namespace;
        }

        inline QXmlName::PrefixCode prefix() const
        {
            return m_prefix;
        }

        inline QXmlName::NamespaceCode namespaceURI() const
        {
            return m_namespace;
        }

        inline bool isNull() const
        {
            return m_prefix == InvalidCode;
        }

        /**
         * @short Constructs a NamespaceBinding whose prefix and namespace is
         * taken from @p qName.
         *
         * The local name in @p qName is ignored. @p qName may not be null.
         */
        static inline NamespaceBinding fromQXmlName(const QXmlName qName)
        {
            Q_ASSERT(!qName.isNull());
            return NamespaceBinding(qName.prefix(), qName.namespaceURI());
        }

    private:
        QXmlName::PrefixCode      m_prefix;
        QXmlName::NamespaceCode   m_namespace;
    };

    /**
     * @relates NamespaceBinding
     */
    static inline uint qHash(const NamespaceBinding nb)
    {
        return (nb.prefix() << 16) + nb.namespaceURI();
    }

}

QT_END_NAMESPACE

#endif
