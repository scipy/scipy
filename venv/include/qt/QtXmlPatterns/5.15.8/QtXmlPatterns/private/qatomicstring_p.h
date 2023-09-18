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

#ifndef Patternist_String_H
#define Patternist_String_H

#include <QUrl>

#include <private/qitem_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Implements the value instance of the @c xs:string type.
     *
     *
     * This class was originally called String, and correspondingly the header
     * file was called String.h. However, this broke building on OS X, which
     * looks up file names case insensitively, and therefore found string.h.
     *
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_xdm
     * @todo Documentation is missing/incomplete
     */
    class Q_AUTOTEST_EXPORT AtomicString : public AtomicValue
    {
    public:
        friend class CommonValues;

        typedef AtomicValue::Ptr Ptr;

        /**
         * Creates an instance representing @p value.
         *
         * @note This function does not remove the string literal escaping allowed in XPath 2.0
         */
        static AtomicString::Ptr fromValue(const QString &value);

        static inline AtomicString::Ptr fromValue(const QUrl &value)
        {
            return fromValue(value.toString());
        }

        /**
         * Get the Effective %Boolean Value of this string. A zero-length
         * string has an effective boolean value of @c false, in all other cases @c true.
         *
         * @returns @c false if the contained string has a zero-length, otherwise @c true.
         */
        virtual bool evaluateEBV(const QExplicitlySharedDataPointer<DynamicContext> &) const;

        /**
         * The string value of a AtomicString instance is the value space.
         */
        virtual QString stringValue() const;

        virtual ItemType::Ptr type() const;

    protected:
        friend class StringComparator;
        friend class CompareFN;
        AtomicString(const QString &value);
        const QString m_value;
    };
}

QT_END_NAMESPACE

#endif
