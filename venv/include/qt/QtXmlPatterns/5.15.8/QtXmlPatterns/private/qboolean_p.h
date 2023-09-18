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

#ifndef Patternist_Boolean_H
#define Patternist_Boolean_H

#include <private/qitem_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Implements the value instance of the @c xs:boolean type.
     *
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_xdm
     */
    class Q_AUTOTEST_EXPORT Boolean : public AtomicValue
    {
    public:
        typedef AtomicValue::Ptr Ptr;

        /**
         * @returns the boolean value this Boolean represents
         */
        static bool evaluateEBV(const Item::Iterator::Ptr &e,
                                const QExplicitlySharedDataPointer<DynamicContext> &);

        static bool evaluateEBV(const Item &first,
                                const Item::Iterator::Ptr &e,
                                const QExplicitlySharedDataPointer<DynamicContext> &);

        static bool evaluateEBV(const Item &item,
                                const QExplicitlySharedDataPointer<DynamicContext> &context);

        virtual QString stringValue() const;

        /**
         * @returns a Boolean object instantiating @p value. Use True() or False()
         * if you already know what value you need.
         */
        static Boolean::Ptr fromValue(const bool value);

        /**
         * Creates a boolean value from a lexical representation. "true" and "1"
         * becomes @c true, while "false" and "0" becomes @c false.
         */
        static AtomicValue::Ptr fromLexical(const QString &val);

        /**
         * Get the Effective %Boolean Value of this boolean value. For <tt>xs:boolean</tt>, this
         * is simply the value.
         */
        virtual bool evaluateEBV(const QExplicitlySharedDataPointer<DynamicContext> &) const;

        virtual ItemType::Ptr type() const;

        inline bool value() const
        {
            return m_value;
        }

    protected:
        friend class CommonValues;
        Boolean(const bool value);

    private:
        const bool m_value;
    };
}

QT_END_NAMESPACE

#endif
