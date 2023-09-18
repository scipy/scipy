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

#ifndef Patternist_CallTargetDescription_H
#define Patternist_CallTargetDescription_H

#include <QSharedData>

#include <private/qexpression_p.h>

QT_BEGIN_NAMESPACE

template<typename Key, typename Value> class QHash;
template<typename T> class QList;

namespace QPatternist
{
    class CallSite;

    /**
     * @short Contains metadata for a callable component, such as a function or
     * template.
     *
     * CallTargetDescription can be used directly and is so for templates, but
     * can also be sub-classed which FunctionSignature do.
     *
     * @ingroup Patternist_expr
     * @author Frans Englich <frans.englich@nokia.com>
     */
    class Q_AUTOTEST_EXPORT CallTargetDescription : public QSharedData
    {
    public:
        typedef QExplicitlySharedDataPointer<CallTargetDescription> Ptr;
        typedef QList<Ptr> List;

        CallTargetDescription(const QXmlName &name);

        /**
         * The function's name. For example, the name of the signature
         * <tt>fn:string() as xs:string</tt> is <tt>fn:string</tt>.
         */
        QXmlName name() const;

        /**
         * Flags callsites to be aware of their recursion by calling
         * UserFunctionCallsite::configureRecursion(), if that is the case.
         *
         * @note We pass @p expr by value here intentionally.
         */
        static void checkCallsiteCircularity(CallTargetDescription::List &signList,
                                             const Expression::Ptr expr);
    private:
        /**
         * Helper function for checkCallsiteCircularity(). If C++ allowed it,
         * it would have been local to it.
         */
        static void checkArgumentsCircularity(CallTargetDescription::List &signList,
                                              const Expression::Ptr callsite);

        Q_DISABLE_COPY(CallTargetDescription)
        const QXmlName m_name;
    };
}

QT_END_NAMESPACE

#endif

