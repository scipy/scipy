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

#ifndef Patternist_CppCastingHelper_H
#define Patternist_CppCastingHelper_H

#include <QtCore/QtGlobal>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Provides convenience methods for performing static casts between C++ classes.
     *
     * In Patternist, it is very common to do up-casts from Expression or Item, which typically
     * involves writing messy code. Such an old-way cast looks like this:
     *
     * @code
     * static_cast<const MyClass *>(myInstance.data())->myClassMember()
     * @endcode
     *
     * CppCastingHelper provides the convenience method as() for this, which is functionally
     * equivalent to the above code, but simpler:
     *
     * @code
     * myInstance->as<MyClass>()->myClassMember()
     * @endcode
     *
     * The as() function performs a static cast.
     *
     * By using CppCastingHelper, this is achieved:
     *
     * - Const correctness is automatically taken care of
     * - Less code to write
     * - When compiling in debug mode, the as() functions uses a @c dynamic_cast to verify that the
     *   static casts are properly done, such that sensible error messages are given when the casts
     *   are invalid. It also traps invalid casts which nevertheless happen to work on a particular
     *   platform/compiler/hardware architecture.
     *
     * CppCastingHelper is a template class where the TSubClass parameter must be the class
     * inheriting CppCastingHelper. See Item or Expression for demonstration.
     *
     * @author Frans Englich <frans.englich@nokia.com>
     */
    template<typename TSubClass>
    class CppCastingHelper
    {
    public:

        /**
         * Casts this instance to:
         *
         * @code
         * const TCastTarget *
         * @endcode
         *
         * and returns the result.
         *
         * When compiled in debug mode, this function perform a @c dynamic_cast, in order to
         * check the correctness of the cast.
         */
        template<typename TCastTarget>
        inline const TCastTarget *as() const
        {
#if defined(Patternist_DEBUG) && !defined(Q_CC_XLC)
/* At least on aix-xlc-64, the compiler cries when it sees dynamic_cast. */
            Q_ASSERT_X(dynamic_cast<const TCastTarget *>(static_cast<const TSubClass *>(this)),
                       Q_FUNC_INFO,
                       "The cast is invalid. This class does not inherit the cast target.");
#endif
            return static_cast<const TCastTarget *>(static_cast<const TSubClass *>(this));
        }

        /**
         * Casts this instance to:
         *
         * @code
         * TCastTarget *
         * @endcode
         *
         * and returns the result.
         *
         * When compiled in debug mode, a @c dynamic_cast is attempted, in order to
         * check the correctness of the cast.
         */
        template<typename TCastTarget>
        inline TCastTarget *as()
        {
#if defined(Patternist_DEBUG) && !defined(Q_CC_XLC)
/* At least on aix-xlc-64, the compiler cries when it sees dynamic_cast. */
            Q_ASSERT_X(dynamic_cast<TCastTarget *>(static_cast<TSubClass *>(this)),
                       Q_FUNC_INFO,
                       "The cast is invalid. This class does not inherit the cast target.");
#endif
            return static_cast<TCastTarget *>(static_cast<TSubClass *>(this));
        }

    protected:
        /**
         * This constructor is protected because this class must be sub-classed.
         */
        inline CppCastingHelper() {}
    };
}

QT_END_NAMESPACE

#endif
