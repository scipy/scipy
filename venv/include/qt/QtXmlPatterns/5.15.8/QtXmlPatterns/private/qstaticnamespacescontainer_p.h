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

#ifndef Patternist_StaticNamespacesContainer_H
#define Patternist_StaticNamespacesContainer_H

#include <private/qfunctioncall_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short A helper subclass that stores a NamespaceResolver for the static
     * namespaces.
     *
     * This is used by functionality which needs to resolve names against the
     * statically known namespaces, at runtime. A good example of this is @c
     * function-available().
     *
     * The resolver is accessed through staticNamespaces(), which will be
     * available after the typeCheck() stage.
     *
     * This class must be subclassed.
     *
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_functions
     */
    class StaticNamespacesContainer : public FunctionCall
    {
    public:
        /**
         * Reimplemented to store data from the @p context.
         */
        virtual Expression::Ptr typeCheck(const StaticContext::Ptr &context,
                                          const SequenceType::Ptr &reqType);

    protected:
        /**
         * Before typeCheck(), behavior of this function is undefined. After
         * typeCheck(), this function guarantees to return a valid pointer.
         */
        inline const NamespaceResolver::Ptr &staticNamespaces() const
        {
            Q_ASSERT(m_resolver);
            return m_resolver;
        }

        /**
         * This constructor only exists to ensure this class is subclassed.
         */
        inline StaticNamespacesContainer()
        {
        }

    private:
        NamespaceResolver::Ptr m_resolver;
    };
}

QT_END_NAMESPACE

#endif
