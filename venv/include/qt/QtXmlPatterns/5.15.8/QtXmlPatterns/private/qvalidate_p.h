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

#ifndef Patternist_Validate_H
#define Patternist_Validate_H

#include <private/qexpression_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Handles XQuery 1.0's <tt>validate</tt> expression.
     *
     * This class is currently not used. The Schema Validation Feature is not supported.
     *
     * @see <a href="http://www.w3.org/TR/xquery/#id-validate">XQuery 1.0: An XML
     * Query Language, 3.13 Validate Expressions</a>
     * @see <a href="http://www.w3.org/TR/xquery/#id-schema-validation-feature">XQuery 1.0: An
     * XML Query Language, 5.2.2 Schema Validation Feature</a>
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_expressions
     */
    class Validate
    {
    public:

        /**
         * Represents the validation mode.
         */
        enum Mode
        {
            Lax = 1,
            Strict
        };

        /**
         * Creates the necessary Expression instances
         * that validates the operand node @p operandNode in mode @p validationMode,
         * and returns it.
         */
        static Expression::Ptr create(const Expression::Ptr &operandNode,
                                      const Mode validationMode,
                                      const StaticContext::Ptr &context);
    private:
        Validate();
        Q_DISABLE_COPY(Validate)
    };
}

QT_END_NAMESPACE

#endif
