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

#ifndef Patternist_AbstractFloatMathematician_H
#define Patternist_AbstractFloatMathematician_H

#include <private/qabstractfloat_p.h>
#include <private/qatomicmathematician_p.h>
#include <private/qinteger_p.h>
#include <private/qschemanumeric_p.h>
#include <private/qpatternistlocale_p.h>
#include <private/qsourcelocationreflection_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Performs arithmetics between AbstractFloat values (Float and Double classes).
     *
     * @ingroup Patternist_xdm
     * @author Vincent Ricard <magic@magicninja.org>
     */
    template <const bool isDouble>
    class AbstractFloatMathematician : public AtomicMathematician
                                     , public DelegatingSourceLocationReflection
    {
    public:
        inline AbstractFloatMathematician(const SourceLocationReflection *const r) : DelegatingSourceLocationReflection(r)
        {
        }

        virtual Item calculate(const Item &o1,
                                    const Operator op,
                                    const Item &o2,
                                    const QExplicitlySharedDataPointer<DynamicContext> &context) const;
    };

#include "qabstractfloatmathematician_tpl_p.h"

    /**
     * An instantiation of AbstractFloatMathematician that handles @c xs:double.
     */
    typedef AbstractFloatMathematician<true> DoubleMathematician;

    /**
     * An instantiation of AbstractFloatMathematician that handles @c xs:float.
     */
    typedef AbstractFloatMathematician<false> FloatMathematician;
}

QT_END_NAMESPACE

#endif
