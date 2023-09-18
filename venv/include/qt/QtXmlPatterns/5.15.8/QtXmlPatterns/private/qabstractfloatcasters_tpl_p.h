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
//

/*
 * NOTE: This file is included by qabstractfloatcasters_p.h
 * if you need some includes, put them in qabstractfloatcasters_p.h (outside of the namespace)
 */

template <const bool isDouble>
Item NumericToAbstractFloatCaster<isDouble>::castFrom(const Item &from,
                                                           const QExplicitlySharedDataPointer<DynamicContext> &) const
{
    // toDouble() returns same thing than toFloat()
    return toItem(AbstractFloat<isDouble>::fromValue(from.template as<Numeric>()->toDouble()));
}

template <const bool isDouble>
Item StringToAbstractFloatCaster<isDouble>::castFrom(const Item &from,
                                                          const QExplicitlySharedDataPointer<DynamicContext> &) const
{
    return toItem(AbstractFloat<isDouble>::fromLexical(from.stringValue()));
}

template <const bool isDouble>
Item BooleanToAbstractFloatCaster<isDouble>::castFrom(const Item &from,
                                                           const QExplicitlySharedDataPointer<DynamicContext> &context) const
{
    // RVCT doesn't like using template parameter in trinary operator when the trinary operator result is
    // passed directly into another constructor.
    bool tempDouble = isDouble;
    if(from.template as<AtomicValue>()->evaluateEBV(context))
        return tempDouble ? toItem(CommonValues::DoubleOne) : toItem(CommonValues::FloatOne);
    else
        return tempDouble ? toItem(CommonValues::DoubleZero) : toItem(CommonValues::FloatZero);
}

