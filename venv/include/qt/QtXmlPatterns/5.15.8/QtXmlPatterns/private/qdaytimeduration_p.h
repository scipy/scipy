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

#ifndef Patternist_DayTimeDuration_H
#define Patternist_DayTimeDuration_H

#include <private/qabstractduration_p.h>
#include <private/qitem_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Implements the value instance of the @c xs:dayTimeDuration type.
     *
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_xdm
     */
    class DayTimeDuration : public AbstractDuration
    {
    public:

        typedef QExplicitlySharedDataPointer<DayTimeDuration> Ptr;

        /**
         * Creates an instance from the lexical representation @p string.
         */
        static DayTimeDuration::Ptr fromLexical(const QString &string);

        static DayTimeDuration::Ptr fromComponents(const bool isPositive,
                                                   const DayCountProperty days,
                                                   const HourProperty hours,
                                                   const MinuteProperty minutes,
                                                   const SecondProperty seconds,
                                                   const MSecondProperty mseconds);
        /**
         * Creates a DayTimeDuration that has the value expressed in seconds @p secs
         * and milli seconds @p msecs. The signedness of @p secs communicates
         * whether this DayTimeDuration is positive or negative. @p msecs must always
         * be positive.
         */
        static DayTimeDuration::Ptr fromSeconds(const SecondCountProperty secs,
                                                const MSecondProperty msecs = 0);

        virtual ItemType::Ptr type() const;
        virtual QString stringValue() const;

        /**
         * @returns always 0.
         */
        virtual YearProperty years() const;

        /**
         * @returns always 0.
         */
        virtual MonthProperty months() const;
        virtual DayCountProperty days() const;
        virtual HourProperty hours() const;
        virtual MinuteProperty minutes() const;
        virtual MSecondProperty mseconds() const;
        virtual SecondProperty seconds() const;

        /**
         * @returns the value of this xs:dayTimeDuration
         * in milli seconds.
         * @see <a href="http://www.w3.org/TR/xpath-functions/#dt-dayTimeDuration">XQuery 1.0
         * and XPath 2.0 Functions and Operators, 10.3.2.2 Calculating the value of a
         * xs:dayTimeDuration from the lexical representation</a>
         */
        virtual Value value() const;

        /**
         * Creates a DayTimeDuration containing the value @p val. @p val is
         * expressed in milli seconds.
         *
         * If @p val is zero, is CommonValues::DayTimeDurationZero returned.
         */
        virtual Item fromValue(const Value val) const;

    protected:
        friend class CommonValues;

        DayTimeDuration(const bool isPositive,
                        const DayCountProperty days,
                        const HourProperty hours,
                        const MinuteProperty minutes,
                        const SecondProperty seconds,
                        const MSecondProperty mseconds);

    private:
        const DayCountProperty  m_days;
        const HourProperty      m_hours;
        const MinuteProperty    m_minutes;
        const SecondProperty    m_seconds;
        const MSecondProperty   m_mseconds;
    };
}

QT_END_NAMESPACE

#endif
