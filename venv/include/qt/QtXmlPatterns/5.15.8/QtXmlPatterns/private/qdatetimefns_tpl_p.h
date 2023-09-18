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

/**
 * @file
 * @short This file is included by qdatetimefns_p.h.
 * If you need includes in this file, put them in qdatetimefns_p.h, outside of the namespace.
 */

template<typename TSubClass>
Item ExtractFromDurationFN<TSubClass>::evaluateSingleton(const DynamicContext::Ptr &context) const
{
    const Item item(m_operands.first()->evaluateSingleton(context));
    if(item)
    {
        return static_cast<const TSubClass *>(this)->
               extract(item.as<AbstractDuration>());
    }
    else
        return Item();
}

Item YearsFromDurationFN::extract(const AbstractDuration *const duration) const
{
    return Integer::fromValue(duration->years() * (duration->isPositive() ? 1 : -1));
}

Item MonthsFromDurationFN::extract(const AbstractDuration *const duration) const
{
    return Integer::fromValue(duration->months() * (duration->isPositive() ? 1 : -1));
}

Item DaysFromDurationFN::extract(const AbstractDuration *const duration) const
{
    return Integer::fromValue(duration->days() * (duration->isPositive() ? 1 : -1));
}

Item HoursFromDurationFN::extract(const AbstractDuration *const duration) const
{
    return Integer::fromValue(duration->hours() * (duration->isPositive() ? 1 : -1));
}

Item MinutesFromDurationFN::extract(const AbstractDuration *const duration) const
{
    return Integer::fromValue(duration->minutes() * (duration->isPositive() ? 1 : -1));
}

Item SecondsFromDurationFN::extract(const AbstractDuration *const duration) const
{
    return toItem(Decimal::fromValue((duration->seconds() + duration->mseconds() / 1000.0) *
                                     (duration->isPositive() ? 1 : -1)));
}

template<typename TSubClass>
Item ExtractFromDateTimeFN<TSubClass>::evaluateSingleton(const DynamicContext::Ptr &context) const
{
    const Item item(m_operands.first()->evaluateSingleton(context));
    if(item)
    {
        return static_cast<const TSubClass *>(this)->
               extract(item.as<AbstractDateTime>()->toDateTime());
    }
    else
        return Item();
}

Item YearFromAbstractDateTimeFN::extract(const QDateTime &dt) const
{
    return Integer::fromValue(dt.date().year());
}

Item DayFromAbstractDateTimeFN::extract(const QDateTime &dt) const
{
    return Integer::fromValue(dt.date().day());
}

Item MinutesFromAbstractDateTimeFN::extract(const QDateTime &dt) const
{
    return Integer::fromValue(dt.time().minute());
}

Item SecondsFromAbstractDateTimeFN::extract(const QDateTime &dt) const
{
    const QTime time(dt.time());
    return toItem(Decimal::fromValue(time.second() + time.msec() / 1000.0));
}

Item TimezoneFromAbstractDateTimeFN::extract(const QDateTime &dt) const
{
    if(dt.timeSpec() == Qt::UTC)
        return toItem(CommonValues::DayTimeDurationZero);
    else if(dt.timeSpec() == Qt::OffsetFromUTC)
        return toItem(DayTimeDuration::fromSeconds(dt.offsetFromUtc()));
    else
        return Item();
}

Item MonthFromAbstractDateTimeFN::extract(const QDateTime &dt) const
{
    return Integer::fromValue(dt.date().month());
}

Item HoursFromAbstractDateTimeFN::extract(const QDateTime &dt) const
{
    return Integer::fromValue(dt.time().hour());
}

