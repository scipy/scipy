/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Toolkit.
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

#ifndef QMEDIATIMERANGE_H
#define QMEDIATIMERANGE_H

#include <QtMultimedia/qtmultimediaglobal.h>
#include <QtMultimedia/qmultimedia.h>
#include <QtCore/qshareddata.h>

QT_BEGIN_NAMESPACE


class QMediaTimeRangePrivate;

class Q_MULTIMEDIA_EXPORT QMediaTimeInterval
{
public:
    QMediaTimeInterval();
    QMediaTimeInterval(qint64 start, qint64 end);
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    QMediaTimeInterval(const QMediaTimeInterval&);
    QMediaTimeInterval &operator=(const QMediaTimeInterval&) = default;
    QMediaTimeInterval(QMediaTimeInterval &&) = default;
    QMediaTimeInterval &operator=(QMediaTimeInterval &&) = default;
#endif

    qint64 start() const;
    qint64 end() const;

    bool contains(qint64 time) const;

    bool isNormal() const;
    QMediaTimeInterval normalized() const;
    QMediaTimeInterval translated(qint64 offset) const;

private:
    friend class QMediaTimeRangePrivate;
    friend class QMediaTimeRange;

    qint64 s;
    qint64 e;
};

Q_MULTIMEDIA_EXPORT bool operator==(const QMediaTimeInterval&, const QMediaTimeInterval&);
Q_MULTIMEDIA_EXPORT bool operator!=(const QMediaTimeInterval&, const QMediaTimeInterval&);

class Q_MULTIMEDIA_EXPORT QMediaTimeRange
{
public:

    QMediaTimeRange();
    QMediaTimeRange(qint64 start, qint64 end);
    QMediaTimeRange(const QMediaTimeInterval&);
    QMediaTimeRange(const QMediaTimeRange &range);
    ~QMediaTimeRange();

    QMediaTimeRange &operator=(const QMediaTimeRange&);
    QMediaTimeRange &operator=(const QMediaTimeInterval&);

    qint64 earliestTime() const;
    qint64 latestTime() const;

    QList<QMediaTimeInterval> intervals() const;
    bool isEmpty() const;
    bool isContinuous() const;

    bool contains(qint64 time) const;

    void addInterval(qint64 start, qint64 end);
    void addInterval(const QMediaTimeInterval &interval);
    void addTimeRange(const QMediaTimeRange&);

    void removeInterval(qint64 start, qint64 end);
    void removeInterval(const QMediaTimeInterval &interval);
    void removeTimeRange(const QMediaTimeRange&);

    QMediaTimeRange& operator+=(const QMediaTimeRange&);
    QMediaTimeRange& operator+=(const QMediaTimeInterval&);
    QMediaTimeRange& operator-=(const QMediaTimeRange&);
    QMediaTimeRange& operator-=(const QMediaTimeInterval&);

    void clear();

private:
    QSharedDataPointer<QMediaTimeRangePrivate> d;
};

Q_MULTIMEDIA_EXPORT bool operator==(const QMediaTimeRange&, const QMediaTimeRange&);
Q_MULTIMEDIA_EXPORT bool operator!=(const QMediaTimeRange&, const QMediaTimeRange&);
Q_MULTIMEDIA_EXPORT QMediaTimeRange operator+(const QMediaTimeRange&, const QMediaTimeRange&);
Q_MULTIMEDIA_EXPORT QMediaTimeRange operator-(const QMediaTimeRange&, const QMediaTimeRange&);

#ifndef QT_NO_DEBUG_STREAM
Q_MULTIMEDIA_EXPORT QDebug operator<<(QDebug, const QMediaTimeRange &);
#endif

QT_END_NAMESPACE


#endif  // QMEDIATIMERANGE_H
