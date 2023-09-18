/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtQml module of the Qt Toolkit.
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

#ifndef QQMLPROFILEREVENTLOCATION_P_H
#define QQMLPROFILEREVENTLOCATION_P_H

#include <QtCore/qstring.h>
#include <QtCore/qhash.h>
#include <QtCore/qdatastream.h>

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

QT_BEGIN_NAMESPACE

class QQmlProfilerEventLocation
{
public:
    QQmlProfilerEventLocation() : m_line(-1),m_column(-1) {}
    QQmlProfilerEventLocation(const QString &file, int lineNumber, int columnNumber) :
        m_filename(file), m_line(lineNumber), m_column(columnNumber)
    {}

    void clear()
    {
        m_filename.clear();
        m_line = m_column = -1;
    }

    bool isValid() const
    {
        return !m_filename.isEmpty();
    }

    QString filename() const { return m_filename; }
    int line() const { return m_line; }
    int column() const { return m_column; }

private:
    friend QDataStream &operator>>(QDataStream &stream, QQmlProfilerEventLocation &location);
    friend QDataStream &operator<<(QDataStream &stream, const QQmlProfilerEventLocation &location);

    QString m_filename;
    int m_line;
    int m_column;
};

inline bool operator==(const QQmlProfilerEventLocation &location1,
                       const QQmlProfilerEventLocation &location2)
{
    // compare filename last as it's expensive.
    return location1.line() == location2.line() && location1.column() == location2.column()
            && location1.filename() == location2.filename();
}

inline bool operator!=(const QQmlProfilerEventLocation &location1,
                       const QQmlProfilerEventLocation &location2)
{
    return !(location1 == location2);
}

inline uint qHash(const QQmlProfilerEventLocation &location)
{
    return qHash(location.filename())
            ^ ((location.line() & 0xfff)                   // 12 bits of line number
               | ((location.column() << 16) & 0xff0000));  // 8 bits of column

}

QDataStream &operator>>(QDataStream &stream, QQmlProfilerEventLocation &location);
QDataStream &operator<<(QDataStream &stream, const QQmlProfilerEventLocation &location);

Q_DECLARE_TYPEINFO(QQmlProfilerEventLocation, Q_MOVABLE_TYPE);

QT_END_NAMESPACE

#endif // QQMLPROFILEREVENTLOCATION_P_H
