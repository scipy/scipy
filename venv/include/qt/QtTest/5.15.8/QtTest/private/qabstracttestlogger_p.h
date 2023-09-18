/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtTest module of the Qt Toolkit.
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

#ifndef QABSTRACTTESTLOGGER_P_H
#define QABSTRACTTESTLOGGER_P_H

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

#include <QtTest/qttestglobal.h>

#include <stdio.h>
#include <stdlib.h>

QT_BEGIN_NAMESPACE

class QBenchmarkResult;
class QTestData;

class Q_TESTLIB_EXPORT QAbstractTestLogger
{
public:
    enum IncidentTypes {
        Pass,
        XFail,
        Fail,
        XPass,
        BlacklistedPass,
        BlacklistedFail,
        BlacklistedXPass,
        BlacklistedXFail
    };

    enum MessageTypes {
        Warn,
        QWarning,
        QDebug,
        QSystem,
        QFatal,
        Skip,
        Info,
        QInfo
    };

    QAbstractTestLogger(const char *filename);
    virtual ~QAbstractTestLogger();

    virtual void startLogging();
    virtual void stopLogging();

    virtual void enterTestFunction(const char *function) = 0;
    virtual void leaveTestFunction() = 0;

    virtual void enterTestData(QTestData *) {}

    virtual void addIncident(IncidentTypes type, const char *description,
                             const char *file = nullptr, int line = 0) = 0;
    virtual void addBenchmarkResult(const QBenchmarkResult &result) = 0;

    virtual void addMessage(QtMsgType, const QMessageLogContext &,
                            const QString &);

    virtual void addMessage(MessageTypes type, const QString &message,
                            const char *file = nullptr, int line = 0) = 0;

    bool isLoggingToStdout() const;

    void outputString(const char *msg);

protected:
    void filterUnprintable(char *str) const;
    FILE *stream;
};

struct QTestCharBuffer
{
    enum { InitialSize = 512 };

    inline QTestCharBuffer() : buf(staticBuf)
    {
        staticBuf[0] = '\0';
    }

    inline ~QTestCharBuffer()
    {
        if (buf != staticBuf)
            free(buf);
    }

    inline char *data()
    {
        return buf;
    }

    inline char **buffer()
    {
        return &buf;
    }

    inline const char* constData() const
    {
        return buf;
    }

    inline int size() const
    {
        return _size;
    }

    inline bool reset(int newSize)
    {
        char *newBuf = nullptr;
        if (buf == staticBuf) {
            // if we point to our internal buffer, we need to malloc first
            newBuf = reinterpret_cast<char *>(malloc(newSize));
        } else {
            // if we already malloc'ed, just realloc
            newBuf = reinterpret_cast<char *>(realloc(buf, newSize));
        }

        // if the allocation went wrong (newBuf == 0), we leave the object as is
        if (!newBuf)
            return false;

        _size = newSize;
        buf = newBuf;
        return true;
    }

private:
    int _size = InitialSize;
    char* buf;
    char staticBuf[InitialSize];
};

namespace QTest
{
    int qt_asprintf(QTestCharBuffer *buf, const char *format, ...);
}

namespace QTestPrivate
{
    enum IdentifierPart { TestObject = 0x1, TestFunction = 0x2, TestDataTag = 0x4, AllParts = 0xFFFF };
    void Q_TESTLIB_EXPORT generateTestIdentifier(QTestCharBuffer *identifier, int parts = AllParts);
}

QT_END_NAMESPACE

#endif
