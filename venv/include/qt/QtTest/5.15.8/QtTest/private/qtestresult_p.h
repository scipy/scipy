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

#ifndef QTESTRESULT_P_H
#define QTESTRESULT_P_H

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

QT_BEGIN_NAMESPACE

class QLatin1String;
class QStringView;

class QTestResultPrivate;
class QTestData;

class Q_TESTLIB_EXPORT QTestResult
{
public:
    static const char *currentTestObjectName();
    static bool currentTestFailed();
    static QTestData *currentTestData();
    static QTestData *currentGlobalTestData();
    static const char *currentTestFunction();
    static const char *currentDataTag();
    static const char *currentGlobalDataTag();
    static void finishedCurrentTestData();
    static void finishedCurrentTestDataCleanup();
    static void finishedCurrentTestFunction();
    static void reset();
    static void setBlacklistCurrentTest(bool b);

    static void addFailure(const char *message, const char *file, int line);
    static bool compare(bool success, const char *failureMsg,
                        char *val1, char *val2,
                        const char *actual, const char *expected,
                        const char *file, int line);
    static bool compare(bool success, const char *failureMsg,
                        double val1, double val2,
                        const char *actual, const char *expected,
                        const char *file, int line);
    static bool compare(bool success, const char *failureMsg,
                        float val1, float val2,
                        const char *actual, const char *expected,
                        const char *file, int line);
    static bool compare(bool success, const char *failureMsg,
                        int val1, int val2,
                        const char *actual, const char *expected,
                        const char *file, int line);
    static bool compare(bool success, const char *failureMsg,
                        unsigned val1, unsigned val2,
                        const char *actual, const char *expected,
                        const char *file, int line);
    static bool compare(bool success, const char *failureMsg,
                        QStringView val1, QStringView val2,
                        const char *actual, const char *expected,
                        const char *file, int line);
    static bool compare(bool success, const char *failureMsg,
                        const QLatin1String &val1, QStringView val2,
                        const char *actual, const char *expected,
                        const char *file, int line);
    static bool compare(bool success, const char *failureMsg,
                        QStringView val1, const QLatin1String &val2,
                        const char *actual, const char *expected,
                        const char *file, int line);
    static void setCurrentGlobalTestData(QTestData *data);
    static void setCurrentTestData(QTestData *data);
    static void setCurrentTestFunction(const char *func);
    static void setCurrentTestObject(const char *name);
    static void addSkip(const char *message, const char *file, int line);
    static bool expectFail(const char *dataIndex, const char *comment,
                           QTest::TestFailMode mode, const char *file, int line);
    static bool verify(bool statement, const char *statementStr, const char *extraInfo,
                       const char *file, int line);
    static void setSkipCurrentTest(bool value);
    static bool skipCurrentTest();

    static void setCurrentAppName(const char *appName);
    static const char *currentAppName();

private:
    Q_DISABLE_COPY(QTestResult)
};

QT_END_NAMESPACE

#endif
