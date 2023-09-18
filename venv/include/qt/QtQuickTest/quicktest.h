/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the test suite of the Qt Toolkit.
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

#ifndef QUICKTEST_H
#define QUICKTEST_H

#include <QtQuickTest/quicktestglobal.h>
#include <QtTest/qtest.h>

QT_BEGIN_NAMESPACE

class QQuickItem;

Q_QUICK_TEST_EXPORT int quick_test_main(int argc, char **argv, const char *name, const char *sourceDir);
Q_QUICK_TEST_EXPORT int quick_test_main_with_setup(int argc, char **argv, const char *name, const char *sourceDir, QObject *setup);

#ifdef QUICK_TEST_SOURCE_DIR

#define QUICK_TEST_MAIN(name) \
    int main(int argc, char **argv) \
    { \
        QTEST_SET_MAIN_SOURCE_PATH \
        return quick_test_main(argc, argv, #name, QUICK_TEST_SOURCE_DIR); \
    }

#define QUICK_TEST_OPENGL_MAIN(name) \
    int main(int argc, char **argv) \
    { \
        QTEST_SET_MAIN_SOURCE_PATH \
        return quick_test_main(argc, argv, #name, QUICK_TEST_SOURCE_DIR); \
    }

#define QUICK_TEST_MAIN_WITH_SETUP(name, QuickTestSetupClass) \
    int main(int argc, char **argv) \
    { \
        QTEST_SET_MAIN_SOURCE_PATH \
        QuickTestSetupClass setup; \
        return quick_test_main_with_setup(argc, argv, #name, QUICK_TEST_SOURCE_DIR, &setup); \
    }

#else

#define QUICK_TEST_MAIN(name) \
    int main(int argc, char **argv) \
    { \
        QTEST_SET_MAIN_SOURCE_PATH \
        return quick_test_main(argc, argv, #name, nullptr); \
    }

#define QUICK_TEST_OPENGL_MAIN(name) \
    int main(int argc, char **argv) \
    { \
        QTEST_SET_MAIN_SOURCE_PATH \
        return quick_test_main(argc, argv, #name, nullptr); \
    }

#define QUICK_TEST_MAIN_WITH_SETUP(name, QuickTestSetupClass) \
    int main(int argc, char **argv) \
    { \
        QTEST_SET_MAIN_SOURCE_PATH \
        QuickTestSetupClass setup; \
        return quick_test_main_with_setup(argc, argv, #name, nullptr, &setup); \
    }

#endif

namespace QQuickTest {
Q_QUICK_TEST_EXPORT bool qIsPolishScheduled(const QQuickItem *item);
Q_QUICK_TEST_EXPORT bool qWaitForItemPolished(const QQuickItem *item, int timeout = 5000);
}

QT_END_NAMESPACE

#endif
