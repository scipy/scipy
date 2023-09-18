/****************************************************************************
**
** Copyright (C) 2020 The Qt Company Ltd.
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

#ifndef QUICKTEST_P_H
#define QUICKTEST_P_H

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

#include <QtQuickTest/quicktest.h>

#include <QtQml/qqmlpropertymap.h>
#include <QtQml/qqml.h>

QT_BEGIN_NAMESPACE

class QTestRootObject : public QObject
{
    Q_OBJECT
    Q_PROPERTY(bool windowShown READ windowShown NOTIFY windowShownChanged)
    Q_PROPERTY(bool hasTestCase READ hasTestCase WRITE setHasTestCase NOTIFY hasTestCaseChanged)
    Q_PROPERTY(QObject *defined READ defined)
    QML_SINGLETON
    QML_ELEMENT

public:
    QTestRootObject(QObject *parent = nullptr)
        : QObject(parent), hasQuit(false), m_windowShown(false), m_hasTestCase(false)  {
        m_defined = new QQmlPropertyMap(this);
#if defined(QT_OPENGL_ES_2_ANGLE)
        m_defined->insert(QLatin1String("QT_OPENGL_ES_2_ANGLE"), QVariant(true));
#endif
    }

    static QTestRootObject *instance() {
        static QPointer<QTestRootObject> object = new QTestRootObject;
        if (!object) {
            // QTestRootObject was deleted when previous test ended, create a new one
            object = new QTestRootObject;
        }
        return object;
    }

    bool hasQuit:1;
    bool hasTestCase() const { return m_hasTestCase; }
    void setHasTestCase(bool value) { m_hasTestCase = value; emit hasTestCaseChanged(); }

    bool windowShown() const { return m_windowShown; }
    void setWindowShown(bool value) { m_windowShown = value; emit windowShownChanged(); }
    QQmlPropertyMap *defined() const { return m_defined; }

    void init() { setWindowShown(false); setHasTestCase(false); hasQuit = false; }

Q_SIGNALS:
    void windowShownChanged();
    void hasTestCaseChanged();

private Q_SLOTS:
    void quit() { hasQuit = true; }

private:
    bool m_windowShown : 1;
    bool m_hasTestCase :1;
    QQmlPropertyMap *m_defined;
};

QT_END_NAMESPACE

#endif // QUICKTEST_P_H
