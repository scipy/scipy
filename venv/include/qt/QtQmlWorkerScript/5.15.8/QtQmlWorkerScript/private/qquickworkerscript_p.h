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

#ifndef QQUICKWORKERSCRIPT_P_H
#define QQUICKWORKERSCRIPT_P_H

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

#include <qqml.h>

#include <QtQmlWorkerScript/private/qtqmlworkerscriptglobal_p.h>
#include <QtQml/qqmlparserstatus.h>
#include <QtCore/qthread.h>
#include <QtQml/qjsvalue.h>
#include <QtCore/qurl.h>

QT_BEGIN_NAMESPACE


class QQuickWorkerScript;
class QQuickWorkerScriptEnginePrivate;
class QQuickWorkerScriptEngine : public QThread
{
Q_OBJECT
public:
    QQuickWorkerScriptEngine(QQmlEngine *parent = nullptr);
    ~QQuickWorkerScriptEngine();

    int registerWorkerScript(QQuickWorkerScript *);
    void removeWorkerScript(int);
    void executeUrl(int, const QUrl &);
    void sendMessage(int, const QByteArray &);

protected:
    void run() override;

private:
    QQuickWorkerScriptEnginePrivate *d;
};

class QQmlV4Function;
class Q_AUTOTEST_EXPORT QQuickWorkerScript : public QObject, public QQmlParserStatus
{
    Q_OBJECT
    Q_PROPERTY(QUrl source READ source WRITE setSource NOTIFY sourceChanged)
    Q_PROPERTY(bool ready READ ready NOTIFY readyChanged REVISION 15)

    QML_NAMED_ELEMENT(WorkerScript);

    Q_INTERFACES(QQmlParserStatus)
public:
    QQuickWorkerScript(QObject *parent = nullptr);
    ~QQuickWorkerScript();

    QUrl source() const;
    void setSource(const QUrl &);

    bool ready() const;

public Q_SLOTS:
    void sendMessage(QQmlV4Function*);

Q_SIGNALS:
    void sourceChanged();
    Q_REVISION(15) void readyChanged();
    void message(const QJSValue &messageObject);

protected:
    void classBegin() override;
    void componentComplete() override;
    bool event(QEvent *) override;

private:
    QQuickWorkerScriptEngine *engine();
    QQuickWorkerScriptEngine *m_engine;
    int m_scriptId;
    QUrl m_source;
    bool m_componentComplete;
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QQuickWorkerScript)

#endif // QQUICKWORKERSCRIPT_P_H
