/****************************************************************************
**
** Copyright (C) 2018 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtSCriptTools module of the Qt Toolkit.
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

#ifndef QSCRIPTDEBUGGERCONSOLEGLOBALOBJECT_P_H
#define QSCRIPTDEBUGGERCONSOLEGLOBALOBJECT_P_H

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

#include <QtCore/qobject.h>

#include <QtCore/qstringlist.h>

#include "qscriptdebuggerconsolecommandgroupdata_p.h"
#include "qscriptdebuggerconsolecommand_p.h"

QT_BEGIN_NAMESPACE

class QScriptDebuggerCommandSchedulerInterface;
class QScriptMessageHandlerInterface;
class QScriptDebuggerResponseHandlerInterface;
class QScriptDebuggerConsole;
class QScriptDebuggerValue;
class QScriptDebuggerCommand;
class QScriptBreakpointData;

class QScriptDebuggerConsoleGlobalObjectPrivate;
class Q_AUTOTEST_EXPORT QScriptDebuggerConsoleGlobalObject
    : public QObject
{
    Q_OBJECT
public:
    QScriptDebuggerConsoleGlobalObject(QObject *parent = 0);
    ~QScriptDebuggerConsoleGlobalObject();

    QScriptDebuggerCommandSchedulerInterface *scheduler() const;
    void setScheduler(QScriptDebuggerCommandSchedulerInterface *scheduler);

    QScriptDebuggerResponseHandlerInterface *responseHandler() const;
    void setResponseHandler(QScriptDebuggerResponseHandlerInterface *responseHandler);

    QScriptMessageHandlerInterface *messageHandler() const;
    void setMessageHandler(QScriptMessageHandlerInterface *messageHandler);

    QScriptDebuggerConsole *console() const;
    void setConsole(QScriptDebuggerConsole *console);

public Q_SLOTS:
    // frontend
    int scheduleInterrupt();
    int scheduleContinue();
    int scheduleStepInto(int count = 1);
    int scheduleStepOver(int count = 1);
    int scheduleStepOut();
    int scheduleRunToLocation(const QString &fileName, int lineNumber);
    int scheduleRunToLocation(qint64 scriptId, int lineNumber);
    int scheduleForceReturn(int contextIndex, const QScriptDebuggerValue &value);

    int scheduleSetBreakpoint(const QScriptBreakpointData &data);
    int scheduleDeleteBreakpoint(int id);
    int scheduleDeleteAllBreakpoints();
    int scheduleGetBreakpoints();
    int scheduleGetBreakpointData(int id);
    int scheduleSetBreakpointData(int id, const QScriptBreakpointData &data);

    int scheduleGetScripts();
    int scheduleGetScriptData(qint64 id);
    int scheduleScriptsCheckpoint();
    int scheduleGetScriptsDelta();
    int scheduleResolveScript(const QString &fileName);

    int scheduleGetBacktrace();
    int scheduleGetThisObject(int contextIndex);
    int scheduleGetActivationObject(int contextIndex);
    int scheduleGetContextCount();
    int scheduleGetContextInfo(int contextIndex);

    int scheduleNewScriptValueIterator(const QScriptDebuggerValue &object);
    int scheduleGetPropertiesByIterator(int id, int count);
    int scheduleDeleteScriptValueIterator(int id);

    int scheduleEvaluate(int contextIndex, const QString &program,
                         const QString &fileName = QString(),
                         int lineNumber = 1);

    int scheduleScriptValueToString(const QScriptDebuggerValue &value);

    int scheduleClearExceptions();

    int scheduleCommand(const QScriptDebuggerCommand &command);

    // message handler
    void message(const QString &text, const QString &fileName = QString(),
                 int lineNumber = -1, int columnNumber = -1);
    void warning(const QString &text, const QString &fileName = QString(),
                 int lineNumber = -1, int columnNumber = -1);
    void error(const QString &text, const QString &fileName = QString(),
               int lineNumber = -1, int columnNumber = -1);

    // console state
    int getCurrentFrameIndex() const;
    void setCurrentFrameIndex(int index);
    qint64 getCurrentScriptId() const;
    void setCurrentScriptId(qint64 id);
    qint64 getSessionId() const;
    int getCurrentLineNumber() const;
    void setCurrentLineNumber(int lineNumber);

    // command introspection
    QScriptDebuggerConsoleCommandGroupMap getCommandGroups() const;
    QScriptDebuggerConsoleCommand *findCommand(const QString &command) const;
    QScriptDebuggerConsoleCommandList getCommandsInGroup(const QString &name) const;
    QStringList getCommandCompletions(const QString &prefix) const;

    bool checkSyntax(const QString &program);

    void setEvaluateAction(int action);

private:
    Q_DECLARE_PRIVATE(QScriptDebuggerConsoleGlobalObject)
    Q_DISABLE_COPY(QScriptDebuggerConsoleGlobalObject)
};

QT_END_NAMESPACE

#endif
