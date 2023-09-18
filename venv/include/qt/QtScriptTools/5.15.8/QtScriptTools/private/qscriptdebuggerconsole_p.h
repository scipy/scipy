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

#ifndef QSCRIPTDEBUGGERCONSOLE_P_H
#define QSCRIPTDEBUGGERCONSOLE_P_H

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

#include <QtCore/qobjectdefs.h>
#include <QtCore/qscopedpointer.h>

#include "qscriptdebuggerconsolehistorianinterface_p.h"

QT_BEGIN_NAMESPACE

class QString;
class QScriptDebuggerConsoleCommandJob;
class QScriptMessageHandlerInterface;
class QScriptDebuggerCommandSchedulerInterface;
class QScriptDebuggerConsoleCommandManager;

class QScriptDebuggerConsolePrivate;
class Q_AUTOTEST_EXPORT QScriptDebuggerConsole
    : public QScriptDebuggerConsoleHistorianInterface
{
public:
    QScriptDebuggerConsole();
    ~QScriptDebuggerConsole();

    void loadScriptedCommands(const QString &scriptsPath,
                              QScriptMessageHandlerInterface *messageHandler);

    void showDebuggerInfoMessage(QScriptMessageHandlerInterface *messageHandler);

    QScriptDebuggerConsoleCommandManager *commandManager() const;

    QScriptDebuggerConsoleCommandJob *consumeInput(
        const QString &input,
        QScriptMessageHandlerInterface *messageHandler,
        QScriptDebuggerCommandSchedulerInterface *commandScheduler);
    bool hasIncompleteInput() const;
    QString incompleteInput() const;
    void setIncompleteInput(const QString &input);
    QString commandPrefix() const;

    int historyCount() const;
    QString historyAt(int index) const;
    void changeHistoryAt(int index, const QString &newHistory);

    int currentFrameIndex() const;
    void setCurrentFrameIndex(int index);

    qint64 currentScriptId() const;
    void setCurrentScriptId(qint64 id);

    int currentLineNumber() const;
    void setCurrentLineNumber(int lineNumber);

    int evaluateAction() const;
    void setEvaluateAction(int action);

    qint64 sessionId() const;
    void bumpSessionId();

private:
    QScopedPointer<QScriptDebuggerConsolePrivate> d_ptr;

    Q_DECLARE_PRIVATE(QScriptDebuggerConsole)
    Q_DISABLE_COPY(QScriptDebuggerConsole)
};

QT_END_NAMESPACE

#endif
