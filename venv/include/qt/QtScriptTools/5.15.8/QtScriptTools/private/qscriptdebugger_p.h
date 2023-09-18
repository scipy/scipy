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

#ifndef QSCRIPTDEBUGGER_P_H
#define QSCRIPTDEBUGGER_P_H

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

QT_BEGIN_NAMESPACE

class QScriptDebuggerFrontend;
class QScriptDebuggerConsoleWidgetInterface;
class QScriptDebuggerScriptsWidgetInterface;
class QScriptDebuggerCodeWidgetInterface;
class QScriptDebuggerCodeFinderWidgetInterface;
class QScriptBreakpointsWidgetInterface;
class QScriptDebuggerStackWidgetInterface;
class QScriptDebuggerLocalsWidgetInterface;
class QScriptDebugOutputWidgetInterface;
class QScriptErrorLogWidgetInterface;
class QScriptDebuggerWidgetFactoryInterface;
class QAction;
class QEvent;
class QMenu;
#ifndef QT_NO_TOOLBAR
class QToolBar;
#endif

class QScriptDebuggerPrivate;
class Q_AUTOTEST_EXPORT QScriptDebugger : public QObject
{
    Q_OBJECT
public:
    // mirrors QScriptEngineDebugger::DebuggerWidget
    enum DebuggerWidget {
        ConsoleWidget,
        StackWidget,
        ScriptsWidget,
        LocalsWidget,
        CodeWidget,
        CodeFinderWidget,
        BreakpointsWidget,
        DebugOutputWidget,
        ErrorLogWidget
    };
    // mirrors QScriptEngineDebugger::DebuggerAction
    enum DebuggerAction {
        InterruptAction,
        ContinueAction,
        StepIntoAction,
        StepOverAction,
        StepOutAction,
        RunToCursorAction,
        RunToNewScriptAction,
        ToggleBreakpointAction,
        ClearDebugOutputAction,
        ClearErrorLogAction,
        ClearConsoleAction,
        FindInScriptAction,
        FindNextInScriptAction,
        FindPreviousInScriptAction,
        GoToLineAction
    };

    QScriptDebugger(QObject *parent = 0);
    ~QScriptDebugger();

    QScriptDebuggerFrontend *frontend() const;
    void setFrontend(QScriptDebuggerFrontend *frontend);

    QWidget *widget(DebuggerWidget widget);
    QAction *action(DebuggerAction action, QObject *parent);

    QScriptDebuggerConsoleWidgetInterface *consoleWidget() const;
    void setConsoleWidget(QScriptDebuggerConsoleWidgetInterface *consoleWidget);

    QScriptDebuggerScriptsWidgetInterface *scriptsWidget() const;
    void setScriptsWidget(QScriptDebuggerScriptsWidgetInterface *scriptsWidget);

    QScriptDebuggerCodeWidgetInterface *codeWidget() const;
    void setCodeWidget(QScriptDebuggerCodeWidgetInterface *codeWidget);

    QScriptDebuggerCodeFinderWidgetInterface *codeFinderWidget() const;
    void setCodeFinderWidget(QScriptDebuggerCodeFinderWidgetInterface *codeFinderWidget);

    QScriptDebuggerStackWidgetInterface *stackWidget() const;
    void setStackWidget(QScriptDebuggerStackWidgetInterface *stackWidget);

    QScriptDebuggerLocalsWidgetInterface *localsWidget() const;
    void setLocalsWidget(QScriptDebuggerLocalsWidgetInterface *localsWidget);

    QScriptBreakpointsWidgetInterface *breakpointsWidget() const;
    void setBreakpointsWidget(QScriptBreakpointsWidgetInterface *breakpointsWidget);

    QScriptDebugOutputWidgetInterface *debugOutputWidget() const;
    void setDebugOutputWidget(QScriptDebugOutputWidgetInterface *debugOutputWidget);

    QScriptErrorLogWidgetInterface *errorLogWidget() const;
    void setErrorLogWidget(QScriptErrorLogWidgetInterface *errorLogWidget);

    QScriptDebuggerWidgetFactoryInterface *widgetFactory() const;
    void setWidgetFactory(QScriptDebuggerWidgetFactoryInterface *factory);

    QAction *interruptAction(QObject *parent) const;
    QAction *continueAction(QObject *parent) const;
    QAction *stepIntoAction(QObject *parent) const;
    QAction *stepOverAction(QObject *parent) const;
    QAction *stepOutAction(QObject *parent) const;
    QAction *runToCursorAction(QObject *parent) const;
    QAction *runToNewScriptAction(QObject *parent) const;

    QAction *toggleBreakpointAction(QObject *parent) const;

    QAction *findInScriptAction(QObject *parent) const;
    QAction *findNextInScriptAction(QObject *parent) const;
    QAction *findPreviousInScriptAction(QObject *parent) const;
    QAction *goToLineAction(QObject *parent) const;

    QAction *clearDebugOutputAction(QObject *parent) const;
    QAction *clearConsoleAction(QObject *parent) const;
    QAction *clearErrorLogAction(QObject *parent) const;

    QMenu *createStandardMenu(QWidget *widgetParent, QObject *actionParent);
#ifndef QT_NO_TOOLBAR
    QToolBar *createStandardToolBar(QWidget *widgetParent, QObject *actionParent);
#endif
    bool eventFilter(QObject *, QEvent *e);

    bool isInteractive() const;

Q_SIGNALS:
    void stopped() const;
    void started() const;

protected:
    void timerEvent(QTimerEvent *e);

protected:
    QScriptDebugger(QScriptDebuggerPrivate &dd, QObject *parent);

private:
    Q_DECLARE_PRIVATE(QScriptDebugger)
    Q_DISABLE_COPY(QScriptDebugger)

    Q_PRIVATE_SLOT(d_func(), void _q_onLineEntered(const QString &))
    Q_PRIVATE_SLOT(d_func(), void _q_onCurrentFrameChanged(int))
    Q_PRIVATE_SLOT(d_func(), void _q_onCurrentScriptChanged(qint64))
    Q_PRIVATE_SLOT(d_func(), void _q_onScriptLocationSelected(int))

    Q_PRIVATE_SLOT(d_func(), void _q_interrupt())
    Q_PRIVATE_SLOT(d_func(), void _q_continue())
    Q_PRIVATE_SLOT(d_func(), void _q_stepInto())
    Q_PRIVATE_SLOT(d_func(), void _q_stepOver())
    Q_PRIVATE_SLOT(d_func(), void _q_stepOut())
    Q_PRIVATE_SLOT(d_func(), void _q_runToCursor())
    Q_PRIVATE_SLOT(d_func(), void _q_runToNewScript())

    Q_PRIVATE_SLOT(d_func(), void _q_toggleBreakpoint())

    Q_PRIVATE_SLOT(d_func(), void _q_clearDebugOutput())
    Q_PRIVATE_SLOT(d_func(), void _q_clearErrorLog())
    Q_PRIVATE_SLOT(d_func(), void _q_clearConsole())

    Q_PRIVATE_SLOT(d_func(), void _q_findInScript())
    Q_PRIVATE_SLOT(d_func(), void _q_findNextInScript())
    Q_PRIVATE_SLOT(d_func(), void _q_findPreviousInScript())
    Q_PRIVATE_SLOT(d_func(), void _q_onFindCodeRequest(const QString &, int))
    Q_PRIVATE_SLOT(d_func(), void _q_goToLine())
};

QT_END_NAMESPACE

#endif
