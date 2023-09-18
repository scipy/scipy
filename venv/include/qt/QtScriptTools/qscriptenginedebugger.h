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

#ifndef QSCRIPTENGINEDEBUGGER_H
#define QSCRIPTENGINEDEBUGGER_H

#include <QtCore/qobject.h>
#include <QtScript/qtscriptglobal.h>

QT_BEGIN_NAMESPACE


class QAction;
class QScriptEngine;
class QWidget;
#ifndef QT_NO_MAINWINDOW
class QMainWindow;
#endif
class QMenu;
class QToolBar;

class QScriptEngineDebuggerPrivate;
class Q_SCRIPTTOOLS_EXPORT QScriptEngineDebugger : public QObject
{
    Q_OBJECT
public:
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

    enum DebuggerState {
        RunningState,
        SuspendedState
    };

    explicit QScriptEngineDebugger(QObject *parent = nullptr);
    ~QScriptEngineDebugger();

    void attachTo(QScriptEngine *engine);
    void detach();

    bool autoShowStandardWindow() const;
    void setAutoShowStandardWindow(bool autoShow);

#ifndef QT_NO_MAINWINDOW
    QMainWindow *standardWindow() const;
#endif
    QToolBar *createStandardToolBar(QWidget *parent = nullptr);
    QMenu *createStandardMenu(QWidget *parent = nullptr);

    QWidget *widget(DebuggerWidget widget) const;
    QAction *action(DebuggerAction action) const;

    DebuggerState state() const;

Q_SIGNALS:
    void evaluationSuspended();
    void evaluationResumed();

private:
    Q_DECLARE_PRIVATE(QScriptEngineDebugger)
    Q_DISABLE_COPY(QScriptEngineDebugger)

    Q_PRIVATE_SLOT(d_func(), void _q_showStandardWindow())
};

QT_END_NAMESPACE

#endif
