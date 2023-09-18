/****************************************************************************
**
** Copyright (C) 2015 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the QtScript module of the Qt Toolkit.
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

#ifndef QSCRIPTENGINEAGENT_P_H
#define QSCRIPTENGINEAGENT_P_H

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
#include "Debugger.h"
#include "qscriptengineagent.h"

#include "CallFrame.h"
#include "SourceCode.h"
#include "UString.h"
#include "DebuggerCallFrame.h"

QT_BEGIN_NAMESPACE

class QScriptEnginePrivate;

class QScriptEngineAgent;
class Q_SCRIPT_EXPORT QScriptEngineAgentPrivate : public JSC::Debugger
{
    Q_DECLARE_PUBLIC(QScriptEngineAgent)
public:
    static QScriptEngineAgent* get(QScriptEngineAgentPrivate* p) {return p->q_func();}
    static QScriptEngineAgentPrivate* get(QScriptEngineAgent* p) {return p->d_func();}

    QScriptEngineAgentPrivate(){}
    virtual ~QScriptEngineAgentPrivate(){}

    void attach();
    void detach();

    //scripts
    virtual void sourceParsed(JSC::ExecState*, const JSC::SourceCode&, int /*errorLine*/, const JSC::UString& /*errorMsg*/) {}
    virtual void scriptUnload(qint64 id)
    {
        q_ptr->scriptUnload(id);
    }
    virtual void scriptLoad(qint64 id, const JSC::UString &program,
                         const JSC::UString &fileName, int baseLineNumber)
    {
        q_ptr->scriptLoad(id,program, fileName, baseLineNumber);
    }

    //exceptions
    virtual void exception(const JSC::DebuggerCallFrame& frame, intptr_t sourceID, int lineno, bool hasHandler)
    {
        Q_UNUSED(frame);
        Q_UNUSED(sourceID);
        Q_UNUSED(lineno);
        Q_UNUSED(hasHandler);
    }
    virtual void exceptionThrow(const JSC::DebuggerCallFrame& frame, intptr_t sourceID, bool hasHandler);
    virtual void exceptionCatch(const JSC::DebuggerCallFrame& frame, intptr_t sourceID);

    //statements
    virtual void atStatement(const JSC::DebuggerCallFrame&, intptr_t sourceID, int lineno/*, int column*/);
    virtual void callEvent(const JSC::DebuggerCallFrame&, intptr_t sourceID, int lineno)
    {
        Q_UNUSED(lineno);
        q_ptr->contextPush();
        q_ptr->functionEntry(sourceID);
    }
    virtual void returnEvent(const JSC::DebuggerCallFrame& frame, intptr_t sourceID, int lineno);
    virtual void willExecuteProgram(const JSC::DebuggerCallFrame& frame, intptr_t sourceID, int lineno)
    {
        Q_UNUSED(frame);
        Q_UNUSED(sourceID);
        Q_UNUSED(lineno);
    }
    virtual void didExecuteProgram(const JSC::DebuggerCallFrame& frame, intptr_t sourceID, int lineno)
    {
        Q_UNUSED(frame);
        Q_UNUSED(sourceID);
        Q_UNUSED(lineno);
    }
    virtual void functionExit(const JSC::JSValue& returnValue, intptr_t sourceID);
    //others
    virtual void didReachBreakpoint(const JSC::DebuggerCallFrame& frame, intptr_t sourceID, int lineno/*, int column*/);

    virtual void evaluateStart(intptr_t sourceID)
    {
        q_ptr->functionEntry(sourceID);
    }
    virtual void evaluateStop(const JSC::JSValue& returnValue, intptr_t sourceID);

    QScriptEnginePrivate *engine;
    QScriptEngineAgent *q_ptr;
};

QT_END_NAMESPACE

#endif
