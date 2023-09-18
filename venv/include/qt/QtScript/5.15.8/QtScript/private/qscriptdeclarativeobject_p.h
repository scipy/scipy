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

#ifndef QSCRIPTDECLARATIVEOBJECT_P_H
#define QSCRIPTDECLARATIVEOBJECT_P_H

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

#include "config.h"
#include "qscriptobject_p.h"
#include "qscriptdeclarativeclass_p.h"

QT_BEGIN_NAMESPACE

class QScriptClass;

class QScriptDeclarativeClassPrivate
{
public:
    QScriptDeclarativeClassPrivate() : engine(0), q_ptr(0), context(0), supportsCall(false) {}

    QScriptEngine *engine;
    QScriptDeclarativeClass *q_ptr;
    QScriptContext *context;
    bool supportsCall:1;

    static QScriptDeclarativeClassPrivate *get(QScriptDeclarativeClass *c) {
        return c->d_ptr.data();
    }
};

namespace QScript
{

class DeclarativeObjectDelegate : public QScriptObjectDelegate
{
public:
    DeclarativeObjectDelegate(QScriptDeclarativeClass *c, QScriptDeclarativeClass::Object *o);
    ~DeclarativeObjectDelegate();

    virtual Type type() const;

    QScriptDeclarativeClass *scriptClass() const { return m_class; }
    QScriptDeclarativeClass::Object *object() const { return m_object; }

    virtual bool getOwnPropertySlot(QScriptObject*, JSC::ExecState*,
                                    const JSC::Identifier& propertyName,
                                    JSC::PropertySlot&);
    virtual void put(QScriptObject*, JSC::ExecState* exec,
                     const JSC::Identifier& propertyName,
                     JSC::JSValue, JSC::PutPropertySlot&);
    virtual bool deleteProperty(QScriptObject*, JSC::ExecState*,
                                const JSC::Identifier& propertyName);
    virtual void getOwnPropertyNames(QScriptObject*, JSC::ExecState*,
                                     JSC::PropertyNameArray&, 
                                     JSC::EnumerationMode mode = JSC::ExcludeDontEnumProperties);

    virtual JSC::CallType getCallData(QScriptObject*, JSC::CallData&);
    static JSC::JSValue JSC_HOST_CALL call(JSC::ExecState*, JSC::JSObject*,
                                           JSC::JSValue, const JSC::ArgList&);

    virtual JSC::ConstructType getConstructData(QScriptObject*, JSC::ConstructData&);

    virtual bool hasInstance(QScriptObject*, JSC::ExecState*,
                             JSC::JSValue value, JSC::JSValue proto);

    bool compareToObject(QScriptObject *, JSC::ExecState *, JSC::JSObject *);

private:
    QScriptDeclarativeClass *m_class;
    QScriptDeclarativeClass::Object *m_object;
};

} // namespace QScript

QT_END_NAMESPACE

#endif
