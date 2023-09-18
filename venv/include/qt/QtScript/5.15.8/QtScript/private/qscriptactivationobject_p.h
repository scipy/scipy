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

#ifndef QSCRIPTACTIVATIONOBJECT_P_H
#define QSCRIPTACTIVATIONOBJECT_P_H

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

#include "JSVariableObject.h"

QT_BEGIN_NAMESPACE

namespace QScript
{

class QScriptActivationObject : public JSC::JSVariableObject {
public:
    QScriptActivationObject(JSC::ExecState *callFrame, JSC::JSObject *delegate = 0);
    virtual ~QScriptActivationObject();
    virtual bool isDynamicScope() const { return true; }

    virtual bool getOwnPropertySlot(JSC::ExecState*, const JSC::Identifier& propertyName, JSC::PropertySlot&);
    virtual bool getOwnPropertyDescriptor(JSC::ExecState*, const JSC::Identifier& propertyName, JSC::PropertyDescriptor&);
    virtual void getOwnPropertyNames(JSC::ExecState*, JSC::PropertyNameArray&, JSC::EnumerationMode mode = JSC::ExcludeDontEnumProperties);

    virtual void putWithAttributes(JSC::ExecState *exec, const JSC::Identifier &propertyName, JSC::JSValue value, unsigned attributes);
    virtual void put(JSC::ExecState*, const JSC::Identifier& propertyName, JSC::JSValue value, JSC::PutPropertySlot&);
    virtual void put(JSC::ExecState*, unsigned propertyName, JSC::JSValue value);

    virtual bool deleteProperty(JSC::ExecState*, const JSC::Identifier& propertyName);

    virtual void defineGetter(JSC::ExecState*, const JSC::Identifier& propertyName, JSC::JSObject* getterFunction);
    virtual void defineSetter(JSC::ExecState*, const JSC::Identifier& propertyName, JSC::JSObject* setterFunction);
    virtual JSC::JSValue lookupGetter(JSC::ExecState*, const JSC::Identifier& propertyName);
    virtual JSC::JSValue lookupSetter(JSC::ExecState*, const JSC::Identifier& propertyName);

    virtual const JSC::ClassInfo* classInfo() const { return &info; }
    static const JSC::ClassInfo info;

    struct QScriptActivationObjectData : public JSVariableObjectData {
        QScriptActivationObjectData(JSC::Register* registers, JSC::JSObject *dlg)
            : JSVariableObjectData(&symbolTable, registers),
              delegate(dlg)
        { }
        JSC::SymbolTable symbolTable;
        JSC::JSObject *delegate;
    };

    JSC::JSObject *delegate() const
    { return d_ptr()->delegate; }
    void setDelegate(JSC::JSObject *delegate)
    { d_ptr()->delegate = delegate; }

    QScriptActivationObjectData *d_ptr() const { return static_cast<QScriptActivationObjectData *>(d); }
};

} // namespace QScript

QT_END_NAMESPACE

#endif
