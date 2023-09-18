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

#ifndef QSCRIPTGLOBALOBJECT_P_H
#define QSCRIPTGLOBALOBJECT_P_H

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

#include "JSGlobalObject.h"

QT_BEGIN_NAMESPACE

namespace QScript
{

class GlobalObject : public JSC::JSGlobalObject
{
public:
    GlobalObject();
    virtual ~GlobalObject();
    virtual JSC::UString className() const { return "global"; }
    virtual void markChildren(JSC::MarkStack&);
    virtual bool getOwnPropertySlot(JSC::ExecState*,
                                    const JSC::Identifier& propertyName,
                                    JSC::PropertySlot&);
    virtual bool getOwnPropertyDescriptor(JSC::ExecState*,
                                          const JSC::Identifier& propertyName,
                                          JSC::PropertyDescriptor&);
    virtual void put(JSC::ExecState* exec, const JSC::Identifier& propertyName,
                     JSC::JSValue, JSC::PutPropertySlot&);
    virtual void putWithAttributes(JSC::ExecState* exec, const JSC::Identifier& propertyName,
                                   JSC::JSValue value, unsigned attributes);
    virtual bool deleteProperty(JSC::ExecState*,
                                const JSC::Identifier& propertyName);
    virtual void getOwnPropertyNames(JSC::ExecState*, JSC::PropertyNameArray&,
                                     JSC::EnumerationMode mode = JSC::ExcludeDontEnumProperties);
    virtual void defineGetter(JSC::ExecState*, const JSC::Identifier& propertyName, JSC::JSObject* getterFunction, unsigned attributes = 0);
    virtual void defineSetter(JSC::ExecState*, const JSC::Identifier& propertyName, JSC::JSObject* setterFunction, unsigned attributes = 0);
    virtual JSC::JSValue lookupGetter(JSC::ExecState*, const JSC::Identifier& propertyName);
    virtual JSC::JSValue lookupSetter(JSC::ExecState*, const JSC::Identifier& propertyName);

public:
    JSC::JSObject *customGlobalObject;
};

class OriginalGlobalObjectProxy : public JSC::JSObject
{
public:
    explicit OriginalGlobalObjectProxy(WTF::PassRefPtr<JSC::Structure> sid,
                                       JSC::JSGlobalObject *object)
        : JSC::JSObject(sid), originalGlobalObject(object)
    {}
    virtual ~OriginalGlobalObjectProxy()
    {}
    virtual JSC::UString className() const
    { return originalGlobalObject->className(); }
    virtual void markChildren(JSC::MarkStack& markStack)
    {
        markStack.append(originalGlobalObject);
        JSC::JSObject::markChildren(markStack);
    }
    virtual bool getOwnPropertySlot(JSC::ExecState* exec,
                                    const JSC::Identifier& propertyName,
                                    JSC::PropertySlot& slot)
    { return originalGlobalObject->JSC::JSGlobalObject::getOwnPropertySlot(exec, propertyName, slot); }
    virtual bool getOwnPropertyDescriptor(JSC::ExecState* exec,
                                          const JSC::Identifier& propertyName,
                                          JSC::PropertyDescriptor& descriptor)
    { return originalGlobalObject->JSC::JSGlobalObject::getOwnPropertyDescriptor(exec, propertyName, descriptor); }
    virtual void put(JSC::ExecState* exec, const JSC::Identifier& propertyName,
                     JSC::JSValue value, JSC::PutPropertySlot& slot)
    { originalGlobalObject->JSC::JSGlobalObject::put(exec, propertyName, value, slot); }
    virtual void putWithAttributes(JSC::ExecState* exec, const JSC::Identifier& propertyName, JSC::JSValue value, unsigned attributes)
    { originalGlobalObject->JSC::JSGlobalObject::putWithAttributes(exec, propertyName, value, attributes); }
    virtual bool deleteProperty(JSC::ExecState* exec,
                                const JSC::Identifier& propertyName)
    { return originalGlobalObject->JSC::JSGlobalObject::deleteProperty(exec, propertyName); }
    virtual void getOwnPropertyNames(JSC::ExecState* exec, JSC::PropertyNameArray& propertyNames, JSC::EnumerationMode mode = JSC::ExcludeDontEnumProperties)
    { originalGlobalObject->JSC::JSGlobalObject::getOwnPropertyNames(exec, propertyNames, mode); }
    virtual void defineGetter(JSC::ExecState* exec, const JSC::Identifier& propertyName, JSC::JSObject* getterFunction, unsigned attributes)
    { originalGlobalObject->JSC::JSGlobalObject::defineGetter(exec, propertyName, getterFunction, attributes); }
    virtual void defineSetter(JSC::ExecState* exec, const JSC::Identifier& propertyName, JSC::JSObject* setterFunction, unsigned attributes)
    { originalGlobalObject->JSC::JSGlobalObject::defineSetter(exec, propertyName, setterFunction, attributes); }
    virtual JSC::JSValue lookupGetter(JSC::ExecState* exec, const JSC::Identifier& propertyName)
    { return originalGlobalObject->JSC::JSGlobalObject::lookupGetter(exec, propertyName); }
    virtual JSC::JSValue lookupSetter(JSC::ExecState* exec, const JSC::Identifier& propertyName)
    { return originalGlobalObject->JSC::JSGlobalObject::lookupSetter(exec, propertyName); }
private:
    JSC::JSGlobalObject *originalGlobalObject;
};

} // namespace QScript

QT_END_NAMESPACE

#endif
