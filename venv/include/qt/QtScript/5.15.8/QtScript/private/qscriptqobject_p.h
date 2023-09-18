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

#ifndef QSCRIPTQOBJECT_P_H
#define QSCRIPTQOBJECT_P_H

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

#include "qscriptobject_p.h"

#include "qscriptengine.h"
#include <QtCore/qpointer.h>

#include "InternalFunction.h"

QT_BEGIN_NAMESPACE

namespace QScript
{

enum AttributeExtension {
    // ### Make sure there's no conflict with JSC::Attribute
    QObjectMemberAttribute = 1 << 12
};

class QObjectDelegate : public QScriptObjectDelegate
{
public:
    struct Data
    {
        QPointer<QObject> value;
        QScriptEngine::ValueOwnership ownership;
        QScriptEngine::QObjectWrapOptions options;

        QHash<QByteArray, JSC::JSValue> cachedMembers;

        Data(QObject *o, QScriptEngine::ValueOwnership own,
             QScriptEngine::QObjectWrapOptions opt)
            : value(o), ownership(own), options(opt) {}
    };

    QObjectDelegate(
        QObject *object, QScriptEngine::ValueOwnership ownership,
        const QScriptEngine::QObjectWrapOptions &options);
    ~QObjectDelegate();

    virtual Type type() const;

    virtual bool getOwnPropertySlot(QScriptObject*, JSC::ExecState*,
                                    const JSC::Identifier& propertyName,
                                    JSC::PropertySlot&);
    virtual bool getOwnPropertyDescriptor(QScriptObject*, JSC::ExecState*,
                                    const JSC::Identifier& propertyName,
                                    JSC::PropertyDescriptor&);

    virtual void put(QScriptObject*, JSC::ExecState* exec,
                     const JSC::Identifier& propertyName,
                     JSC::JSValue, JSC::PutPropertySlot&);
    virtual bool deleteProperty(QScriptObject*, JSC::ExecState*,
                                const JSC::Identifier& propertyName);
    virtual void getOwnPropertyNames(QScriptObject*, JSC::ExecState*,
                                     JSC::PropertyNameArray&,
                                     JSC::EnumerationMode mode = JSC::ExcludeDontEnumProperties);
    virtual void markChildren(QScriptObject*, JSC::MarkStack& markStack);
    virtual bool compareToObject(QScriptObject*, JSC::ExecState*, JSC::JSObject*);

    inline QObject *value() const { return data->value; }
    inline void setValue(QObject* value) { data->value = value; }
    inline bool hasParent() const { return data->value && data->value->parent(); }

    inline QScriptEngine::ValueOwnership ownership() const
        { return data->ownership; }
    inline void setOwnership(QScriptEngine::ValueOwnership ownership)
        { data->ownership = ownership; }

    inline QScriptEngine::QObjectWrapOptions options() const
        { return data->options; }
    inline void setOptions(QScriptEngine::QObjectWrapOptions options)
        { data->options = options; }

protected:
    Data *data;
};

class QObjectPrototypeObject : public QObject
{
    Q_OBJECT
public:
    QObjectPrototypeObject(QObject *parent = 0)
        : QObject(parent) { }
    ~QObjectPrototypeObject() { }
};

class QObjectPrototype : public QScriptObject
{
public:
    QObjectPrototype(JSC::ExecState*, WTF::PassRefPtr<JSC::Structure>,
                     JSC::Structure* prototypeFunctionStructure);
};

class QObjectConnectionManager;

struct QObjectWrapperInfo
{
    QObjectWrapperInfo(QScriptObject *obj,
                       QScriptEngine::ValueOwnership own,
                       const QScriptEngine::QObjectWrapOptions &opt)
        : object(obj), ownership(own), options(opt) {}

    QScriptObject *object;
    QScriptEngine::ValueOwnership ownership;
    QScriptEngine::QObjectWrapOptions options;

    // Returns true if this wrapper can be garbage-collected when there are no
    // other references to it in the JS environment (weak reference), otherwise
    // returns false (should not be collected).
    bool isCollectableWhenWeaklyReferenced() const
    {
        switch (ownership) {
        case QScriptEngine::ScriptOwnership:
            return true;
        case QScriptEngine::AutoOwnership: {
            QScriptObjectDelegate *delegate = object->delegate();
            Q_ASSERT(delegate && (delegate->type() == QScriptObjectDelegate::QtObject));
            return !static_cast<QObjectDelegate *>(delegate)->hasParent();
        }
        default:
            break;
        }
        return false;
    }
};

class QObjectData // : public QObjectUserData
{
public:
    QObjectData(QScriptEnginePrivate *engine);
    ~QObjectData();

    bool addSignalHandler(QObject *sender,
                          int signalIndex,
                          JSC::JSValue receiver,
                          JSC::JSValue slot,
                          JSC::JSValue senderWrapper,
                          Qt::ConnectionType type);
    bool removeSignalHandler(QObject *sender,
                             int signalIndex,
                             JSC::JSValue receiver,
                             JSC::JSValue slot);

    QScriptObject *findWrapper(QScriptEngine::ValueOwnership ownership,
                               const QScriptEngine::QObjectWrapOptions &options) const;
    void registerWrapper(QScriptObject *wrapper,
                         QScriptEngine::ValueOwnership ownership,
                         const QScriptEngine::QObjectWrapOptions &options);

    void clearConnectionMarkBits();
    int markConnections(JSC::MarkStack&);
    void markWrappers(JSC::MarkStack&);

private:
    QScriptEnginePrivate *engine;
    QScript::QObjectConnectionManager *connectionManager;
    QList<QScript::QObjectWrapperInfo> wrappers;
};

class QtFunction: public JSC::InternalFunction
{
public:
    // work around CELL_SIZE limitation
    struct Data
    {
        JSC::JSValue object;
        int initialIndex;
        bool maybeOverloaded;

        Data(JSC::JSValue o, int ii, bool mo)
            : object(o), initialIndex(ii), maybeOverloaded(mo) {}
    };

    QtFunction(JSC::JSValue object, int initialIndex, bool maybeOverloaded,
               JSC::JSGlobalData*, WTF::PassRefPtr<JSC::Structure>, const JSC::Identifier&);
    virtual ~QtFunction();

    virtual JSC::CallType getCallData(JSC::CallData&);
    virtual void markChildren(JSC::MarkStack&);

    virtual const JSC::ClassInfo* classInfo() const { return &info; }
    static const JSC::ClassInfo info;

    static JSC::JSValue JSC_HOST_CALL call(JSC::ExecState*, JSC::JSObject*,
                                           JSC::JSValue, const JSC::ArgList&);

    JSC::JSValue execute(JSC::ExecState *exec, JSC::JSValue thisValue,
                         const JSC::ArgList &args);

    QScriptObject *wrapperObject() const;
    QObject *qobject() const;
    const QMetaObject *metaObject() const;
    int initialIndex() const;
    int specificIndex(const QScriptContext *context) const;
    bool maybeOverloaded() const;
    int mostGeneralMethod(QMetaMethod *out = 0) const;
    QList<int> overloadedIndexes() const;

private:
    Data *data;
};

class QtPropertyFunction: public JSC::InternalFunction
{
public:
    // work around CELL_SIZE limitation
    struct Data
    {
        const QMetaObject *meta;
        int index;

        Data(const QMetaObject *m, int i)
            : meta(m), index(i) {}
    };

    QtPropertyFunction(const QMetaObject *meta, int index,
                       JSC::JSGlobalData*, WTF::PassRefPtr<JSC::Structure>,
                       const JSC::Identifier&);
    virtual ~QtPropertyFunction();

    virtual JSC::CallType getCallData(JSC::CallData&);

    virtual const JSC::ClassInfo* classInfo() const { return &info; }
    static const JSC::ClassInfo info;

    static JSC::JSValue JSC_HOST_CALL call(JSC::ExecState*, JSC::JSObject*,
                                           JSC::JSValue, const JSC::ArgList&);

    JSC::JSValue execute(JSC::ExecState *exec, JSC::JSValue thisValue,
                         const JSC::ArgList &args);

    const QMetaObject *metaObject() const;
    int propertyIndex() const;

private:
    Data *data;
};

class QMetaObjectWrapperObject : public JSC::JSObject
{
public:
    // work around CELL_SIZE limitation
    struct Data
    {
        const QMetaObject *value;
        JSC::JSValue ctor;
        JSC::JSValue prototype;

        Data(const QMetaObject *mo, JSC::JSValue c)
            : value(mo), ctor(c) {}
    };

    explicit QMetaObjectWrapperObject(
        JSC::ExecState *, const QMetaObject *metaobject, JSC::JSValue ctor,
        WTF::PassRefPtr<JSC::Structure> sid);
    ~QMetaObjectWrapperObject();

    virtual bool getOwnPropertySlot(JSC::ExecState*,
                                    const JSC::Identifier& propertyName,
                                    JSC::PropertySlot&);
    virtual bool getOwnPropertyDescriptor(JSC::ExecState*,
                                          const JSC::Identifier& propertyName,
                                          JSC::PropertyDescriptor&);
    virtual void put(JSC::ExecState* exec, const JSC::Identifier& propertyName,
                     JSC::JSValue, JSC::PutPropertySlot&);
    virtual bool deleteProperty(JSC::ExecState*,
                                const JSC::Identifier& propertyName);
    virtual void getOwnPropertyNames(JSC::ExecState*, JSC::PropertyNameArray&,
                                     JSC::EnumerationMode mode = JSC::ExcludeDontEnumProperties);
    virtual void markChildren(JSC::MarkStack& markStack);

    virtual JSC::CallType getCallData(JSC::CallData&);
    virtual JSC::ConstructType getConstructData(JSC::ConstructData&);

    virtual const JSC::ClassInfo* classInfo() const { return &info; }
    static const JSC::ClassInfo info;

    static JSC::JSValue JSC_HOST_CALL call(JSC::ExecState*, JSC::JSObject*,
                                           JSC::JSValue, const JSC::ArgList&);
    static JSC::JSObject* construct(JSC::ExecState *, JSC::JSObject *, const JSC::ArgList &);

    JSC::JSValue execute(JSC::ExecState *exec, const JSC::ArgList &args);

    inline const QMetaObject *value() const { return data->value; }
    inline void setValue(const QMetaObject* value) { data->value = value; }

    static WTF::PassRefPtr<JSC::Structure> createStructure(JSC::JSValue prototype)
    {
        return JSC::Structure::create(prototype, JSC::TypeInfo(JSC::ObjectType, StructureFlags));
    }

protected:
    static const unsigned StructureFlags = JSC::OverridesGetOwnPropertySlot | JSC::OverridesMarkChildren | JSC::OverridesGetPropertyNames | JSC::ImplementsHasInstance | JSObject::StructureFlags;

    Data *data;
};

class QMetaObjectPrototype : public QMetaObjectWrapperObject
{
public:
    QMetaObjectPrototype(JSC::ExecState*, WTF::PassRefPtr<JSC::Structure>,
                         JSC::Structure* prototypeFunctionStructure);
};

} // namespace QScript

QT_END_NAMESPACE

#endif
