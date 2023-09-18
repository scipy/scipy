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

#ifndef QSCRIPTENGINE_P_H
#define QSCRIPTENGINE_P_H

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

#include "private/qobject_p.h"

#include <QtCore/qdatetime.h>
#include <QtCore/qhash.h>
#include <QtCore/qnumeric.h>
#include <QtCore/qregexp.h>
#include <QtCore/qset.h>
#include <QtCore/qstringlist.h>
#include "qscriptvalue_p.h"
#include "qscriptstring_p.h"
#include "bridge/qscriptclassobject_p.h"
#include "bridge/qscriptdeclarativeclass_p.h"
#include "bridge/qscriptdeclarativeobject_p.h"
#include "bridge/qscriptobject_p.h"
#include "bridge/qscriptqobject_p.h"
#include "bridge/qscriptvariant_p.h"
#include "bridge/qscriptactivationobject_p.h"

#include "DateConstructor.h"
#include "DateInstance.h"
#include "Debugger.h"
#include "ErrorInstance.h"
#include "JSArray.h"
#include "Executable.h"
#include "Lexer.h"
#include "RefPtr.h"
#include "RegExpConstructor.h"
#include "RegExpObject.h"
#include "SourceProvider.h"
#include "Structure.h"
#include "UString.h"
#include "JSGlobalObject.h"
#include "JSValue.h"

#include <stdlib.h>

namespace JSC
{
    class EvalExecutable;
    class ExecState;
    typedef ExecState CallFrame;
    class JSCell;
    class JSGlobalObject;
}


QT_BEGIN_NAMESPACE

class QString;
class QStringList;
class QScriptContext;
class QScriptValue;
class QScriptTypeInfo;
class QScriptEngineAgent;
class QScriptEnginePrivate;
class QScriptSyntaxCheckResult;
class QScriptEngine;
class QScriptProgramPrivate;

namespace QScript
{
    class QObjectPrototype;
    class QMetaObjectPrototype;
    class QVariantPrototype;
#ifndef QT_NO_QOBJECT
    class QObjectData;
#endif
    class TimeoutCheckerProxy;

    qint32 ToInt32(qsreal);
    quint32 ToUInt32(qsreal);
    quint16 ToUInt16(qsreal);
    qsreal ToInteger(qsreal);

    inline bool ToBool(qsreal);
    inline bool ToBool(const QString &);
    inline qint32 ToInt32(const QString &);
    inline quint32 ToUInt32(const QString &);
    inline quint16 ToUInt16(const QString &);
    inline qsreal ToInteger(const QString &);
#ifdef Q_CC_MSVC
    // MSVC2008 crashes if these are inlined.
    qsreal ToNumber(const QString &);
    QString ToString(qsreal);
#else
    inline qsreal ToNumber(const QString &);
    inline QString ToString(qsreal);
#endif

    QDateTime MsToDateTime(JSC::ExecState *, qsreal);
    qsreal DateTimeToMs(JSC::ExecState *, const QDateTime &);

    //some conversion helper functions
    inline QScriptEnginePrivate *scriptEngineFromExec(const JSC::ExecState *exec);
    bool isFunction(JSC::JSValue value);

    inline void convertToLatin1_helper(const UChar *i, int length, char *s);
    inline QByteArray convertToLatin1(const JSC::UString &str);

    class UStringSourceProviderWithFeedback;

struct GlobalClientData : public JSC::JSGlobalData::ClientData
{
    GlobalClientData(QScriptEnginePrivate *e)
        : engine(e) {}
    virtual ~GlobalClientData() {}
    virtual void mark(JSC::MarkStack& markStack);
    virtual void uncaughtException(JSC::ExecState*, unsigned bytecodeOffset,
                                   JSC::JSValue);

    QScriptEnginePrivate *engine;
};

} // namespace QScript

class QScriptEnginePrivate
#ifndef QT_NO_QOBJECT
    : public QObjectPrivate
#endif
{
    Q_DECLARE_PUBLIC(QScriptEngine)
public:
    QScriptEnginePrivate();
    virtual ~QScriptEnginePrivate();

    static QScriptEnginePrivate *get(QScriptEngine *q) { return q ? q->d_func() : 0; }
    static QScriptEngine *get(QScriptEnginePrivate *d) { return d ? d->q_func() : 0; }

    static inline bool isArray(JSC::JSValue);
    static inline bool isDate(JSC::JSValue);
    static inline bool isError(JSC::JSValue);
    static inline bool isObject(JSC::JSValue);
    static inline bool isRegExp(JSC::JSValue);
    static inline bool isVariant(JSC::JSValue);
    static inline bool isQObject(JSC::JSValue);
    static inline bool isQMetaObject(JSC::JSValue);

    static inline bool toBool(JSC::ExecState *, JSC::JSValue);
    static inline qsreal toInteger(JSC::ExecState *, JSC::JSValue);
    static inline qsreal toNumber(JSC::ExecState *, JSC::JSValue);
    static inline qint32 toInt32(JSC::ExecState *, JSC::JSValue);
    static inline quint32 toUInt32(JSC::ExecState *, JSC::JSValue);
    static inline quint16 toUInt16(JSC::ExecState *, JSC::JSValue);
    static inline JSC::UString toString(JSC::ExecState *, JSC::JSValue);

    static inline QDateTime toDateTime(JSC::ExecState *, JSC::JSValue);
#ifndef QT_NO_REGEXP
    static QRegExp toRegExp(JSC::ExecState*, JSC::JSValue);
#endif
    static QVariant toVariant(JSC::ExecState *, JSC::JSValue);
    static inline QObject *toQObject(JSC::ExecState *, JSC::JSValue);
    static inline const QMetaObject *toQMetaObject(JSC::ExecState *, JSC::JSValue);

    static inline JSC::JSValue property(JSC::ExecState*, JSC::JSValue, const JSC::Identifier &id,
                                 int resolveMode = QScriptValue::ResolvePrototype);
    static JSC::JSValue propertyHelper(JSC::ExecState*, JSC::JSValue, const JSC::Identifier &id, int resolveMode);
    static inline JSC::JSValue property(JSC::ExecState*, JSC::JSValue, quint32 index,
                                 int resolveMode = QScriptValue::ResolvePrototype);
    static JSC::JSValue propertyHelper(JSC::ExecState*, JSC::JSValue, quint32, int resolveMode);
    static inline JSC::JSValue property(JSC::ExecState*, JSC::JSValue, const JSC::UString &, int resolveMode);
    static inline void setProperty(JSC::ExecState*, JSC::JSValue object, const JSC::UString &name, JSC::JSValue,
                     const QScriptValue::PropertyFlags &flags = QScriptValue::KeepExistingFlags);
    static void setProperty(JSC::ExecState*, JSC::JSValue object, const JSC::Identifier &id, JSC::JSValue,
                     const QScriptValue::PropertyFlags &flags = QScriptValue::KeepExistingFlags);
    static void setProperty(JSC::ExecState*, JSC::JSValue object, quint32 index, JSC::JSValue,
                     const QScriptValue::PropertyFlags &flags = QScriptValue::KeepExistingFlags);
    static QScriptValue::PropertyFlags propertyFlags(JSC::ExecState*, JSC::JSValue value,
                                              const JSC::Identifier &id, const QScriptValue::ResolveFlags &mode);
    static inline QScriptValue::PropertyFlags propertyFlags(JSC::ExecState*, JSC::JSValue value,
                                              const JSC::UString &name, const QScriptValue::ResolveFlags &mode);

    static bool convertValue(JSC::ExecState*, JSC::JSValue value,
                             int type, void *ptr);
    static bool convertNumber(qsreal, int type, void *ptr);
    static bool convertString(const QString &, int type, void *ptr);
    static JSC::JSValue create(JSC::ExecState*, int type, const void *ptr);
    bool hasDemarshalFunction(int type) const;

    inline QScriptValue scriptValueFromJSCValue(JSC::JSValue value);
    inline JSC::JSValue scriptValueToJSCValue(const QScriptValue &value);
    static inline unsigned propertyFlagsToJSCAttributes(const QScriptValue::PropertyFlags &flags);

    static inline JSC::JSValue jscValueFromVariant(JSC::ExecState*, const QVariant &value);
    static QVariant jscValueToVariant(JSC::ExecState*, JSC::JSValue value, int targetType);
    static inline QVariant &variantValue(JSC::JSValue value);
    static inline void setVariantValue(JSC::JSValue objectValue, const QVariant &value);

    static JSC::JSValue arrayFromStringList(JSC::ExecState*, const QStringList &lst);
    static QStringList stringListFromArray(JSC::ExecState*, JSC::JSValue arr);

    static JSC::JSValue arrayFromVariantList(JSC::ExecState*, const QVariantList &lst);
    static QVariantList variantListFromArray(JSC::ExecState*, JSC::JSArray *arr);

    static JSC::JSValue objectFromVariantMap(JSC::ExecState*, const QVariantMap &vmap);
    static QVariantMap variantMapFromObject(JSC::ExecState*, JSC::JSObject *obj);

    JSC::JSValue defaultPrototype(int metaTypeId) const;
    void setDefaultPrototype(int metaTypeId, JSC::JSValue prototype);

    static inline QScriptContext *contextForFrame(JSC::ExecState *frame);
    static inline JSC::ExecState *frameForContext(QScriptContext *context);
    static inline const JSC::ExecState *frameForContext(const QScriptContext *context);

    static inline bool hasValidCodeBlockRegister(JSC::ExecState *frame);

    JSC::JSGlobalObject *originalGlobalObject() const;
    JSC::JSObject *getOriginalGlobalObjectProxy();
    JSC::JSObject *customGlobalObject() const;
    JSC::JSObject *globalObject() const;
    void setGlobalObject(JSC::JSObject *object);
    inline JSC::ExecState *globalExec() const;
    JSC::JSValue toUsableValue(JSC::JSValue value);
    static JSC::JSValue thisForContext(JSC::ExecState *frame);
    static JSC::Register *thisRegisterForFrame(JSC::ExecState *frame);

    JSC::CallFrame *pushContext(JSC::CallFrame *exec, JSC::JSValue thisObject, const JSC::ArgList& args,
                                JSC::JSObject *callee, bool calledAsConstructor = false);
    void popContext();

    void mark(JSC::MarkStack& markStack);
    bool isCollecting() const;
    void collectGarbage();
    void reportAdditionalMemoryCost(int size);

    //flags that we set on the return value register for native function. (ie when codeBlock is 0)
    enum ContextFlags {
        NativeContext = 1,
        CalledAsConstructorContext = 2,
        HasScopeContext = 4, // Specifies that the is a QScriptActivationObject
        ShouldRestoreCallFrame = 8
    };
    static uint contextFlags(JSC::ExecState *);
    static void setContextFlags(JSC::ExecState *, uint);

    QScript::TimeoutCheckerProxy *timeoutChecker() const;

    void agentDeleted(QScriptEngineAgent *agent);

    static bool isLikelyStackOverflowError(JSC::ExecState *, JSC::JSValue);
    void uncaughtException(JSC::ExecState *, unsigned bytecodeOffset, JSC::JSValue);

    static inline void saveException(JSC::ExecState *, JSC::JSValue *);
    static inline void restoreException(JSC::ExecState *, JSC::JSValue);

    void setCurrentException(QScriptValue exception) { m_currentException = exception; }
    QScriptValue currentException() const { return m_currentException; }
    void clearCurrentException()
    {
        m_currentException.d_ptr.reset();
        uncaughtExceptionBacktrace.clear();
        uncaughtExceptionLineNumber = -1;
    }

    static QScriptSyntaxCheckResult checkSyntax(const QString &program);
    static bool canEvaluate(const QString &program);

    inline void registerScriptProgram(QScriptProgramPrivate *program);
    inline void unregisterScriptProgram(QScriptProgramPrivate *program);
    void detachAllRegisteredScriptPrograms();

    inline QScriptValuePrivate *allocateScriptValuePrivate(size_t);
    inline void freeScriptValuePrivate(QScriptValuePrivate *p);

    inline void registerScriptValue(QScriptValuePrivate *value);
    inline void unregisterScriptValue(QScriptValuePrivate *value);
    void detachAllRegisteredScriptValues();

    inline void registerScriptString(QScriptStringPrivate *value);
    inline void unregisterScriptString(QScriptStringPrivate *value);
    void detachAllRegisteredScriptStrings();
    QScriptString toStringHandle(const JSC::Identifier &name);

    static inline JSC::JSValue newArray(JSC::ExecState *, uint length);
    static inline JSC::JSValue newDate(JSC::ExecState *, qsreal value);
    static inline JSC::JSValue newDate(JSC::ExecState *, const QDateTime &);
    inline JSC::JSValue newObject();

#ifndef QT_NO_REGEXP
    static JSC::JSValue newRegExp(JSC::ExecState *, const QRegExp &);
#endif

    static JSC::JSValue newRegExp(JSC::ExecState *, const QString &pattern, const QString &flags);
    JSC::JSValue newVariant(const QVariant &);
    JSC::JSValue newVariant(JSC::JSValue objectValue, const QVariant &);

    static inline QScriptDeclarativeClass *declarativeClass(JSC::JSValue);
    static inline QScriptDeclarativeClass::Object *declarativeObject(JSC::JSValue);

    JSC::UString translationContextFromUrl(const JSC::UString &);

#ifndef QT_NO_QOBJECT
    void markQObjectData(JSC::MarkStack&);
    JSC::JSValue newQObject(QObject *object,
        QScriptEngine::ValueOwnership ownership = QScriptEngine::QtOwnership,
                 const QScriptEngine:: QObjectWrapOptions &options = {});
    JSC::JSValue newQMetaObject(const QMetaObject *metaObject,
                                JSC::JSValue ctor);

    static bool convertToNativeQObject(JSC::ExecState*, JSC::JSValue,
                                       const QByteArray &targetType,
                                       void **result);

    JSC::JSValue evaluateHelper(JSC::ExecState *exec, intptr_t sourceId,
                                JSC::EvalExecutable *executable,
                                bool &compile);

    QScript::QObjectData *qobjectData(QObject *object);
    void disposeQObject(QObject *object);
    void emitSignalHandlerException();

    bool scriptConnect(QObject *sender, const char *signal,
                       JSC::JSValue receiver, JSC::JSValue function,
                       Qt::ConnectionType type);
    bool scriptDisconnect(QObject *sender, const char *signal,
                          JSC::JSValue receiver, JSC::JSValue function);

    bool scriptConnect(QObject *sender, int index,
                       JSC::JSValue receiver, JSC::JSValue function,
                       JSC::JSValue senderWrapper,
                       Qt::ConnectionType type);
    bool scriptDisconnect(QObject *sender, int index,
                          JSC::JSValue receiver, JSC::JSValue function);

    bool scriptConnect(JSC::JSValue signal, JSC::JSValue receiver,
                       JSC::JSValue function, Qt::ConnectionType type);
    bool scriptDisconnect(JSC::JSValue signal, JSC::JSValue receiver,
                          JSC::JSValue function);

    // private slots
    void _q_objectDestroyed(QObject *);
#endif

    JSC::JSGlobalData *globalData;
    JSC::JSObject *originalGlobalObjectProxy;
    JSC::ExecState *currentFrame;

    WTF::RefPtr<JSC::Structure> scriptObjectStructure;
    WTF::RefPtr<JSC::Structure> staticScopeObjectStructure;

    QScript::QObjectPrototype *qobjectPrototype;
    WTF::RefPtr<JSC::Structure> qobjectWrapperObjectStructure;

    QScript::QMetaObjectPrototype *qmetaobjectPrototype;
    WTF::RefPtr<JSC::Structure> qmetaobjectWrapperObjectStructure;

    QScript::QVariantPrototype *variantPrototype;
    WTF::RefPtr<JSC::Structure> variantWrapperObjectStructure;

    QList<QScriptEngineAgent*> ownedAgents;
    QScriptEngineAgent *activeAgent;
    int agentLineNumber;
    QScriptValuePrivate *registeredScriptValues;
    QScriptValuePrivate *freeScriptValues;
    static const int maxFreeScriptValues = 256;
    int freeScriptValuesCount;
    QScriptStringPrivate *registeredScriptStrings;
    QSet<QScriptProgramPrivate*> registeredScriptPrograms;
    QHash<int, QScriptTypeInfo*> m_typeInfos;
    int processEventsInterval;
    QScriptValue abortResult;
    bool inEval;

    JSC::UString cachedTranslationUrl;
    JSC::UString cachedTranslationContext;

    QSet<QString> importedExtensions;
    QSet<QString> extensionsBeingImported;
    
    QHash<intptr_t, QScript::UStringSourceProviderWithFeedback*> loadedScripts;
    QScriptValue m_currentException;
    QStringList uncaughtExceptionBacktrace;
    int uncaughtExceptionLineNumber;

    QSet<JSC::JSObject*> visitedConversionObjects;

#ifndef QT_NO_QOBJECT
    QHash<QObject*, QScript::QObjectData*> m_qobjectData;
#endif

#ifdef QT_NO_QOBJECT
    QScriptEngine *q_ptr;
#endif
};

namespace QScript
{

class APIShim
{
public:
    APIShim(QScriptEnginePrivate *engine)
        : m_engine(engine), m_oldTable(JSC::setCurrentIdentifierTable(engine->globalData->identifierTable))
    {
    }
    ~APIShim()
    {
        JSC::setCurrentIdentifierTable(m_oldTable);
    }

private:
    QScriptEnginePrivate *m_engine;
    JSC::IdentifierTable *m_oldTable;
};

/*Helper class. Main purpose is to give debugger feedback about unloading and loading scripts.
  It keeps pointer to JSGlobalObject assuming that it is always the same - there is no way to update
  this data. Class is internal and used as an implementation detail in and only in QScriptEngine::evaluate.*/
class UStringSourceProviderWithFeedback: public JSC::UStringSourceProvider
{
public:
    static PassRefPtr<UStringSourceProviderWithFeedback> create(
        const JSC::UString& source, const JSC::UString& url,
        int lineNumber, QScriptEnginePrivate* engine)
    {
        return adoptRef(new UStringSourceProviderWithFeedback(source, url, lineNumber, engine));
    }

    /* Destruction means that there is no more copies of script so create scriptUnload event
       and unregister script in QScriptEnginePrivate::loadedScripts */
    virtual ~UStringSourceProviderWithFeedback()
    {
        if (m_ptr) {
            if (JSC::Debugger* debugger = this->debugger())
                debugger->scriptUnload(asID());
            m_ptr->loadedScripts.remove(asID());
        }
    }

    /* set internal QScriptEnginePrivate pointer to null and create unloadScript event, should be called
       only if QScriptEnginePrivate is about to be  destroyed.*/
    void disconnectFromEngine()
    {
        if (JSC::Debugger* debugger = this->debugger())
            debugger->scriptUnload(asID());
        m_ptr = 0;
    }

    int columnNumberFromOffset(int offset) const
    {
        for (const UChar *c = m_source.data() + offset; c >= m_source.data(); --c) {
            if (JSC::Lexer::isLineTerminator(*c))
                return offset - static_cast<int>(c - data());
        }
        return offset + 1;
    }

protected:
    UStringSourceProviderWithFeedback(const JSC::UString& source, const JSC::UString& url,
                                      int lineNumber, QScriptEnginePrivate* engine)
        : UStringSourceProvider(source, url),
          m_ptr(engine)
    {
        if (JSC::Debugger* debugger = this->debugger())
            debugger->scriptLoad(asID(), source, url, lineNumber);
        if (m_ptr)
            m_ptr->loadedScripts.insert(asID(), this);
    }

    JSC::Debugger* debugger()
    {
        //if m_ptr is null it mean that QScriptEnginePrivate was destroyed and scriptUnload was called
        //else m_ptr is stable and we can use it as normal pointer without hesitation
        if(!m_ptr)
            return 0; //we are in ~QScriptEnginePrivate
        else
            return m_ptr->originalGlobalObject()->debugger(); //QScriptEnginePrivate is still alive
    }

    //trace global object and debugger instance
    QScriptEnginePrivate* m_ptr;
};

class SaveFrameHelper
{
public:
    SaveFrameHelper(QScriptEnginePrivate *eng,
                    JSC::ExecState *newFrame)
        : engine(eng), oldFrame(eng->currentFrame)
    {
        eng->currentFrame = newFrame;
    }
    ~SaveFrameHelper()
    {
        engine->currentFrame = oldFrame;
    }
private:
    QScriptEnginePrivate *engine;
    JSC::ExecState *oldFrame;
};

inline QScriptEnginePrivate *scriptEngineFromExec(const JSC::ExecState *exec)
{
    return static_cast<GlobalClientData*>(exec->globalData().clientData)->engine;
}

#ifndef Q_CC_MSVC
// MSVC2008 crashes if these are inlined.

inline QString ToString(qsreal value)
{
    return JSC::UString::from(value);
}

inline qsreal ToNumber(const QString &value)
{
    return ((JSC::UString)value).toDouble();
}

#endif

inline qint32 ToInt32(const QString &value)
{
    return ToInt32(ToNumber(value));
}

inline quint32 ToUInt32(const QString &value)
{
    return ToUInt32(ToNumber(value));
}

inline quint16 ToUInt16(const QString &value)
{
    return ToUInt16(ToNumber(value));
}

inline qsreal ToInteger(const QString &value)
{
    return ToInteger(ToNumber(value));
}

inline bool ToBool(qsreal value)
{
    return (value != 0) && !qIsNaN(value);
}

inline bool ToBool(const QString &value)
{
     return !value.isEmpty();
}

inline void convertToLatin1_helper(const UChar *i, int length, char *s)
{
    const UChar *e = i + length;
    while (i != e)
        *(s++) = (uchar) *(i++);
    *s = '\0';
}

inline QByteArray convertToLatin1(const JSC::UString &str)
{
    QByteArray ba(str.size(), Qt::Uninitialized);
    convertToLatin1_helper(str.data(), str.size(), ba.data());
    return ba;
}

} // namespace QScript

inline void QScriptEnginePrivate::registerScriptProgram(QScriptProgramPrivate *program)
{
    Q_ASSERT(!registeredScriptPrograms.contains(program));
    registeredScriptPrograms.insert(program);
}

inline void QScriptEnginePrivate::unregisterScriptProgram(QScriptProgramPrivate *program)
{
    Q_ASSERT(registeredScriptPrograms.contains(program));
    registeredScriptPrograms.remove(program);
}

inline QScriptValuePrivate *QScriptEnginePrivate::allocateScriptValuePrivate(size_t size)
{
    if (freeScriptValues) {
        QScriptValuePrivate *p = freeScriptValues;
        freeScriptValues = p->next;
        --freeScriptValuesCount;
        return p;
    }
    return reinterpret_cast<QScriptValuePrivate*>(malloc(size));
}

inline void QScriptEnginePrivate::freeScriptValuePrivate(QScriptValuePrivate *p)
{
    if (freeScriptValuesCount < maxFreeScriptValues) {
        p->next = freeScriptValues;
        freeScriptValues = p;
        ++freeScriptValuesCount;
    } else {
        free(p);
    }
}

inline void QScriptEnginePrivate::registerScriptValue(QScriptValuePrivate *value)
{
    value->prev = 0;
    value->next = registeredScriptValues;
    if (registeredScriptValues)
        registeredScriptValues->prev = value;
    registeredScriptValues = value;
}

inline void QScriptEnginePrivate::unregisterScriptValue(QScriptValuePrivate *value)
{
    if (value->prev)
        value->prev->next = value->next;
    if (value->next)
        value->next->prev = value->prev;
    if (value == registeredScriptValues)
        registeredScriptValues = value->next;
    value->prev = 0;
    value->next = 0;
}

inline JSC::JSValue QScriptEnginePrivate::jscValueFromVariant(JSC::ExecState *exec, const QVariant &v)
{
    JSC::JSValue result = create(exec, v.userType(), v.data());
    Q_ASSERT(result);
    return result;
}

inline QScriptValue QScriptEnginePrivate::scriptValueFromJSCValue(JSC::JSValue value)
{
    if (!value)
        return QScriptValue();

    QScriptValuePrivate *p_value = new (this)QScriptValuePrivate(this);
    p_value->initFrom(value);
    return QScriptValuePrivate::toPublic(p_value);
}

inline JSC::JSValue QScriptEnginePrivate::scriptValueToJSCValue(const QScriptValue &value)
{
    QScriptValuePrivate *vv = QScriptValuePrivate::get(value);
    if (!vv)
        return JSC::JSValue();
    if (vv->type != QScriptValuePrivate::JavaScriptCore) {
        Q_ASSERT(!vv->engine || vv->engine == this);
        vv->engine = this;
        if (vv->type == QScriptValuePrivate::Number) {
            vv->initFrom(JSC::jsNumber(currentFrame, vv->numberValue));
        } else { //QScriptValuePrivate::String
            vv->initFrom(JSC::jsString(currentFrame, vv->stringValue));
        }
    }
    return vv->jscValue;
}

inline unsigned QScriptEnginePrivate::propertyFlagsToJSCAttributes(const QScriptValue::PropertyFlags &flags)
{
    unsigned attribs = 0;
    if (flags & QScriptValue::ReadOnly)
        attribs |= JSC::ReadOnly;
    if (flags & QScriptValue::SkipInEnumeration)
        attribs |= JSC::DontEnum;
    if (flags & QScriptValue::Undeletable)
        attribs |= JSC::DontDelete;
    attribs |= flags & QScriptValue::UserRange;
    return attribs;
}

inline QScriptValuePrivate::~QScriptValuePrivate()
{
    if (engine)
        engine->unregisterScriptValue(this);
}

inline void QScriptValuePrivate::initFrom(JSC::JSValue value)
{
    if (value.isCell()) {
        Q_ASSERT(engine != 0);
        value = engine->toUsableValue(value);
    }
    type = JavaScriptCore;
    jscValue = value;
    if (engine)
        engine->registerScriptValue(this);
}

inline void QScriptValuePrivate::initFrom(qsreal value)
{
    type = Number;
    numberValue = value;
    if (engine)
        engine->registerScriptValue(this);
}

inline void QScriptValuePrivate::initFrom(const QString &value)
{
    type = String;
    stringValue = value;
    if (engine)
        engine->registerScriptValue(this);
}

inline JSC::JSValue QScriptEnginePrivate::property(JSC::ExecState *exec, JSC::JSValue value, const JSC::UString &name, int resolveMode)
{
    return property(exec, value, JSC::Identifier(exec, name), resolveMode);
}

inline JSC::JSValue QScriptEnginePrivate::property(JSC::ExecState *exec, JSC::JSValue value, const JSC::Identifier &id, int resolveMode)
{
    Q_ASSERT(isObject(value));
    JSC::JSObject *object = JSC::asObject(value);
    JSC::PropertySlot slot(object);
    if ((resolveMode & QScriptValue::ResolvePrototype) && object->getPropertySlot(exec, id, slot))
        return slot.getValue(exec, id);
    return propertyHelper(exec, value, id, resolveMode);
}

inline JSC::JSValue QScriptEnginePrivate::property(JSC::ExecState *exec, JSC::JSValue value, quint32 index, int resolveMode)
{
    Q_ASSERT(isObject(value));
    JSC::JSObject *object = JSC::asObject(value);
    JSC::PropertySlot slot(object);
    if ((resolveMode & QScriptValue::ResolvePrototype) && object->getPropertySlot(exec, index, slot))
        return slot.getValue(exec, index);
    return propertyHelper(exec, value, index, resolveMode);
}

inline QScriptValue::PropertyFlags QScriptEnginePrivate::propertyFlags(JSC::ExecState *exec, JSC::JSValue value,
                                                                       const JSC::UString &name,
                                                                       const QScriptValue::ResolveFlags &mode)
{
    return propertyFlags(exec, value, JSC::Identifier(exec, name), mode);
}

inline void QScriptEnginePrivate::setProperty(JSC::ExecState *exec, JSC::JSValue objectValue, const JSC::UString &name,
                                              JSC::JSValue value, const QScriptValue::PropertyFlags &flags)
{
    setProperty(exec, objectValue, JSC::Identifier(exec, name), value, flags);
}

inline JSC::JSValue QScriptValuePrivate::property(const JSC::Identifier &id, const QScriptValue::ResolveFlags &resolveMode) const
{
    return QScriptEnginePrivate::property(engine->currentFrame, jscValue, id, resolveMode);
}

inline JSC::JSValue QScriptValuePrivate::property(quint32 index, const QScriptValue::ResolveFlags &resolveMode) const
{
    return QScriptEnginePrivate::property(engine->currentFrame, jscValue, index, resolveMode);
}

inline JSC::JSValue QScriptValuePrivate::property(const JSC::UString &name, const QScriptValue::ResolveFlags &resolveMode) const
{
    JSC::ExecState *exec = engine->currentFrame;
    return QScriptEnginePrivate::property(exec, jscValue, JSC::Identifier(exec, name), resolveMode);
}

inline QScriptValue::PropertyFlags QScriptValuePrivate::propertyFlags(
        const JSC::Identifier &id, const QScriptValue::ResolveFlags &mode) const
{
    return QScriptEnginePrivate::propertyFlags(engine->currentFrame, jscValue, id, mode);
}

inline void QScriptValuePrivate::setProperty(const JSC::Identifier &id, const JSC::JSValue &value,
                                             const QScriptValue::PropertyFlags &flags)
{
    QScriptEnginePrivate::setProperty(engine->currentFrame, jscValue, id, value, flags);
}

inline void QScriptValuePrivate::setProperty(quint32 index, const JSC::JSValue &value,
                                             const QScriptValue::PropertyFlags &flags)
{
    QScriptEnginePrivate::setProperty(engine->currentFrame, jscValue, index, value, flags);
}

inline void QScriptValuePrivate::setProperty(const JSC::UString &name, const JSC::JSValue &value,
                                             const QScriptValue::PropertyFlags &flags)
{
    JSC::ExecState *exec = engine->currentFrame;
    QScriptEnginePrivate::setProperty(exec, jscValue, JSC::Identifier(exec, name), value, flags);
}

inline void* QScriptValuePrivate::operator new(size_t size, QScriptEnginePrivate *engine)
{
    if (engine)
        return engine->allocateScriptValuePrivate(size);
    return malloc(size);
}

inline void QScriptValuePrivate::operator delete(void *ptr)
{
    QScriptValuePrivate *d = reinterpret_cast<QScriptValuePrivate*>(ptr);
    if (d->engine)
        d->engine->freeScriptValuePrivate(d);
    else
        free(d);
}

inline void QScriptEnginePrivate::saveException(JSC::ExecState *exec, JSC::JSValue *val)
{
    if (exec) {
        *val = exec->exception();
        exec->clearException();
    } else {
        *val = JSC::JSValue();
    }
}

inline void QScriptEnginePrivate::restoreException(JSC::ExecState *exec, JSC::JSValue val)
{
    if (exec && val)
        exec->setException(val);
}

inline void QScriptEnginePrivate::registerScriptString(QScriptStringPrivate *value)
{
    Q_ASSERT(value->type == QScriptStringPrivate::HeapAllocated);
    value->prev = 0;
    value->next = registeredScriptStrings;
    if (registeredScriptStrings)
        registeredScriptStrings->prev = value;
    registeredScriptStrings = value;
}

inline void QScriptEnginePrivate::unregisterScriptString(QScriptStringPrivate *value)
{
    Q_ASSERT(value->type == QScriptStringPrivate::HeapAllocated);
    if (value->prev)
        value->prev->next = value->next;
    if (value->next)
        value->next->prev = value->prev;
    if (value == registeredScriptStrings)
        registeredScriptStrings = value->next;
    value->prev = 0;
    value->next = 0;
}

inline QScriptContext *QScriptEnginePrivate::contextForFrame(JSC::ExecState *frame)
{
    if (frame && frame->callerFrame()->hasHostCallFrameFlag() && !frame->callee()
        && frame->callerFrame()->removeHostCallFrameFlag() == QScript::scriptEngineFromExec(frame)->globalExec()) {
        //skip the "fake" context created in Interpreter::execute.
        frame = frame->callerFrame()->removeHostCallFrameFlag();
    }
    return reinterpret_cast<QScriptContext *>(frame);
}

inline JSC::ExecState *QScriptEnginePrivate::frameForContext(QScriptContext *context)
{
    return reinterpret_cast<JSC::ExecState*>(context);
}

inline const JSC::ExecState *QScriptEnginePrivate::frameForContext(const QScriptContext *context)
{
    return reinterpret_cast<const JSC::ExecState*>(context);
}

inline bool QScriptEnginePrivate::hasValidCodeBlockRegister(JSC::ExecState *frame)
{
#if ENABLE(JIT)
    // Frames created by the VM don't have their CodeBlock register
    // initialized. We can detect such frames by checking if the
    // callee is a host JSFunction.
    JSC::JSObject *callee = frame->callee();
    return !(callee && callee->inherits(&JSC::JSFunction::info)
             && JSC::asFunction(callee)->isHostFunction());
#else
    Q_UNUSED(frame);
    return true;
#endif
}

inline JSC::ExecState *QScriptEnginePrivate::globalExec() const
{
    return originalGlobalObject()->globalExec();
}

inline JSC::JSValue QScriptEnginePrivate::newArray(JSC::ExecState *exec, uint length)
{
    return JSC::constructEmptyArray(exec, length);
}

inline JSC::JSValue QScriptEnginePrivate::newDate(JSC::ExecState *exec, qsreal value)
{
    JSC::JSValue val = JSC::jsNumber(exec, value);
    JSC::ArgList args(&val, 1);
    return JSC::constructDate(exec, args);
}

inline JSC::JSValue QScriptEnginePrivate::newDate(JSC::ExecState *exec, const QDateTime &value)
{
    return newDate(exec, QScript::DateTimeToMs(exec, value));
}

inline JSC::JSValue QScriptEnginePrivate::newObject()
{
    return new (currentFrame)QScriptObject(scriptObjectStructure);
}

inline bool QScriptEnginePrivate::isObject(JSC::JSValue value)
{
    return value && value.isObject();
}

inline bool QScriptEnginePrivate::isArray(JSC::JSValue value)
{
    return isObject(value) && value.inherits(&JSC::JSArray::info);
}

inline bool QScriptEnginePrivate::isDate(JSC::JSValue value)
{
    return isObject(value) && value.inherits(&JSC::DateInstance::info);
}

inline bool QScriptEnginePrivate::isError(JSC::JSValue value)
{
    return isObject(value) && value.inherits(&JSC::ErrorInstance::info);
}

inline bool QScriptEnginePrivate::isRegExp(JSC::JSValue value)
{
    return isObject(value) && value.inherits(&JSC::RegExpObject::info);
}

inline bool QScriptEnginePrivate::isVariant(JSC::JSValue value)
{
    if (!isObject(value) || !value.inherits(&QScriptObject::info))
        return false;
    QScriptObject *object = static_cast<QScriptObject*>(JSC::asObject(value));
    QScriptObjectDelegate *delegate = object->delegate();
    return (delegate && (delegate->type() == QScriptObjectDelegate::Variant));
}

inline bool QScriptEnginePrivate::isQObject(JSC::JSValue value)
{
#ifndef QT_NO_QOBJECT
    if (!isObject(value) || !value.inherits(&QScriptObject::info))
        return false;
    QScriptObject *object = static_cast<QScriptObject*>(JSC::asObject(value));
    QScriptObjectDelegate *delegate = object->delegate();

    if (delegate) {
        if (delegate->type() == QScriptObjectDelegate::QtObject
            || (delegate->type() == QScriptObjectDelegate::DeclarativeClassObject
            && static_cast<QScript::DeclarativeObjectDelegate*>(delegate)->scriptClass()->isQObject()))
            return true;

        if (delegate->type() == QScriptObjectDelegate::Variant) {
            QVariant var = variantValue(value);
            int type = var.userType();
            if ((QMetaType::typeFlags(type) & QMetaType::PointerToQObject))
                return true;
        }
    }
#endif
    return false;
}

inline bool QScriptEnginePrivate::isQMetaObject(JSC::JSValue value)
{
#ifndef QT_NO_QOBJECT
    return isObject(value) && JSC::asObject(value)->inherits(&QScript::QMetaObjectWrapperObject::info);
#else
    return false;
#endif
}

inline bool QScriptEnginePrivate::toBool(JSC::ExecState *exec, JSC::JSValue value)
{
    JSC::JSValue savedException;
    saveException(exec, &savedException);
    bool result = value.toBoolean(exec);
    restoreException(exec, savedException);
    return result;
}

inline qsreal QScriptEnginePrivate::toInteger(JSC::ExecState *exec, JSC::JSValue value)
{
    JSC::JSValue savedException;
    saveException(exec, &savedException);
    qsreal result = value.toInteger(exec);
    restoreException(exec, savedException);
    return result;
}

inline qsreal QScriptEnginePrivate::toNumber(JSC::ExecState *exec, JSC::JSValue value)
{
    JSC::JSValue savedException;
    saveException(exec, &savedException);
    qsreal result = value.toNumber(exec);
    restoreException(exec, savedException);
    return result;
}

inline qint32 QScriptEnginePrivate::toInt32(JSC::ExecState *exec, JSC::JSValue value)
{
    JSC::JSValue savedException;
    saveException(exec, &savedException);
    qint32 result = value.toInt32(exec);
    restoreException(exec, savedException);
    return result;
}

inline quint32 QScriptEnginePrivate::toUInt32(JSC::ExecState *exec, JSC::JSValue value)
{
    JSC::JSValue savedException;
    saveException(exec, &savedException);
    quint32 result = value.toUInt32(exec);
    restoreException(exec, savedException);
    return result;
}

inline quint16 QScriptEnginePrivate::toUInt16(JSC::ExecState *exec, JSC::JSValue value)
{
    // ### no equivalent function in JSC
    return QScript::ToUInt16(toNumber(exec, value));
}

inline JSC::UString QScriptEnginePrivate::toString(JSC::ExecState *exec, JSC::JSValue value)
{
    if (!value)
        return JSC::UString();
    JSC::JSValue savedException;
    saveException(exec, &savedException);
    JSC::UString str = value.toString(exec);
    if (exec && exec->hadException() && !str.size()) {
        JSC::JSValue savedException2;
        saveException(exec, &savedException2);
        str = savedException2.toString(exec);
        restoreException(exec, savedException2);
    }
    if (savedException)
        restoreException(exec, savedException);
    return str;
}

inline QDateTime QScriptEnginePrivate::toDateTime(JSC::ExecState *exec, JSC::JSValue value)
{
    if (!isDate(value))
        return QDateTime();
    qsreal t = static_cast<JSC::DateInstance*>(JSC::asObject(value))->internalNumber();
    return QScript::MsToDateTime(exec, t);
}

inline QObject *QScriptEnginePrivate::toQObject(JSC::ExecState *exec, JSC::JSValue value)
{
#ifndef QT_NO_QOBJECT
    if (isObject(value) && value.inherits(&QScriptObject::info)) {
        QScriptObject *object = static_cast<QScriptObject*>(JSC::asObject(value));
        QScriptObjectDelegate *delegate = object->delegate();
        if (!delegate)
            return 0;
        if (delegate->type() == QScriptObjectDelegate::QtObject)
            return static_cast<QScript::QObjectDelegate*>(delegate)->value();
        if (delegate->type() == QScriptObjectDelegate::DeclarativeClassObject)
            return static_cast<QScript::DeclarativeObjectDelegate*>(delegate)->scriptClass()->toQObject(declarativeObject(value));
        if (delegate->type() == QScriptObjectDelegate::Variant) {
            QVariant var = variantValue(value);
            int type = var.userType();
            if (QMetaType::typeFlags(type) & QMetaType::PointerToQObject)
                return *reinterpret_cast<QObject* const *>(var.constData());
        }
    } else if (isObject(value) && value.inherits(&QScript::QScriptActivationObject::info)) {
        QScript::QScriptActivationObject *proxy = static_cast<QScript::QScriptActivationObject *>(JSC::asObject(value));
        return toQObject(exec, proxy->delegate());
    }
#endif
    return 0;
}

inline const QMetaObject *QScriptEnginePrivate::toQMetaObject(JSC::ExecState*, JSC::JSValue value)
{
#ifndef QT_NO_QOBJECT
    if (isQMetaObject(value))
        return static_cast<QScript::QMetaObjectWrapperObject*>(JSC::asObject(value))->value();
#endif
    return 0;
}

inline QVariant &QScriptEnginePrivate::variantValue(JSC::JSValue value)
{
    Q_ASSERT(value.inherits(&QScriptObject::info));
    QScriptObjectDelegate *delegate = static_cast<QScriptObject*>(JSC::asObject(value))->delegate();
    Q_ASSERT(delegate && (delegate->type() == QScriptObjectDelegate::Variant));
    return static_cast<QScript::QVariantDelegate*>(delegate)->value();
}

inline void QScriptEnginePrivate::setVariantValue(JSC::JSValue objectValue, const QVariant &value)
{
    Q_ASSERT(objectValue.inherits(&QScriptObject::info));
    QScriptObjectDelegate *delegate = static_cast<QScriptObject*>(JSC::asObject(objectValue))->delegate();
    Q_ASSERT(delegate && (delegate->type() == QScriptObjectDelegate::Variant));
    static_cast<QScript::QVariantDelegate*>(delegate)->setValue(value);
}

inline QScriptDeclarativeClass *QScriptEnginePrivate::declarativeClass(JSC::JSValue v)
{
    if (!QScriptEnginePrivate::isObject(v) || !v.inherits(&QScriptObject::info))
        return 0;
    QScriptObject *scriptObject = static_cast<QScriptObject*>(JSC::asObject(v));
    QScriptObjectDelegate *delegate = scriptObject->delegate();
    if (!delegate || (delegate->type() != QScriptObjectDelegate::DeclarativeClassObject))
        return 0;
    return static_cast<QScript::DeclarativeObjectDelegate*>(delegate)->scriptClass();
}

inline QScriptDeclarativeClass::Object *QScriptEnginePrivate::declarativeObject(JSC::JSValue v)
{
    if (!QScriptEnginePrivate::isObject(v) || !v.inherits(&QScriptObject::info))
        return 0;
    QScriptObject *scriptObject = static_cast<QScriptObject*>(JSC::asObject(v));
    QScriptObjectDelegate *delegate = scriptObject->delegate();
    if (!delegate || (delegate->type() != QScriptObjectDelegate::DeclarativeClassObject))
        return 0;
    return static_cast<QScript::DeclarativeObjectDelegate*>(delegate)->object();
}

QT_END_NAMESPACE

#endif
