/****************************************************************************
**
** Copyright (C) 2017 Ford Motor Company
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtRemoteObjects module of the Qt Toolkit.
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

#ifndef QREMOTEOBJECTSOURCE_H
#define QREMOTEOBJECTSOURCE_H

#include <QtCore/qscopedpointer.h>
#include <QtRemoteObjects/qtremoteobjectglobal.h>
#include <QtCore/qmetaobject.h>

QT_BEGIN_NAMESPACE

namespace QtPrivate {

//Based on compile time checks for static connect() from qobjectdefs_impl.h
template <class ObjectType, typename Func1, typename Func2>
static inline int qtro_property_index(Func1, Func2, const char *propName)
{
    typedef QtPrivate::FunctionPointer<Func1> Type1;
    typedef QtPrivate::FunctionPointer<Func2> Type2;

    //compilation error if the arguments do not match.
    Q_STATIC_ASSERT_X(int(Type1::ArgumentCount) >= int(Type2::ArgumentCount),
                      "Argument counts are not compatible.");
    Q_STATIC_ASSERT_X((QtPrivate::CheckCompatibleArguments<typename Type1::Arguments, typename Type2::Arguments>::value),
                      "Arguments are not compatible.");
    Q_STATIC_ASSERT_X((QtPrivate::AreArgumentsCompatible<typename Type1::ReturnType, typename Type2::ReturnType>::value),
                      "Return types are not compatible.");
    return ObjectType::staticMetaObject.indexOfProperty(propName);
}

template <class ObjectType, typename Func1, typename Func2>
static inline int qtro_signal_index(Func1 func, Func2, int *count, int const **types)
{
    typedef QtPrivate::FunctionPointer<Func1> Type1;
    typedef QtPrivate::FunctionPointer<Func2> Type2;

    //compilation error if the arguments do not match.
    Q_STATIC_ASSERT_X(int(Type1::ArgumentCount) >= int(Type2::ArgumentCount),
                      "Argument counts are not compatible.");
    Q_STATIC_ASSERT_X((QtPrivate::CheckCompatibleArguments<typename Type1::Arguments, typename Type2::Arguments>::value),
                      "Arguments are not compatible.");
    Q_STATIC_ASSERT_X((QtPrivate::AreArgumentsCompatible<typename Type1::ReturnType, typename Type2::ReturnType>::value),
                      "Return types are not compatible.");
    const QMetaMethod sig = QMetaMethod::fromSignal(func);
    *count = Type2::ArgumentCount;
    *types = QtPrivate::ConnectionTypes<typename Type2::Arguments>::types();
    return sig.methodIndex();
}

template <class ObjectType, typename Func1, typename Func2>
static inline void qtro_method_test(Func1, Func2)
{
    typedef QtPrivate::FunctionPointer<Func1> Type1;
    typedef QtPrivate::FunctionPointer<Func2> Type2;

    //compilation error if the arguments do not match.
    Q_STATIC_ASSERT_X(int(Type1::ArgumentCount) >= int(Type2::ArgumentCount),
                      "Argument counts are not compatible.");
    Q_STATIC_ASSERT_X((QtPrivate::CheckCompatibleArguments<typename Type1::Arguments, typename Type2::Arguments>::value),
                      "Arguments are not compatible.");
    Q_STATIC_ASSERT_X((QtPrivate::AreArgumentsCompatible<typename Type1::ReturnType, typename Type2::ReturnType>::value),
                      "Return types are not compatible.");
}

// The stringData, methodMatch and QMetaObjectPrivate methods are modified versions of the code
// from qmetaobject_p.h/qmetaobject.cpp.  The modifications are based on our custom need to match
// a method name that comes from the .rep file.
// The QMetaObjectPrivate struct should only have members appended to maintain binary compatibility,
// so we should be fine with only the listed version with the fields we use.
inline const QByteArray apiStringData(const QMetaObject *mo, int index)
{
    const QByteArrayDataPtr data = { const_cast<QByteArrayData*>(&mo->d.stringdata[index]) };
    return data;
}

inline bool apiMethodMatch(const QMetaObject *m, int handle,
                        const QByteArray &name, int argc,
                        const int *types)
{
    if (int(m->d.data[handle + 1]) != argc)
        return false;
    if (apiStringData(m, m->d.data[handle]) != name)
        return false;
    int paramsIndex = m->d.data[handle + 2] + 1;
    for (int i = 0; i < argc; ++i) {
        uint typeInfo = m->d.data[paramsIndex + i];
        if (typeInfo & 0x80000000) { // Custom/named type, compare names
            const char *t = QMetaType::typeName(types[i]);
            const auto type = QByteArray::fromRawData(t, qstrlen(t));
            if (type != apiStringData(m, typeInfo & 0x7FFFFFFF))
                return false;
        } else if (types[i] != int(typeInfo))
            return false;
    }
    return true;
}

struct QMetaObjectPrivate
{
    // revision 7 is Qt 5.0 everything lower is not supported
    // revision 8 is Qt 5.12: It adds the enum name to QMetaEnum
    enum { OutputRevision = 8 }; // Used by moc, qmetaobjectbuilder and qdbus

    int revision;
    int className;
    int classInfoCount, classInfoData;
    int methodCount, methodData;
    int propertyCount, propertyData;
    int enumeratorCount, enumeratorData;
    int constructorCount, constructorData;
    int flags;
    int signalCount;
};

template <class ObjectType, typename Func1, typename Func2>
static inline int qtro_method_index(Func1, Func2, const char *methodName, int *count, int const **types)
{
    typedef QtPrivate::FunctionPointer<Func1> Type1;
    typedef QtPrivate::FunctionPointer<Func2> Type2;

    //compilation error if the arguments do not match.
    Q_STATIC_ASSERT_X(int(Type1::ArgumentCount) >= int(Type2::ArgumentCount),
                      "Argument counts are not compatible.");
    Q_STATIC_ASSERT_X((QtPrivate::CheckCompatibleArguments<typename Type1::Arguments, typename Type2::Arguments>::value),
                      "Arguments are not compatible.");
    Q_STATIC_ASSERT_X((QtPrivate::AreArgumentsCompatible<typename Type1::ReturnType, typename Type2::ReturnType>::value),
                      "Return types are not compatible.");
    *count = Type2::ArgumentCount;
    *types = QtPrivate::ConnectionTypes<typename Type2::Arguments>::types();

    int result = ObjectType::staticMetaObject.indexOfMethod(methodName);
    if (result >= 0)
        return result;
    // We can have issues, specifically with enums, since the compiler can infer the class.  Since
    // indexOfMethod() is doing string comparisons for registered types, "MyEnum" and "MyClass::MyEnum"
    // won't match.
    // Below is similar to QMetaObject->indexOfMethod, but template magic has already matched parameter
    // types, so we need to find a match for the API method name + parameters.  Neither approach works
    // 100%, as the below code doesn't match a parameter of type "size_t" (which the template match
    // identifies as "ulong").  These subtleties can cause the below string comparison fails.
    // There is no known case that would fail both methods.
    // TODO: is there a way to make this a constexpr so a failure is detected at compile time?
    int nameLength = strchr(methodName, '(') - methodName;
    const auto name = QByteArray::fromRawData(methodName, nameLength);
    for (const QMetaObject *m = &ObjectType::staticMetaObject; m; m = m->d.superdata) {
        const auto priv = reinterpret_cast<const QMetaObjectPrivate*>(m->d.data);
        int i = (priv->methodCount - 1);
        const int end = priv->signalCount;
        for (; i >= end; --i) {
            int handle = priv->methodData + 5*i;
            if (apiMethodMatch(m, handle, name, *count, *types))
                return i + m->methodOffset();
        }
    }
    qWarning() << "No matching method for" << methodName << "in the provided metaclass" << ObjectType::staticMetaObject.className();
    return -1;
}

template <class ObjectType>
static inline QByteArray qtro_enum_signature(const char *enumName)
{
    const auto qme = ObjectType::staticMetaObject.enumerator(ObjectType::staticMetaObject.indexOfEnumerator(enumName));
    return QByteArrayLiteral("1::2").replace("1", qme.scope()).replace("2", qme.name());
}

QByteArray qtro_classinfo_signature(const QMetaObject *metaObject);

}

// TODO ModelInfo just needs roles, and no need for SubclassInfo
class QAbstractItemModel;

struct ModelInfo
{
    QAbstractItemModel *ptr;
    QString name;
    QByteArray roles;
};

class SourceApiMap
{
protected:
    SourceApiMap() {}
public:
    virtual ~SourceApiMap() {}
    virtual QString name() const = 0;
    virtual QString typeName() const = 0;
    virtual QByteArray className() const { return typeName().toLatin1().append("Source"); }
    virtual int enumCount() const = 0;
    virtual int propertyCount() const = 0;
    virtual int signalCount() const = 0;
    virtual int methodCount() const = 0;
    virtual int sourceEnumIndex(int index) const = 0;
    virtual int sourcePropertyIndex(int index) const = 0;
    virtual int sourceSignalIndex(int index) const = 0;
    virtual int sourceMethodIndex(int index) const = 0;
    virtual int signalParameterCount(int index) const = 0;
    virtual int signalParameterType(int sigIndex, int paramIndex) const = 0;
    virtual const QByteArray signalSignature(int index) const = 0;
    virtual QList<QByteArray> signalParameterNames(int index) const = 0;
    virtual int methodParameterCount(int index) const = 0;
    virtual int methodParameterType(int methodIndex, int paramIndex) const = 0;
    virtual const QByteArray methodSignature(int index) const = 0;
    virtual QMetaMethod::MethodType methodType(int index) const = 0;
    virtual const QByteArray typeName(int index) const = 0;
    virtual QList<QByteArray> methodParameterNames(int index) const = 0;
    virtual int propertyIndexFromSignal(int index) const = 0;
    virtual int propertyRawIndexFromSignal(int index) const = 0;
    virtual QByteArray objectSignature() const = 0;
    virtual bool isDynamic() const { return false; }
    virtual bool isAdapterSignal(int) const { return false; }
    virtual bool isAdapterMethod(int) const { return false; }
    virtual bool isAdapterProperty(int) const { return false; }
    QVector<ModelInfo> m_models;
    QVector<SourceApiMap *> m_subclasses;
};

QT_END_NAMESPACE

#endif
