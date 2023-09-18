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

#ifndef QSCRIPTVALUE_H
#define QSCRIPTVALUE_H

#include <QtCore/qstring.h>

#include <QtCore/qlist.h>
#include <QtCore/qsharedpointer.h>
#include <QtScript/qtscriptglobal.h>

QT_BEGIN_NAMESPACE


class QScriptClass;
class QScriptValue;
class QScriptEngine;
class QScriptString;
class QVariant;
class QObject;
struct QMetaObject;
class QDateTime;
#ifndef QT_NO_REGEXP
class QRegExp;
#endif

typedef QList<QScriptValue> QScriptValueList;

typedef double qsreal;

class QScriptValuePrivate;
class QScriptEnginePrivate;
struct QScriptValuePrivatePointerDeleter;
class Q_SCRIPT_EXPORT QScriptValue
{
public:
    enum ResolveFlag {
        ResolveLocal        = 0x00,
        ResolvePrototype    = 0x01,
        ResolveScope        = 0x02,
        ResolveFull         = ResolvePrototype | ResolveScope
    };

    Q_DECLARE_FLAGS(ResolveFlags, ResolveFlag)

    enum PropertyFlag {
        ReadOnly            = 0x00000001,
        Undeletable         = 0x00000002,
        SkipInEnumeration   = 0x00000004,

        PropertyGetter      = 0x00000008,
        PropertySetter      = 0x00000010,

        QObjectMember       = 0x00000020,

        KeepExistingFlags   = 0x00000800,

        UserRange           = 0xff000000            // Users may use these as they see fit.
    };
    Q_DECLARE_FLAGS(PropertyFlags, PropertyFlag)

    enum SpecialValue {
        NullValue,
        UndefinedValue
    };

public:
    QScriptValue();
    ~QScriptValue();
    QScriptValue(const QScriptValue &other);
    QScriptValue(QScriptEngine *engine, SpecialValue val);
    QScriptValue(QScriptEngine *engine, bool val);
    QScriptValue(QScriptEngine *engine, int val);
    QScriptValue(QScriptEngine *engine, uint val);
    QScriptValue(QScriptEngine *engine, qsreal val);
    QScriptValue(QScriptEngine *engine, const QString &val);
#ifndef QT_NO_CAST_FROM_ASCII
    QT_ASCII_CAST_WARN QScriptValue(QScriptEngine *engine, const char *val);
#endif

    QScriptValue(SpecialValue value);
    QScriptValue(bool value);
    QScriptValue(int value);
    QScriptValue(uint value);
    QScriptValue(qsreal value);
    QScriptValue(const QString &value);
    QScriptValue(const QLatin1String &value);
#ifndef QT_NO_CAST_FROM_ASCII
    QT_ASCII_CAST_WARN QScriptValue(const char *value);
#endif

    QScriptValue &operator=(const QScriptValue &other);

    QScriptEngine *engine() const;

    bool isValid() const;
    bool isBool() const;
    bool isBoolean() const;
    bool isNumber() const;
    bool isFunction() const;
    bool isNull() const;
    bool isString() const;
    bool isUndefined() const;
    bool isVariant() const;
    bool isQObject() const;
    bool isQMetaObject() const;
    bool isObject() const;
    bool isDate() const;
    bool isRegExp() const;
    bool isArray() const;
    bool isError() const;

    QString toString() const;
    qsreal toNumber() const;
    bool toBool() const;
    bool toBoolean() const;
    qsreal toInteger() const;
    qint32 toInt32() const;
    quint32 toUInt32() const;
    quint16 toUInt16() const;
    QVariant toVariant() const;
    QObject *toQObject() const;
    const QMetaObject *toQMetaObject() const;
    QScriptValue toObject() const;
    QDateTime toDateTime() const;
#ifndef QT_NO_REGEXP
    QRegExp toRegExp() const;
#endif

    bool instanceOf(const QScriptValue &other) const;

    bool lessThan(const QScriptValue &other) const;
    bool equals(const QScriptValue &other) const;
    bool strictlyEquals(const QScriptValue &other) const;

    QScriptValue prototype() const;
    void setPrototype(const QScriptValue &prototype);

    QScriptValue scope() const;
    void setScope(const QScriptValue &scope);

    QScriptValue property(const QString &name,
                          const ResolveFlags &mode = ResolvePrototype) const;
    void setProperty(const QString &name, const QScriptValue &value,
                     const PropertyFlags &flags = KeepExistingFlags);

    QScriptValue property(quint32 arrayIndex,
                          const ResolveFlags &mode = ResolvePrototype) const;
    void setProperty(quint32 arrayIndex, const QScriptValue &value,
                     const PropertyFlags &flags = KeepExistingFlags);

    QScriptValue property(const QScriptString &name,
                          const ResolveFlags &mode = ResolvePrototype) const;
    void setProperty(const QScriptString &name, const QScriptValue &value,
                     const PropertyFlags &flags = KeepExistingFlags);

    QScriptValue::PropertyFlags propertyFlags(
        const QString &name, const ResolveFlags &mode = ResolvePrototype) const;
    QScriptValue::PropertyFlags propertyFlags(
        const QScriptString &name, const ResolveFlags &mode = ResolvePrototype) const;

    QScriptValue call(const QScriptValue &thisObject = QScriptValue(),
                      const QScriptValueList &args = QScriptValueList());
    QScriptValue call(const QScriptValue &thisObject,
                      const QScriptValue &arguments);
    QScriptValue construct(const QScriptValueList &args = QScriptValueList());
    QScriptValue construct(const QScriptValue &arguments);

    QScriptValue data() const;
    void setData(const QScriptValue &data);

    QScriptClass *scriptClass() const;
    void setScriptClass(QScriptClass *scriptClass);

    qint64 objectId() const;

private:
    // force compile error, prevent QScriptValue(bool) to be called
    QScriptValue(void *);
    // force compile error, prevent QScriptValue(QScriptEngine*, bool) to be called
    QScriptValue(QScriptEngine *, void *);

    QScriptValue(QScriptValuePrivate*);

private:
    QExplicitlySharedDataPointer<QScriptValuePrivate> d_ptr;

    Q_DECLARE_PRIVATE(QScriptValue)

    friend class QScriptEnginePrivate;
};

Q_DECLARE_OPERATORS_FOR_FLAGS(QScriptValue::ResolveFlags)
Q_DECLARE_OPERATORS_FOR_FLAGS(QScriptValue::PropertyFlags)

QT_END_NAMESPACE

#endif
