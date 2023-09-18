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

#ifndef QSCRIPTDEBUGGERVALUE_P_H
#define QSCRIPTDEBUGGERVALUE_P_H

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
#include <QtCore/private/qscopedpointer_p.h>
#include <QtCore/qlist.h>

QT_BEGIN_NAMESPACE

class QScriptValue;
class QScriptEngine;
class QDataStream;

class QScriptDebuggerValuePrivate;
class Q_AUTOTEST_EXPORT QScriptDebuggerValue
{
public:
    enum ValueType {
        NoValue,
        UndefinedValue,
        NullValue,
        BooleanValue,
        StringValue,
        NumberValue,
        ObjectValue
    };

    QScriptDebuggerValue();
    QScriptDebuggerValue(const QScriptValue &value);
    QScriptDebuggerValue(double value);
    QScriptDebuggerValue(bool value);
    QScriptDebuggerValue(const QString &value);
    QScriptDebuggerValue(qint64 objectId);
    QScriptDebuggerValue(ValueType type);
    QScriptDebuggerValue(const QScriptDebuggerValue &other);
    ~QScriptDebuggerValue();

    QScriptDebuggerValue &operator=(const QScriptDebuggerValue &other);

    ValueType type() const;

    double numberValue() const;
    bool booleanValue() const;
    QString stringValue() const;
    qint64 objectId() const;

    QScriptValue toScriptValue(QScriptEngine *engine) const;
    QString toString() const;

    bool operator==(const QScriptDebuggerValue &other) const;
    bool operator!=(const QScriptDebuggerValue &other) const;

private:
    QScopedSharedPointer<QScriptDebuggerValuePrivate> d_ptr;

    Q_DECLARE_PRIVATE(QScriptDebuggerValue)
};

typedef QList<QScriptDebuggerValue> QScriptDebuggerValueList;

Q_AUTOTEST_EXPORT QDataStream &operator<<(QDataStream &, const QScriptDebuggerValue &);
Q_AUTOTEST_EXPORT QDataStream &operator>>(QDataStream &, QScriptDebuggerValue &);

QT_END_NAMESPACE

#endif
