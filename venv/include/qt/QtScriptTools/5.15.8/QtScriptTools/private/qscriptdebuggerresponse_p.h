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

#ifndef QSCRIPTDEBUGGERRESPONSE_P_H
#define QSCRIPTDEBUGGERRESPONSE_P_H

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
#include <QtCore/qscopedpointer.h>
#include <QtCore/qmap.h>
#include <QtCore/qvariant.h>

#include "qscriptbreakpointdata_p.h"
#include "qscriptscriptdata_p.h"
#include "qscriptdebuggervalue_p.h"
#include "qscriptdebuggervalueproperty_p.h"

QT_BEGIN_NAMESPACE

class QDataStream;
class QScriptContextInfo;

class QScriptDebuggerResponsePrivate;
class Q_AUTOTEST_EXPORT QScriptDebuggerResponse
{
public:
    friend Q_AUTOTEST_EXPORT QDataStream &operator<<(QDataStream &, const QScriptDebuggerResponse &);
    friend Q_AUTOTEST_EXPORT QDataStream &operator>>(QDataStream &, QScriptDebuggerResponse &);

    enum Error {
        NoError,
        InvalidContextIndex,
        InvalidArgumentIndex,
        InvalidScriptID,
        InvalidBreakpointID,
        UserError = 1000,
        MaxUserError = 32767
    };

    QScriptDebuggerResponse();
    QScriptDebuggerResponse(const QScriptDebuggerResponse &other);
    ~QScriptDebuggerResponse();

    Error error() const;
    void setError(Error error);

    QVariant result() const;
    void setResult(const QVariant &value);

    void setResult(int value);
    void setResult(const QString &value);
    void setResult(const QScriptBreakpointData &data);
    void setResult(const QScriptBreakpointMap &breakpoints);
    void setResult(const QScriptScriptMap &scripts);
    void setResult(const QScriptScriptData &data);
    void setResult(const QScriptDebuggerValue &value);
    void setResult(const QScriptDebuggerValueList &value);
    void setResult(const QScriptDebuggerValuePropertyList &props);
    void setResult(const QScriptContextInfo &info);

    int resultAsInt() const;
    qint64 resultAsLongLong() const;
    QString resultAsString() const;
    QScriptBreakpointData resultAsBreakpointData() const;
    QScriptBreakpointMap resultAsBreakpoints() const;
    QScriptScriptMap resultAsScripts() const;
    QScriptScriptData resultAsScriptData() const;
    QScriptDebuggerValue resultAsScriptValue() const;
    QScriptDebuggerValueList resultAsScriptValueList() const;
    QScriptDebuggerValuePropertyList resultAsScriptValuePropertyList() const;
    QScriptContextInfo resultAsContextInfo() const;

    bool async() const;
    void setAsync(bool async);

    QScriptDebuggerResponse &operator=(const QScriptDebuggerResponse &other);

    bool operator==(const QScriptDebuggerResponse &other) const;
    bool operator!=(const QScriptDebuggerResponse &other) const;

private:
    QScopedPointer<QScriptDebuggerResponsePrivate> d_ptr;

    Q_DECLARE_PRIVATE(QScriptDebuggerResponse)
};

Q_AUTOTEST_EXPORT QDataStream &operator<<(QDataStream &, const QScriptDebuggerResponse &);
Q_AUTOTEST_EXPORT QDataStream &operator>>(QDataStream &, QScriptDebuggerResponse &);

QT_END_NAMESPACE

#endif
