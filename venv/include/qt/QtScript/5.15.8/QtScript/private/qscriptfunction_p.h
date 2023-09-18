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

#ifndef QSCRIPTFUNCTION_P_H
#define QSCRIPTFUNCTION_P_H

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

#include <QtCore/qglobal.h>

#include "qscriptengine.h"

#include "PrototypeFunction.h"

QT_BEGIN_NAMESPACE

namespace QScript
{

class FunctionWrapper : public JSC::PrototypeFunction // ### subclass InternalFunction instead
{
public:
    // work around CELL_SIZE limitation
    struct Data
    {
        QScriptEngine::FunctionSignature function;
    };

    FunctionWrapper(JSC::ExecState*, int length, const JSC::Identifier&,
                    QScriptEngine::FunctionSignature);
    ~FunctionWrapper();

    virtual const JSC::ClassInfo* classInfo() const { return &info; }
    static const JSC::ClassInfo info;

    QScriptEngine::FunctionSignature function() const
    { return data->function; }

private:
    virtual JSC::ConstructType getConstructData(JSC::ConstructData&);

    static JSC::JSValue JSC_HOST_CALL proxyCall(JSC::ExecState *, JSC::JSObject *,
                                                JSC::JSValue, const JSC::ArgList &);
    static JSC::JSObject* proxyConstruct(JSC::ExecState *, JSC::JSObject *,
                                         const JSC::ArgList &);

private:
    Data *data;
};

class FunctionWithArgWrapper : public JSC::PrototypeFunction
{
public:
    // work around CELL_SIZE limitation
    struct Data
    {
        QScriptEngine::FunctionWithArgSignature function;
        void *arg;
    };

    FunctionWithArgWrapper(JSC::ExecState*, int length, const JSC::Identifier&,
                           QScriptEngine::FunctionWithArgSignature, void *);
    ~FunctionWithArgWrapper();

    virtual const JSC::ClassInfo* classInfo() const { return &info; }
    static const JSC::ClassInfo info;

    QScriptEngine::FunctionWithArgSignature function() const
    { return data->function; }

    void *arg() const
    { return data->arg; }

private:
    virtual JSC::ConstructType getConstructData(JSC::ConstructData&);

    static JSC::JSValue JSC_HOST_CALL proxyCall(JSC::ExecState *, JSC::JSObject *,
                                                JSC::JSValue , const JSC::ArgList &);
    static JSC::JSObject* proxyConstruct(JSC::ExecState *, JSC::JSObject *,
                                         const JSC::ArgList &);

private:
    Data *data;
};

} // namespace QScript

QT_END_NAMESPACE

#endif
