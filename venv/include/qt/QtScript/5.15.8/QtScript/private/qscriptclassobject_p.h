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

#ifndef QSCRIPTCLASSOBJECT_P_H
#define QSCRIPTCLASSOBJECT_P_H

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

#include "qscriptobject_p.h"

QT_BEGIN_NAMESPACE

class QScriptClass;

namespace QScript
{

class ClassObjectDelegate : public QScriptObjectDelegate
{
public:
    ClassObjectDelegate(QScriptClass *scriptClass);
    ~ClassObjectDelegate();

    inline QScriptClass *scriptClass() const;
    inline void setScriptClass(QScriptClass *scriptClass);

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

    virtual JSC::CallType getCallData(QScriptObject*, JSC::CallData&);
    static JSC::JSValue JSC_HOST_CALL call(JSC::ExecState*, JSC::JSObject*,
                                           JSC::JSValue, const JSC::ArgList&);
    virtual JSC::ConstructType getConstructData(QScriptObject*, JSC::ConstructData&);
    static JSC::JSObject* construct(JSC::ExecState*, JSC::JSObject*,
                                    const JSC::ArgList&);

    virtual bool hasInstance(QScriptObject*, JSC::ExecState*,
                             JSC::JSValue value, JSC::JSValue proto);

private:
    QScriptClass *m_scriptClass;
};

inline QScriptClass *ClassObjectDelegate::scriptClass() const
{
    return m_scriptClass;
}

inline void ClassObjectDelegate::setScriptClass(QScriptClass *scriptClass)
{
    Q_ASSERT(scriptClass != 0);
    m_scriptClass = scriptClass;
}

} // namespace QScript

QT_END_NAMESPACE

#endif
