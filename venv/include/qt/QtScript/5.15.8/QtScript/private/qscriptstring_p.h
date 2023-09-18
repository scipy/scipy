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

#ifndef QSCRIPTSTRING_P_H
#define QSCRIPTSTRING_P_H

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
#include <QtCore/qshareddata.h>

#include "Identifier.h"

QT_BEGIN_NAMESPACE

class QScriptEnginePrivate;
class QScriptStringPrivate : public QSharedData
{
public:
    enum AllocationType {
        StackAllocated,
        HeapAllocated
    };

    inline QScriptStringPrivate(QScriptEnginePrivate *engine, const JSC::Identifier &id,
                                AllocationType type);
    inline ~QScriptStringPrivate();
    static inline void init(QScriptString &q, QScriptStringPrivate *d);

    static inline QScriptStringPrivate *get(const QScriptString &q);

    inline void detachFromEngine();

    static inline bool isValid(const QScriptString &q);

    QScriptEnginePrivate *engine;
    JSC::Identifier identifier;
    AllocationType type;

    // linked list of engine's script values
    QScriptStringPrivate *prev;
    QScriptStringPrivate *next;
};

inline QScriptStringPrivate::QScriptStringPrivate(QScriptEnginePrivate *e, const JSC::Identifier &id,
                                                  AllocationType tp)
    : engine(e), identifier(id), type(tp), prev(0), next(0)
{
}

inline QScriptStringPrivate::~QScriptStringPrivate()
{
}

inline void QScriptStringPrivate::init(QScriptString &q, QScriptStringPrivate *d)
{
    q.d_ptr = d;
}

inline QScriptStringPrivate *QScriptStringPrivate::get(const QScriptString &q)
{
    return const_cast<QScriptStringPrivate*>(q.d_func());
}

inline void QScriptStringPrivate::detachFromEngine()
{
    engine = 0;
    identifier = JSC::Identifier();
}

inline bool QScriptStringPrivate::isValid(const QScriptString &q)
{
    return (q.d_ptr && q.d_ptr->engine);
}

QT_END_NAMESPACE

#endif
