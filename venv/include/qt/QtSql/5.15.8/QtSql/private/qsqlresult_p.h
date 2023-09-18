/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtSql module of the Qt Toolkit.
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

#ifndef QSQLRESULT_P_H
#define QSQLRESULT_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of qsql*model.h .  This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.
//

#include <QtSql/private/qtsqlglobal_p.h>
#include <QtCore/qpointer.h>
#include "qsqlerror.h"
#include "qsqlresult.h"
#include "qsqldriver.h"

QT_BEGIN_NAMESPACE

// convenience method Q*ResultPrivate::drv_d_func() returns pointer to private driver. Compare to Q_DECLARE_PRIVATE in qglobal.h.
#define Q_DECLARE_SQLDRIVER_PRIVATE(Class) \
    inline const Class##Private* drv_d_func() const { return !sqldriver ? nullptr : reinterpret_cast<const Class *>(static_cast<const QSqlDriver*>(sqldriver))->d_func(); } \
    inline Class##Private* drv_d_func()  { return !sqldriver ? nullptr : reinterpret_cast<Class *>(static_cast<QSqlDriver*>(sqldriver))->d_func(); }

struct QHolder {
    QHolder(const QString &hldr = QString(), int index = -1): holderName(hldr), holderPos(index) { }
    bool operator==(const QHolder &h) const { return h.holderPos == holderPos && h.holderName == holderName; }
    bool operator!=(const QHolder &h) const { return h.holderPos != holderPos || h.holderName != holderName; }
    QString holderName;
    int holderPos;
};

class Q_SQL_EXPORT QSqlResultPrivate
{
    Q_DECLARE_PUBLIC(QSqlResult)

public:
    QSqlResultPrivate(QSqlResult *q, const QSqlDriver *drv)
      : q_ptr(q),
        sqldriver(const_cast<QSqlDriver*>(drv))
    { }
    virtual ~QSqlResultPrivate() { }

    void clearValues()
    {
        values.clear();
        bindCount = 0;
    }

    void resetBindCount()
    {
        bindCount = 0;
    }

    void clearIndex()
    {
        indexes.clear();
        holders.clear();
        types.clear();
    }

    void clear()
    {
        clearValues();
        clearIndex();;
    }

    virtual QString fieldSerial(int) const;
    QString positionalToNamedBinding(const QString &query) const;
    QString namedToPositionalBinding(const QString &query);
    QString holderAt(int index) const;

    QSqlResult *q_ptr = nullptr;
    QPointer<QSqlDriver> sqldriver;
    QString sql;
    QSqlError error;
    QSql::NumericalPrecisionPolicy precisionPolicy = QSql::LowPrecisionDouble;
    QSqlResult::BindingSyntax binds = QSqlResult::PositionalBinding;
    int idx = QSql::BeforeFirstRow;
    int bindCount = 0;
    bool active = false;
    bool isSel = false;
    bool forwardOnly = false;

    QString executedQuery;
    QHash<int, QSql::ParamType> types;
    QVector<QVariant> values;
    typedef QHash<QString, QVector<int> > IndexMap;
    IndexMap indexes;

    typedef QVector<QHolder> QHolderVector;
    QHolderVector holders;
};

QT_END_NAMESPACE

#endif // QSQLRESULT_P_H
