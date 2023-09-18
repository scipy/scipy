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

#ifndef QSQLCACHEDRESULT_P_H
#define QSQLCACHEDRESULT_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of other Qt classes.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <QtSql/private/qtsqlglobal_p.h>
#include "QtSql/qsqlresult.h"
#include "QtSql/private/qsqlresult_p.h"

QT_BEGIN_NAMESPACE

class QVariant;
template <typename T> class QVector;

class QSqlCachedResultPrivate;

class Q_SQL_EXPORT QSqlCachedResult: public QSqlResult
{
    Q_DECLARE_PRIVATE(QSqlCachedResult)

public:
    typedef QVector<QVariant> ValueCache;

protected:
    QSqlCachedResult(QSqlCachedResultPrivate &d);

    void init(int colCount);
    void cleanup();
    void clearValues();

    virtual bool gotoNext(ValueCache &values, int index) = 0;

    QVariant data(int i) override;
    bool isNull(int i) override;
    bool fetch(int i) override;
    bool fetchNext() override;
    bool fetchPrevious() override;
    bool fetchFirst() override;
    bool fetchLast() override;

    int colCount() const;
    ValueCache &cache();

    void virtual_hook(int id, void *data) override;
    void detachFromResultSet() override;
    void setNumericalPrecisionPolicy(QSql::NumericalPrecisionPolicy policy) override;
private:
    bool cacheNext();
};

class Q_SQL_EXPORT QSqlCachedResultPrivate: public QSqlResultPrivate
{
    Q_DECLARE_PUBLIC(QSqlCachedResult)

public:
    QSqlCachedResultPrivate(QSqlCachedResult *q, const QSqlDriver *drv);
    bool canSeek(int i) const;
    inline int cacheCount() const;
    void init(int count, bool fo);
    void cleanup();
    int nextIndex();
    void revertLast();

    QSqlCachedResult::ValueCache cache;
    int rowCacheEnd;
    int colCount;
    bool atEnd;
};

QT_END_NAMESPACE

#endif // QSQLCACHEDRESULT_P_H
