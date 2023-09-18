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

#ifndef QSQLTABLEMODEL_P_H
#define QSQLTABLEMODEL_P_H

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
#include "private/qsqlquerymodel_p.h"
#include "QtSql/qsqlindex.h"
#include "QtCore/qmap.h"

QT_REQUIRE_CONFIG(sqlmodel);

QT_BEGIN_NAMESPACE

class Q_AUTOTEST_EXPORT QSqlTableModelPrivate: public QSqlQueryModelPrivate
{
    Q_DECLARE_PUBLIC(QSqlTableModel)

public:
    QSqlTableModelPrivate()
        : sortColumn(-1),
          sortOrder(Qt::AscendingOrder),
          strategy(QSqlTableModel::OnRowChange),
          busyInsertingRows(false)
    {}
    ~QSqlTableModelPrivate();

    void clear();
    virtual void clearCache();
    QSqlRecord record(const QVector<QVariant> &values) const;

    bool exec(const QString &stmt, bool prepStatement,
              const QSqlRecord &rec, const QSqlRecord &whereValues);
    virtual void revertCachedRow(int row);
    virtual int nameToIndex(const QString &name) const;
    QString strippedFieldName(const QString &name) const;
    int insertCount(int maxRow = -1) const;
    void initRecordAndPrimaryIndex();

    QSqlDatabase db;

    int sortColumn;
    Qt::SortOrder sortOrder;

    QSqlTableModel::EditStrategy strategy;
    bool busyInsertingRows;

    QSqlQuery editQuery = { QSqlQuery(nullptr) };
    QSqlIndex primaryIndex;
    QString tableName;
    QString filter;
    QString autoColumn;

    enum Op { None, Insert, Update, Delete };

    class ModifiedRow
    {
    public:
        inline ModifiedRow(Op o = None, const QSqlRecord &r = QSqlRecord())
            : m_op(None), m_db_values(r), m_insert(o == Insert)
        { setOp(o); }
        inline Op op() const { return m_op; }
        inline void setOp(Op o)
        {
            if (o == None)
                m_submitted = true;
            if (o == m_op)
                return;
            m_submitted = (o != Insert && o != Delete);
            m_op = o;
            m_rec = m_db_values;
            setGenerated(m_rec, m_op == Delete);
        }
        inline const QSqlRecord &rec() const { return m_rec; }
        inline QSqlRecord& recRef() { return m_rec; }
        inline void setValue(int c, const QVariant &v)
        {
            m_submitted = false;
            m_rec.setValue(c, v);
            m_rec.setGenerated(c, true);
        }
        inline bool submitted() const { return m_submitted; }
        inline void setSubmitted()
        {
            m_submitted = true;
            setGenerated(m_rec, false);
            if (m_op == Delete) {
                m_rec.clearValues();
            }
            else {
                m_op = Update;
                m_db_values = m_rec;
                setGenerated(m_db_values, true);
            }
        }
        inline void refresh(bool exists, const QSqlRecord& newvals)
        {
            m_submitted = true;
            if (exists) {
                m_op = Update;
                m_db_values = newvals;
                m_rec = newvals;
                setGenerated(m_rec, false);
            } else {
                m_op = Delete;
                m_rec.clear();
                m_db_values.clear();
            }
        }
        inline bool insert() const { return m_insert; }
        inline void revert()
        {
            if (m_submitted)
                return;
            if (m_op == Delete)
                m_op = Update;
            m_rec = m_db_values;
            setGenerated(m_rec, false);
            m_submitted = true;
        }
        inline QSqlRecord primaryValues(const QSqlRecord& pi) const
        {
            if (m_op == None || m_op == Insert)
                return QSqlRecord();

            return m_db_values.keyValues(pi);
        }
    private:
        inline static void setGenerated(QSqlRecord& r, bool g)
        {
            for (int i = r.count() - 1; i >= 0; --i)
                r.setGenerated(i, g);
        }
        Op m_op;
        QSqlRecord m_rec;
        QSqlRecord m_db_values;
        bool m_submitted;
        bool m_insert;
    };

    typedef QMap<int, ModifiedRow> CacheMap;
    CacheMap cache;
};

class QSqlTableModelSql: public QSqlQueryModelSql
{
public:
};

QT_END_NAMESPACE

#endif // QSQLTABLEMODEL_P_H
