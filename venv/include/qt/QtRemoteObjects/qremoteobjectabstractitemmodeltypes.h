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

#ifndef QREMOTEOBJECTS_ABSTRACT_ITEM_MODEL_TYPES_H
#define QREMOTEOBJECTS_ABSTRACT_ITEM_MODEL_TYPES_H

#include <QtCore/qdatastream.h>
#include <QtCore/qlist.h>
#include <QtCore/qvector.h>
#include <QtCore/qpair.h>
#include <QtCore/qvariant.h>
#include <QtCore/qabstractitemmodel.h>
#include <QtCore/qitemselectionmodel.h>
#include <QtCore/qsize.h>
#include <QtCore/qdebug.h>
#include <QtCore/qnamespace.h>
#include <QtRemoteObjects/qtremoteobjectglobal.h>

QT_BEGIN_NAMESPACE

struct ModelIndex
{
    ModelIndex() : row(-1), column(-1) {}
    ModelIndex(int row_, int column_)
        : row(row_)
        , column(column_)
    {}

    inline bool operator==(const ModelIndex &other) const { return row == other.row && column == other.column; }
    inline bool operator!=(const ModelIndex &other) const { return !(*this == other); }
    int row;
    int column;
};

typedef QList<ModelIndex> IndexList;

struct IndexValuePair
{
    explicit IndexValuePair(const IndexList index_ = IndexList(), const QVariantList &data_ = QVariantList(),
                            bool hasChildren_ = false, const Qt::ItemFlags &flags_ = Qt::ItemFlags(), const QSize &size_ = {})
        : index(index_)
        , data(data_)
        , flags(flags_)
        , hasChildren(hasChildren_)
        , size(size_)
    {}

    inline bool operator==(const IndexValuePair &other) const { return index == other.index && data == other.data && hasChildren == other.hasChildren && flags == other.flags; }
    inline bool operator!=(const IndexValuePair &other) const { return !(*this == other); }

    IndexList index;
    QVariantList data;
    Qt::ItemFlags flags;
    bool hasChildren;
    QVector<IndexValuePair> children;
    QSize size;
};

struct DataEntries
{
    inline bool operator==(const DataEntries &other) const { return data == other.data; }
    inline bool operator!=(const DataEntries &other) const { return !(*this == other); }

    QVector<IndexValuePair> data;
};

struct MetaAndDataEntries : DataEntries
{
    QVector<int> roles;
    QSize size;
};

inline QDebug operator<<(QDebug stream, const ModelIndex &index)
{
    return stream.nospace() << "ModelIndex[row=" << index.row << ", column=" << index.column << "]";
}

inline QDataStream& operator<<(QDataStream &stream, const ModelIndex &index)
{
    return stream << index.row << index.column;
}

inline QDataStream& operator>>(QDataStream &stream, ModelIndex &index)
{
    return stream >> index.row >> index.column;
}

inline QDataStream& operator<<(QDataStream &stream, Qt::Orientation orient)
{
    return stream << static_cast<int>(orient);
}

inline QDataStream& operator>>(QDataStream &stream, Qt::Orientation &orient)
{
    int val;
    QDataStream &ret = stream >> val;
    orient = static_cast<Qt::Orientation>(val);
    return ret;
}

inline QDataStream& operator<<(QDataStream &stream, QItemSelectionModel::SelectionFlags command)
{
    return stream << static_cast<int>(command);
}

inline QDataStream& operator>>(QDataStream &stream, QItemSelectionModel::SelectionFlags &command)
{
    int val;
    QDataStream &ret = stream >> val;
    command = static_cast<QItemSelectionModel::SelectionFlags>(val);
    return ret;
}

inline QDebug operator<<(QDebug stream, const IndexValuePair &pair)
{
    return stream.nospace() << "IndexValuePair[index=" << pair.index << ", data=" << pair.data << ", hasChildren=" << pair.hasChildren << ", flags=" << pair.flags << "]";
}

inline QDebug operator<<(QDebug stream, const DataEntries &entries)
{
    return stream.nospace() << "DataEntries[" << entries.data << "]";
}

inline QDataStream& operator<<(QDataStream &stream, const DataEntries &entries)
{
    return stream << entries.data;
}

inline QDataStream& operator>>(QDataStream &stream, DataEntries &entries)
{
    return stream >> entries.data;
}

inline QDataStream& operator<<(QDataStream &stream, const MetaAndDataEntries &entries)
{
    return stream << entries.data << entries.roles << entries.size;
}

inline QDataStream& operator>>(QDataStream &stream, MetaAndDataEntries &entries)
{
    return stream >> entries.data >> entries.roles >> entries.size;
}

inline QDataStream& operator<<(QDataStream &stream, const IndexValuePair &pair)
{
    return stream << pair.index << pair.data << pair.hasChildren << static_cast<int>(pair.flags) << pair.children << pair.size;
}

inline QDataStream& operator>>(QDataStream &stream, IndexValuePair &pair)
{
    int flags;
    QDataStream &ret = stream >> pair.index >> pair.data >> pair.hasChildren >> flags >> pair.children >> pair.size;
    pair.flags = static_cast<Qt::ItemFlags>(flags);
    return ret;
}

inline QString modelIndexToString(const IndexList &list)
{
    QString s;
    QDebug(&s) << list;
    return s;
}

inline QString modelIndexToString(const ModelIndex &index)
{
    QString s;
    QDebug(&s) << index;
    return s;
}

inline QModelIndex toQModelIndex(const IndexList &list, const QAbstractItemModel *model, bool *ok = nullptr, bool ensureItem = false)
{
    if (ok)
        *ok = true;
    QModelIndex result;
    for (int i = 0; i < list.count(); ++i) {
        const ModelIndex &index = list[i];
        if (ensureItem)
            const_cast<QAbstractItemModel *>(model)->setData(result, index.row, Qt::UserRole - 1);

        result = model->index(index.row, index.column, result);
        if (!result.isValid()) {
            if (ok) {
                *ok = false;
            } else {
                qFatal("Internal error: invalid index=%s in indexList=%s",
                       qPrintable(modelIndexToString(list[i])), qPrintable(modelIndexToString(list)));
            }
            return QModelIndex();
        }
    }
    return result;
}

inline IndexList toModelIndexList(const QModelIndex &index, const QAbstractItemModel *model)
{
    IndexList list;
    if (index.isValid()) {
        list << ModelIndex(index.row(), index.column());
        for (QModelIndex curIndex = model->parent(index); curIndex.isValid(); curIndex = model->parent(curIndex))
            list.prepend(ModelIndex(curIndex.row(), curIndex.column()));
    }
    return list;
}

QT_END_NAMESPACE

Q_DECLARE_METATYPE(ModelIndex)
Q_DECLARE_METATYPE(IndexList)
Q_DECLARE_METATYPE(DataEntries)
Q_DECLARE_METATYPE(MetaAndDataEntries)
Q_DECLARE_METATYPE(IndexValuePair)
Q_DECLARE_METATYPE(Qt::Orientation)
Q_DECLARE_METATYPE(QItemSelectionModel::SelectionFlags)

#endif // QREMOTEOBJECTS_ABSTRACT_ITEM_MODEL_TYPES_H
