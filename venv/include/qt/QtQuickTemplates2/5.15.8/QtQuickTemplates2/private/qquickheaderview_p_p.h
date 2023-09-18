/****************************************************************************
**
** Copyright (C) 2020 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the Qt Quick Templates 2 module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL3$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see http://www.qt.io/terms-conditions. For further
** information use the contact form at http://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPLv3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or later as published by the Free
** Software Foundation and appearing in the file LICENSE.GPL included in
** the packaging of this file. Please review the following information to
** ensure the GNU General Public License version 2.0 requirements will be
** met: http://www.gnu.org/licenses/gpl-2.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QQUICKHEADERVIEW_P_P_H
#define QQUICKHEADERVIEW_P_P_H

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

#include <QtCore/QAbstractItemModel>
#include <QtCore/QPointer>
#include <QtCore/QTransposeProxyModel>
#include <QtQuick/private/qquicktableview_p_p.h>
#include <private/qquickheaderview_p.h>

QT_BEGIN_NAMESPACE

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QHeaderDataProxyModel : public QAbstractItemModel
{
    Q_OBJECT
    Q_DISABLE_COPY(QHeaderDataProxyModel)
    Q_PROPERTY(QAbstractItemModel *sourceModel READ sourceModel)
public:
    explicit QHeaderDataProxyModel(QObject *parent = nullptr);
    ~QHeaderDataProxyModel();

    void setSourceModel(QAbstractItemModel *newSourceModel);
    QPointer<QAbstractItemModel> sourceModel() const;
    QModelIndex index(int row, int column, const QModelIndex &parent = QModelIndex()) const override;
    QModelIndex parent(const QModelIndex &child) const override;
    QModelIndex sibling(int row, int column, const QModelIndex &idx) const override;
    int rowCount(const QModelIndex &parent = QModelIndex()) const override;
    int columnCount(const QModelIndex &parent = QModelIndex()) const override;
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const override;
    bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole) override;
    bool hasChildren(const QModelIndex &parent = QModelIndex()) const override;

    inline QVariant variantValue() const;
    inline Qt::Orientation orientation() const;
    inline void setOrientation(Qt::Orientation o);

private:
    inline void connectToModel();
    inline void disconnectFromModel();
    QPointer<QAbstractItemModel> m_model = nullptr;
    Qt::Orientation m_orientation = Qt::Horizontal;
};

class QQuickHeaderViewBasePrivate : public QQuickTableViewPrivate
{
    Q_DECLARE_PUBLIC(QQuickHeaderViewBase)
public:
    QQuickHeaderViewBasePrivate();
    ~QQuickHeaderViewBasePrivate();

    Qt::Orientation orientation() const;
    void setOrientation(Qt::Orientation orientation);
    const QPointer<QQuickItem> delegateItemAt(int row, int col) const;
    QVariant modelImpl() const override;
    void setModelImpl(const QVariant &newModel) override;
    void syncModel() override;
    void syncSyncView() override;

protected:
    QHeaderDataProxyModel m_headerDataProxyModel;
    QTransposeProxyModel m_transposeProxyModel;
    struct SectionSize
    {
        int section;
        qreal previousSize;
    };
    QStack<SectionSize> m_hiddenSectionSizes;
    bool m_modelExplicitlySetByUser = false;
    QString m_textRole;
};

class QQuickHorizontalHeaderViewPrivate : public QQuickHeaderViewBasePrivate
{
    Q_DECLARE_PUBLIC(QQuickHorizontalHeaderView)
public:
    QQuickHorizontalHeaderViewPrivate();
    ~QQuickHorizontalHeaderViewPrivate();
};

class QQuickVerticalHeaderViewPrivate : public QQuickHeaderViewBasePrivate
{
    Q_DECLARE_PUBLIC(QQuickVerticalHeaderView)
public:
    QQuickVerticalHeaderViewPrivate();
    ~QQuickVerticalHeaderViewPrivate();
};

QT_END_NAMESPACE

#endif // QQUICKHEADERVIEW_P_P_H
