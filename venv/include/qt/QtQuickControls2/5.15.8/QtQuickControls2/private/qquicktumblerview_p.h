/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the Qt Quick Controls 2 module of the Qt Toolkit.
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

#ifndef QQUICKTUMBLERVIEW_P_H
#define QQUICKTUMBLERVIEW_P_H

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

#include <QQuickItem>
#include <QtQuickControls2/private/qtquickcontrols2global_p.h>

QT_BEGIN_NAMESPACE

class QQuickListView;
class QQuickPath;
class QQuickPathView;

class QQuickTumbler;

class Q_QUICKCONTROLS2_PRIVATE_EXPORT QQuickTumblerView : public QQuickItem
{
    Q_OBJECT
    Q_PROPERTY(QVariant model READ model WRITE setModel NOTIFY modelChanged)
    Q_PROPERTY(QQmlComponent *delegate READ delegate WRITE setDelegate NOTIFY delegateChanged)
    Q_PROPERTY(QQuickPath *path READ path WRITE setPath NOTIFY pathChanged)

public:
    QQuickTumblerView(QQuickItem *parent = nullptr);

    QVariant model() const;
    void setModel(const QVariant &model);

    QQmlComponent *delegate() const;
    void setDelegate(QQmlComponent *delegate);

    QQuickPath *path() const;
    void setPath(QQuickPath *path);

Q_SIGNALS:
    void modelChanged();
    void delegateChanged();
    void pathChanged();

protected:
    void geometryChanged(const QRectF &newGeometry, const QRectF &oldGeometry) override;
    void componentComplete() override;
    void itemChange(ItemChange change, const ItemChangeData &data) override;

private:
    QQuickItem *view();
    void createView();
    void updateView();
    void updateModel();

    void wrapChange();

    QQuickTumbler *m_tumbler = nullptr;
    QVariant m_model;
    QQmlComponent *m_delegate = nullptr;
    QQuickPathView *m_pathView = nullptr;
    QQuickListView *m_listView = nullptr;
    QQuickPath *m_path = nullptr;
};

QT_END_NAMESPACE

#endif // TUMBLERVIEW_H
