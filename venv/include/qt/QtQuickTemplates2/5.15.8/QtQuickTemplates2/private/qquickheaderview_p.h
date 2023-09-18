/****************************************************************************
**
** Copyright (C) 2020 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtQuick module of the Qt Toolkit.
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

#ifndef QQUICKHEADERVIEW_P_H
#define QQUICKHEADERVIEW_P_H

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

#include <private/qquicktableview_p.h>
#include <private/qtquicktemplates2global_p.h>

QT_BEGIN_NAMESPACE

class QQuickHeaderViewBase;
class QQuickHeaderViewBasePrivate;
class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickHeaderViewBase : public QQuickTableView
{
    Q_OBJECT
    Q_DECLARE_PRIVATE(QQuickHeaderViewBase)
    Q_PROPERTY(QString textRole READ textRole WRITE setTextRole NOTIFY textRoleChanged FINAL)

public:
    explicit QQuickHeaderViewBase(Qt::Orientation orient, QQuickItem *parent = nullptr);
    ~QQuickHeaderViewBase();

    QString textRole() const;
    void setTextRole(const QString &role);

protected:
    QQuickHeaderViewBase(QQuickHeaderViewBasePrivate &dd, QQuickItem *parent);

Q_SIGNALS:
    void textRoleChanged();

private:
    Q_DISABLE_COPY(QQuickHeaderViewBase)
    friend class QQuickHorizontalHeaderView;
    friend class QQuickVerticalHeaderView;
};

class QQuickHorizontalHeaderViewPrivate;
class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickHorizontalHeaderView : public QQuickHeaderViewBase
{
    Q_OBJECT
    Q_DECLARE_PRIVATE(QQuickHorizontalHeaderView)

public:
    QQuickHorizontalHeaderView(QQuickItem *parent = nullptr);
    ~QQuickHorizontalHeaderView() override;

protected:
    QQuickHorizontalHeaderView(QQuickHorizontalHeaderViewPrivate &dd, QQuickItem *parent);

private:
    Q_DISABLE_COPY(QQuickHorizontalHeaderView)
};

class QQuickVerticalHeaderViewPrivate;
class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickVerticalHeaderView : public QQuickHeaderViewBase
{
    Q_OBJECT
    Q_DECLARE_PRIVATE(QQuickVerticalHeaderView)

public:
    QQuickVerticalHeaderView(QQuickItem *parent = nullptr);
    ~QQuickVerticalHeaderView() override;

protected:
    QQuickVerticalHeaderView(QQuickVerticalHeaderViewPrivate &dd, QQuickItem *parent);

private:
    Q_DISABLE_COPY(QQuickVerticalHeaderView)
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QQuickHorizontalHeaderView)
QML_DECLARE_TYPE(QQuickVerticalHeaderView)

#endif // QQUICKHEADERVIEW_P_H
