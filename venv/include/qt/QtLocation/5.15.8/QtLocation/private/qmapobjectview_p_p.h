/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the QtLocation module of the Qt Toolkit.
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

#ifndef QMAPOBJECTVIEW_P_P_H
#define QMAPOBJECTVIEW_P_P_H


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

#include <QtLocation/private/qlocationglobal_p.h>
#include <QtLocation/private/qgeomapobject_p_p.h>
#include <QPointer>
#include <QVector>
#include <QQmlComponent>

QT_BEGIN_NAMESPACE

class QQmlDelegateModel;
class QGeoMap;
class Q_LOCATION_PRIVATE_EXPORT QMapObjectViewPrivate : public QGeoMapObjectPrivate
{
public:
    QMapObjectViewPrivate(QGeoMapObject *q);
    ~QMapObjectViewPrivate() override;

    virtual QVariant model() const = 0;
    virtual void setModel(const QVariant &model) = 0;
    virtual QQmlComponent *delegate() const = 0;
    virtual void setDelegate(QQmlComponent *delegate) = 0;

    virtual QGeoMapObject::Type type() const override final;
};

class Q_LOCATION_PRIVATE_EXPORT QMapObjectViewPrivateDefault : public QMapObjectViewPrivate
{
public:
    QMapObjectViewPrivateDefault(QGeoMapObject *q);
    QMapObjectViewPrivateDefault(const QMapObjectViewPrivate &other);
    ~QMapObjectViewPrivateDefault() override;

    virtual QVariant model() const override;
    virtual void setModel(const QVariant &model) override;
    virtual QQmlComponent *delegate() const override;
    virtual void setDelegate(QQmlComponent *delegate) override;

    // QGeoMapObjectPrivate interface
public:
    QGeoMapObjectPrivate *clone() override;
    bool equals(const QGeoMapObjectPrivate &other) const override;
    virtual QGeoShape geoShape() const override;
    virtual void setGeoShape(const QGeoShape &shape) override;

    QVariant m_model;
    QQmlComponent *m_delegate = nullptr;
};

QT_END_NAMESPACE

#endif // QMAPOBJECTVIEW_P_P_H
