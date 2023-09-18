/****************************************************************************
**
** Copyright (C) 2018 The Qt Company Ltd.
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

#ifndef QMAPCIRCLEOBJECT_P_P_H
#define QMAPCIRCLEOBJECT_P_P_H

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
#include <QGeoCoordinate>
#include <QColor>

QT_BEGIN_NAMESPACE

class Q_LOCATION_PRIVATE_EXPORT QMapCircleObjectPrivate : public QGeoMapObjectPrivate
{
public:
    QMapCircleObjectPrivate(QGeoMapObject *q);
    ~QMapCircleObjectPrivate() override;

    virtual QGeoMapObject::Type type() const override final;

    virtual QGeoCoordinate center() const = 0;
    virtual void setCenter(const QGeoCoordinate &center) = 0;
    virtual qreal radius() const = 0;
    virtual void setRadius(qreal radius) = 0;
    virtual QColor color() const = 0;
    virtual void setColor(const QColor &color) = 0;
    virtual QColor borderColor() const = 0;
    virtual void setBorderColor(const QColor &color) = 0;
    virtual qreal borderWidth() const = 0;
    virtual void setBorderWidth(qreal width) = 0;

    // QGeoMapObjectPrivate interface
    bool equals(const QGeoMapObjectPrivate &other) const override;
    virtual QGeoShape geoShape() const override;
    virtual void setGeoShape(const QGeoShape &shape) override;
};


class Q_LOCATION_PRIVATE_EXPORT QMapCircleObjectPrivateDefault : public QMapCircleObjectPrivate
{
public:
    QMapCircleObjectPrivateDefault(QGeoMapObject *q);
    QMapCircleObjectPrivateDefault(const QMapCircleObjectPrivate &other);
    ~QMapCircleObjectPrivateDefault() override;

    // QMapCircleObjectPrivate interface
    QGeoCoordinate center() const override;
    void setCenter(const QGeoCoordinate &center) override;
    qreal radius() const override;
    void setRadius(qreal radius) override;
    QColor color() const override;
    void setColor(const QColor &color) override;
    QColor borderColor() const override;
    void setBorderColor(const QColor &color) override;
    qreal borderWidth() const override;
    void setBorderWidth(qreal width) override;

    // QGeoMapObjectPrivate interface
    QGeoMapObjectPrivate *clone() override;

public:
    QGeoCoordinate m_center;
    qreal m_radius = 0;
    QColor m_fillColor = Qt::transparent;
    QColor m_borderColor;
    qreal m_borderWidth = 1.0;

private:
    QMapCircleObjectPrivateDefault(const QMapCircleObjectPrivateDefault &other) = delete;
};

QT_END_NAMESPACE

#endif // QMAPCIRCLEOBJECT_P_P_H
