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

#ifndef QGEOMAPICON_P_P_H
#define QGEOMAPICON_P_P_H

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

QT_BEGIN_NAMESPACE

class Q_LOCATION_PRIVATE_EXPORT QMapIconObjectPrivate : public QGeoMapObjectPrivate
{
public:
    QMapIconObjectPrivate(QGeoMapObject *q);
    ~QMapIconObjectPrivate() override;

    virtual QGeoMapObject::Type type() const override final;

    virtual QGeoCoordinate coordinate() const = 0;
    virtual void setCoordinate(const QGeoCoordinate &coordinate) = 0;
    virtual QVariant content() const = 0;
    virtual void setContent(const QVariant &content) = 0;
    virtual QSizeF iconSize() const = 0;
    virtual void setIconSize(const QSizeF &size) = 0;

    // QGeoMapObjectPrivate interface
    bool equals(const QGeoMapObjectPrivate &other) const override;
};

class Q_LOCATION_PRIVATE_EXPORT QMapIconObjectPrivateDefault : public QMapIconObjectPrivate
{
public:
    QMapIconObjectPrivateDefault(QGeoMapObject *q);
    QMapIconObjectPrivateDefault(const QMapIconObjectPrivate &other);
    ~QMapIconObjectPrivateDefault() override;

    // QGeoMapIconPrivate interface
    QGeoCoordinate coordinate() const override;
    void setCoordinate(const QGeoCoordinate &coordinate) override;
    QVariant content() const override;
    void setContent(const QVariant &content) override;
    virtual QSizeF iconSize() const override;
    virtual void setIconSize(const QSizeF &size) override;

    // QGeoMapObjectPrivate interface
    QGeoMapObjectPrivate *clone() override;
    QGeoShape geoShape() const override;
    void setGeoShape(const QGeoShape &shape) override;

public:
    QVariant m_content;
    QGeoCoordinate m_coordinate;
    QSizeF m_iconSize;
    qreal m_radius = 100.0; // can be set only via setGeoShape, will be used only by

private:
    QMapIconObjectPrivateDefault(const QMapIconObjectPrivateDefault &other) = delete;
};

QT_END_NAMESPACE

#endif // QGEOMAPICON_P_P_H
