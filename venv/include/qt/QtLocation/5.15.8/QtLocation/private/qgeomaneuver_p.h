/****************************************************************************
**
** Copyright (C) 2015 The Qt Company Ltd.
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

#ifndef QGEOMANEUVER_P_H
#define QGEOMANEUVER_P_H

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
#include <QtLocation/qgeomaneuver.h>
#include <QtPositioning/qgeocoordinate.h>

#include <QSharedData>
#include <QString>

QT_BEGIN_NAMESPACE

class Q_LOCATION_PRIVATE_EXPORT QGeoManeuverPrivate : public QSharedData
{
public:
    QGeoManeuverPrivate();
    QGeoManeuverPrivate(const QGeoManeuverPrivate &other);
    virtual ~QGeoManeuverPrivate();
    virtual QGeoManeuverPrivate *clone() = 0;

    bool operator== (const QGeoManeuverPrivate &other) const;

    virtual bool valid() const;
    virtual void setValid(bool valid);

    virtual QString id() const;
    virtual void setId(const QString id);

    virtual QGeoCoordinate position() const;
    virtual void setPosition(const QGeoCoordinate &position);

    virtual QString text() const;
    virtual void setText(const QString &text);

    virtual QGeoManeuver::InstructionDirection direction() const;
    virtual void setDirection(QGeoManeuver::InstructionDirection direction);

    virtual int timeToNextInstruction() const;
    virtual void setTimeToNextInstruction(int timeToNextInstruction);

    virtual qreal distanceToNextInstruction() const;
    virtual void setDistanceToNextInstruction(qreal distanceToNextInstruction);

    virtual QGeoCoordinate waypoint() const;
    virtual void setWaypoint(const QGeoCoordinate &waypoint);

    virtual QVariantMap extendedAttributes() const;
    virtual void setExtendedAttributes(const QVariantMap &extendedAttributes);

protected:
    virtual bool equals(const QGeoManeuverPrivate &other) const;
};

class Q_LOCATION_PRIVATE_EXPORT QGeoManeuverPrivateDefault : public QGeoManeuverPrivate
{
public:
    QGeoManeuverPrivateDefault();
    QGeoManeuverPrivateDefault(const QGeoManeuverPrivateDefault &other);
    ~QGeoManeuverPrivateDefault();
    virtual QGeoManeuverPrivate *clone() override;

    virtual bool valid() const override;
    virtual void setValid(bool valid) override;

    virtual QString id() const override;
    virtual void setId(const QString id) override;

    virtual QGeoCoordinate position() const override;
    virtual void setPosition(const QGeoCoordinate &position) override;

    virtual QString text() const override;
    virtual void setText(const QString &text) override;

    virtual QGeoManeuver::InstructionDirection direction() const override;
    virtual void setDirection(QGeoManeuver::InstructionDirection direction) override;

    virtual int timeToNextInstruction() const override;
    virtual void setTimeToNextInstruction(int timeToNextInstruction) override;

    virtual qreal distanceToNextInstruction() const override;
    virtual void setDistanceToNextInstruction(qreal distanceToNextInstruction) override;

    virtual QGeoCoordinate waypoint() const override;
    virtual void setWaypoint(const QGeoCoordinate &waypoint) override;

    virtual QVariantMap extendedAttributes() const override;
    virtual void setExtendedAttributes(const QVariantMap &extendedAttributes) override;

    bool m_valid;
    QString m_id;
    QGeoCoordinate m_position;
    QString m_text;
    QGeoManeuver::InstructionDirection m_direction;
    int m_timeToNextInstruction;
    qreal m_distanceToNextInstruction;
    QGeoCoordinate m_waypoint;
    QVariantMap m_extendedAttributes;
};

QT_END_NAMESPACE

#endif
