/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtPositioning module of the Qt Toolkit.
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

#ifndef QGEOAREAMONITORINFO_H
#define QGEOAREAMONITORINFO_H

#include <QtPositioning/QGeoCoordinate>
#include <QtPositioning/QGeoShape>
#include <QtCore/QSharedDataPointer>
#include <QtCore/QVariantMap>

QT_BEGIN_NAMESPACE

class QDataStream;
class QGeoAreaMonitorInfo;

#ifndef QT_NO_DATASTREAM
Q_POSITIONING_EXPORT QDataStream &operator<<(QDataStream &, const QGeoAreaMonitorInfo &);
Q_POSITIONING_EXPORT QDataStream &operator>>(QDataStream &, QGeoAreaMonitorInfo &);
#endif

class QGeoAreaMonitorInfoPrivate;
class Q_POSITIONING_EXPORT QGeoAreaMonitorInfo
{
public:
    explicit QGeoAreaMonitorInfo(const QString &name = QString());
    QGeoAreaMonitorInfo(const QGeoAreaMonitorInfo &other);
    ~QGeoAreaMonitorInfo();

    QGeoAreaMonitorInfo &operator=(const QGeoAreaMonitorInfo &other);

    bool operator==(const QGeoAreaMonitorInfo &other) const;
    bool operator!=(const QGeoAreaMonitorInfo &other) const;

    QString name() const;
    void setName(const QString &name);

    QString identifier() const;
    bool isValid() const;

    QGeoShape area() const;
    void setArea(const QGeoShape &newShape);

    QDateTime expiration() const;
    void setExpiration(const QDateTime &expiry);

    bool isPersistent() const;
    void setPersistent(bool isPersistent);

    QVariantMap notificationParameters() const;
    void setNotificationParameters(const QVariantMap &parameters);
private:
    QSharedDataPointer<QGeoAreaMonitorInfoPrivate> d;

#ifndef QT_NO_DATASTREAM
    friend Q_POSITIONING_EXPORT QDataStream &operator<<(QDataStream &, const QGeoAreaMonitorInfo &);
    friend Q_POSITIONING_EXPORT QDataStream &operator>>(QDataStream &, QGeoAreaMonitorInfo &);
#endif
};

#ifndef QT_NO_DEBUG_STREAM
Q_POSITIONING_EXPORT QDebug operator<<(QDebug, const QGeoAreaMonitorInfo &);
#endif

QT_END_NAMESPACE

#endif // QGEOAREAMONITORINFO_H
