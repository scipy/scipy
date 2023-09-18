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
#ifndef QGEOSATELLITEINFOSOURCE_H
#define QGEOSATELLITEINFOSOURCE_H

#include <QtPositioning/QGeoSatelliteInfo>

#include <QtCore/QObject>
#include <QtCore/QList>

QT_BEGIN_NAMESPACE

class QGeoSatelliteInfoSourcePrivate;
class Q_POSITIONING_EXPORT QGeoSatelliteInfoSource : public QObject
{
    Q_OBJECT
    Q_PROPERTY(int updateInterval READ updateInterval WRITE setUpdateInterval)
    Q_PROPERTY(int minimumUpdateInterval READ minimumUpdateInterval)

public:
    enum Error {
        AccessError = 0,
        ClosedError = 1,
        NoError = 2,
        UnknownSourceError = -1
    };
    Q_ENUMS(Error)

    explicit QGeoSatelliteInfoSource(QObject *parent);
    virtual ~QGeoSatelliteInfoSource();

    static QGeoSatelliteInfoSource *createDefaultSource(QObject *parent);
    static QGeoSatelliteInfoSource *createSource(const QString &sourceName, QObject *parent);
    static QGeoSatelliteInfoSource *createDefaultSource(const QVariantMap &parameters, QObject *parent);
    static QGeoSatelliteInfoSource *createSource(const QString &sourceName, const QVariantMap &parameters, QObject *parent);
    static QStringList availableSources();

    QString sourceName() const;

    virtual void setUpdateInterval(int msec);
    int updateInterval() const;
    virtual int minimumUpdateInterval() const = 0;
    virtual Error error() const = 0;

public Q_SLOTS:
    virtual void startUpdates() = 0;
    virtual void stopUpdates() = 0;

    virtual void requestUpdate(int timeout = 0) = 0;

Q_SIGNALS:
    void satellitesInViewUpdated(const QList<QGeoSatelliteInfo> &satellites);
    void satellitesInUseUpdated(const QList<QGeoSatelliteInfo> &satellites);
    void requestTimeout();
    void error(QGeoSatelliteInfoSource::Error);

protected:
    explicit QGeoSatelliteInfoSource(QGeoSatelliteInfoSourcePrivate &dd, QObject *parent);

private:
    Q_DISABLE_COPY(QGeoSatelliteInfoSource)
    QGeoSatelliteInfoSourcePrivate *d;
    friend class QGeoSatelliteInfoSourcePrivate;
};

QT_END_NAMESPACE

#endif
