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
***************************************************************************/

#ifndef QDECLARATIVEGEOADDRESS_P_H
#define QDECLARATIVEGEOADDRESS_P_H

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

#include <QtCore/QObject>
#include <QtPositioning/QGeoAddress>
#include <QtPositioningQuick/private/qpositioningquickglobal_p.h>

QT_BEGIN_NAMESPACE

class Q_POSITIONINGQUICK_PRIVATE_EXPORT QDeclarativeGeoAddress : public QObject
{
    Q_OBJECT

    Q_PROPERTY(QGeoAddress address READ address WRITE setAddress)
    Q_PROPERTY(QString text READ text WRITE setText NOTIFY textChanged)
    Q_PROPERTY(QString country READ country WRITE setCountry NOTIFY countryChanged)
    Q_PROPERTY(QString countryCode READ countryCode WRITE setCountryCode NOTIFY countryCodeChanged)
    Q_PROPERTY(QString state READ state WRITE setState NOTIFY stateChanged)
    Q_PROPERTY(QString county READ county WRITE setCounty NOTIFY countyChanged)
    Q_PROPERTY(QString city READ city WRITE setCity NOTIFY cityChanged)
    Q_PROPERTY(QString district READ district WRITE setDistrict NOTIFY districtChanged)
    Q_PROPERTY(QString street READ street WRITE setStreet NOTIFY streetChanged)
    Q_PROPERTY(QString postalCode READ postalCode WRITE setPostalCode NOTIFY postalCodeChanged)
    Q_PROPERTY(bool isTextGenerated READ isTextGenerated NOTIFY isTextGeneratedChanged)

public:
    explicit QDeclarativeGeoAddress(QObject *parent = 0);
    QDeclarativeGeoAddress(const QGeoAddress &address, QObject *parent = 0);
    QGeoAddress address() const;
    void setAddress(const QGeoAddress &address);

    QString text() const;
    void setText(const QString &address);

    QString country() const;
    void setCountry(const QString &country);
    QString countryCode() const;
    void setCountryCode(const QString &countryCode);
    QString state() const;
    void setState(const QString &state);
    QString county() const;
    void setCounty(const QString &county);
    QString city() const;
    void setCity(const QString &city);
    QString district() const;
    void setDistrict(const QString &district);
    QString street() const;
    void setStreet(const QString &street);
    QString postalCode() const;
    void setPostalCode(const QString &postalCode);
    bool isTextGenerated() const;

Q_SIGNALS:
    void textChanged();
    void countryChanged();
    void countryCodeChanged();
    void stateChanged();
    void countyChanged();
    void cityChanged();
    void districtChanged();
    void streetChanged();
    void postalCodeChanged();
    void isTextGeneratedChanged();

private:
    QGeoAddress m_address;
};

QT_END_NAMESPACE

#endif // QDECLARATIVEGEOADDRESS_P_H
