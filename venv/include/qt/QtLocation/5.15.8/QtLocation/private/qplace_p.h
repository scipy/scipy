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

#ifndef QPLACE_P_H
#define QPLACE_P_H

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

#include <QSharedData>
#include <QUrl>

#include <QtLocation/private/qlocationglobal_p.h>
#include <QtLocation/qplace.h>
#include <QtPositioning/qgeoaddress.h>
#include <QtPositioning/qgeorectangle.h>
#include <QtPositioning/qgeocoordinate.h>
#include <QtLocation/qplacesupplier.h>
#include <QtLocation/QPlaceIcon>

QT_BEGIN_NAMESPACE

class Q_LOCATION_PRIVATE_EXPORT QPlacePrivate : public QSharedData
{
public:
    QPlacePrivate();
    QPlacePrivate(const QPlacePrivate &other);
    virtual ~QPlacePrivate();
    virtual QPlacePrivate *clone() = 0;

    bool operator==(const QPlacePrivate &other) const;

    virtual bool isEmpty() const;
    virtual QList<QPlaceCategory> categories() const = 0;
    virtual void setCategories(const QList<QPlaceCategory> &categories) = 0;
    virtual QGeoLocation location() const = 0;
    virtual void setLocation(const QGeoLocation &location) = 0;
    virtual QPlaceRatings ratings() const = 0;
    virtual void setRatings(const QPlaceRatings &ratings) = 0;
    virtual QPlaceSupplier supplier() const = 0;
    virtual void setSupplier(const QPlaceSupplier &supplier) = 0;
    virtual QString name() const = 0;
    virtual void setName(const QString &name) = 0;
    virtual QString placeId() const = 0;
    virtual void setPlaceId(const QString &placeIdentifier) = 0;
    virtual QString attribution() const = 0;
    virtual void setAttribution(const QString &attribution) = 0;
    virtual QLocation::Visibility visibility() const = 0;
    virtual void setVisibility(QLocation::Visibility visibility) = 0;
    virtual QPlaceIcon icon() const = 0;
    virtual void setIcon(const QPlaceIcon &icon) = 0;
    virtual bool detailsFetched() const = 0;
    virtual void setDetailsFetched(bool fetched) = 0;

    virtual QMap<QString, QPlaceAttribute> extendedAttributes() const = 0;
    virtual QMap<QString, QPlaceAttribute> &extendedAttributes() = 0;
    virtual QMap<QString, QList<QPlaceContactDetail> > contacts() const = 0;
    virtual QMap<QString, QList<QPlaceContactDetail> > &contacts() = 0;
    virtual QPlaceAttribute extendedAttribute(const QString &attributeType) const;


    // The place content, that has to be manually retrieved from the place manager and manually added to the place.
    // Currently, place content types can be:
    //    ImageType,
    //    ReviewType,
    //    EditorialType,
    //    CustomType = 0x0100
    QMap<QPlaceContent::Type, QPlaceContent::Collection> m_contentCollections;
    QMap<QPlaceContent::Type, int> m_contentCounts;
};


class Q_LOCATION_PRIVATE_EXPORT QPlacePrivateDefault : public QPlacePrivate
{
public:
    QPlacePrivateDefault();
    QPlacePrivateDefault(const QPlacePrivateDefault &other);
    virtual ~QPlacePrivateDefault();
    virtual QPlacePrivate *clone() override;

    virtual QList<QPlaceCategory> categories() const override;
    virtual void setCategories(const QList<QPlaceCategory> &categories) override;
    virtual QGeoLocation location() const override;
    virtual void setLocation(const QGeoLocation &location) override;
    virtual QPlaceRatings ratings() const override;
    virtual void setRatings(const QPlaceRatings &ratings) override;
    virtual QPlaceSupplier supplier() const override;
    virtual void setSupplier(const QPlaceSupplier &supplier) override;
    virtual QString name() const override;
    virtual void setName(const QString &name) override;
    virtual QString placeId() const override;
    virtual void setPlaceId(const QString &placeIdentifier) override;
    virtual QString attribution() const override;
    virtual void setAttribution(const QString &attribution) override;
    virtual QLocation::Visibility visibility() const override;
    virtual void setVisibility(QLocation::Visibility visibility) override;
    virtual QPlaceIcon icon() const override;
    virtual void setIcon(const QPlaceIcon &icon) override;
    virtual bool detailsFetched() const override;
    virtual void setDetailsFetched(bool fetched) override;

    virtual QMap<QString, QPlaceAttribute> extendedAttributes() const override;
    virtual QMap<QString, QPlaceAttribute> &extendedAttributes() override;
    virtual QMap<QString, QList<QPlaceContactDetail> > contacts() const override;
    virtual QMap<QString, QList<QPlaceContactDetail> > &contacts() override;


    // data members

    QList<QPlaceCategory> m_categories;
    QGeoLocation m_location;
    QPlaceRatings m_ratings;
    QPlaceSupplier m_supplier;
    QString m_name;
    QString m_placeId;
    QString m_attribution;

    QMap<QString, QPlaceAttribute> m_extendedAttributes;
    QMap<QString, QList<QPlaceContactDetail> > m_contacts;

    QLocation::Visibility m_visibility;
    QPlaceIcon m_icon;
    bool m_detailsFetched;
};

QT_END_NAMESPACE

#endif

