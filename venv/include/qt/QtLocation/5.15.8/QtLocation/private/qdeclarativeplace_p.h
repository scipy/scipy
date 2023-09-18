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

#ifndef QDECLARATIVEPLACE_P_H
#define QDECLARATIVEPLACE_P_H

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
#include <QtCore/QObject>
#include <QtQml/QQmlListProperty>
#include <QtQml/QQmlParserStatus>
#include <QtQml/QQmlPropertyMap>
#include <QtLocation/QPlace>

#include <QtPositioningQuick/private/qdeclarativegeolocation_p.h>
#include <QtLocation/private/qdeclarativecategory_p.h>
#include <QtLocation/private/qdeclarativecontactdetail_p.h>
#include <QtLocation/private/qdeclarativesupplier_p.h>
#include <QtLocation/private/qdeclarativeratings_p.h>
#include <QtLocation/private/qdeclarativereviewmodel_p.h>
#include <QtLocation/private/qdeclarativeplaceimagemodel_p.h>
#include <QtLocation/private/qdeclarativeplaceeditorialmodel_p.h>

QT_BEGIN_NAMESPACE

class QPlaceReply;

class QPlaceManager;
class QDeclarativePlaceIcon;

class Q_LOCATION_PRIVATE_EXPORT QDeclarativePlace : public QObject, public QQmlParserStatus
{
    Q_OBJECT

    Q_ENUMS(Status Visibility)

    Q_PROPERTY(QPlace place READ place WRITE setPlace)
    Q_PROPERTY(QDeclarativeGeoServiceProvider *plugin READ plugin WRITE setPlugin NOTIFY pluginChanged)
    Q_PROPERTY(QQmlListProperty<QDeclarativeCategory> categories READ categories NOTIFY categoriesChanged)
    Q_PROPERTY(QDeclarativeGeoLocation *location READ location WRITE setLocation NOTIFY locationChanged)
    Q_PROPERTY(QDeclarativeRatings *ratings READ ratings WRITE setRatings NOTIFY ratingsChanged)
    Q_PROPERTY(QDeclarativeSupplier *supplier READ supplier WRITE setSupplier NOTIFY supplierChanged)
    Q_PROPERTY(QDeclarativePlaceIcon *icon READ icon WRITE setIcon NOTIFY iconChanged)
    Q_PROPERTY(QString name READ name WRITE setName NOTIFY nameChanged)
    Q_PROPERTY(QString placeId READ placeId WRITE setPlaceId NOTIFY placeIdChanged)
    Q_PROPERTY(QString attribution READ attribution WRITE setAttribution NOTIFY attributionChanged)

    Q_PROPERTY(QDeclarativeReviewModel *reviewModel READ reviewModel NOTIFY reviewModelChanged)
    Q_PROPERTY(QDeclarativePlaceImageModel *imageModel READ imageModel NOTIFY imageModelChanged)
    Q_PROPERTY(QDeclarativePlaceEditorialModel *editorialModel READ editorialModel NOTIFY editorialModelChanged)

    Q_PROPERTY(QObject *extendedAttributes READ extendedAttributes NOTIFY extendedAttributesChanged)
    Q_PROPERTY(QObject *contactDetails READ contactDetails NOTIFY contactDetailsChanged)
    Q_PROPERTY(bool detailsFetched READ detailsFetched NOTIFY detailsFetchedChanged)
    Q_PROPERTY(Status status READ status NOTIFY statusChanged)

    Q_PROPERTY(QString primaryPhone READ primaryPhone NOTIFY primaryPhoneChanged)
    Q_PROPERTY(QString primaryFax READ primaryFax NOTIFY primaryFaxChanged)
    Q_PROPERTY(QString primaryEmail READ primaryEmail NOTIFY primaryEmailChanged)
    Q_PROPERTY(QUrl primaryWebsite READ primaryWebsite NOTIFY primaryWebsiteChanged)

    Q_PROPERTY(Visibility visibility READ visibility WRITE setVisibility NOTIFY visibilityChanged)
    Q_PROPERTY(QDeclarativePlace *favorite READ favorite WRITE setFavorite NOTIFY favoriteChanged)

    Q_INTERFACES(QQmlParserStatus)

public:
    explicit QDeclarativePlace(QObject *parent = 0);
    QDeclarativePlace(const QPlace &src, QDeclarativeGeoServiceProvider *plugin, QObject *parent = 0);
    ~QDeclarativePlace();

    enum Status {Ready, Saving, Fetching, Removing, Error};
    enum Visibility {
        UnspecifiedVisibility = QLocation::UnspecifiedVisibility,
        DeviceVisibility = QLocation::DeviceVisibility,
        PrivateVisibility = QLocation::PrivateVisibility,
        PublicVisibility = QLocation::PublicVisibility
    };

    //From QQmlParserStatus
    virtual void classBegin() {}
    virtual void componentComplete();

    void setPlugin(QDeclarativeGeoServiceProvider *plugin);
    QDeclarativeGeoServiceProvider *plugin() const;

    QDeclarativeReviewModel *reviewModel();
    QDeclarativePlaceImageModel *imageModel();
    QDeclarativePlaceEditorialModel *editorialModel();

    QPlace place();
    void setPlace(const QPlace &src);

    QQmlListProperty<QDeclarativeCategory> categories();
    static void category_append(QQmlListProperty<QDeclarativeCategory> *prop,
                                  QDeclarativeCategory *value);
    static int category_count(QQmlListProperty<QDeclarativeCategory> *prop);
    static QDeclarativeCategory *category_at(QQmlListProperty<QDeclarativeCategory> *prop, int index);
    static void category_clear(QQmlListProperty<QDeclarativeCategory> *prop);

    QDeclarativeGeoLocation *location();
    void setLocation(QDeclarativeGeoLocation *location);
    QDeclarativeRatings *ratings();
    void setRatings(QDeclarativeRatings *ratings);
    QDeclarativeSupplier *supplier() const;
    void setSupplier(QDeclarativeSupplier *supplier);
    QDeclarativePlaceIcon *icon() const;
    void setIcon(QDeclarativePlaceIcon *icon);
    QString name() const;
    void setName(const QString &name);
    QString placeId() const;
    void setPlaceId(const QString &placeId);
    QString attribution() const;
    void setAttribution(const QString &attribution);
    bool detailsFetched() const;

    Status status() const;
    void setStatus(Status status, const QString &errorString = QString());

    Q_INVOKABLE void getDetails();
    Q_INVOKABLE void save();
    Q_INVOKABLE void remove();
    Q_INVOKABLE QString errorString() const;

    QString primaryPhone() const;
    QString primaryFax() const;
    QString primaryEmail() const;
    QUrl primaryWebsite() const;

    QQmlPropertyMap *extendedAttributes() const;

    QDeclarativeContactDetails *contactDetails() const;

    Visibility visibility() const;
    void setVisibility(Visibility visibility);

    QDeclarativePlace *favorite() const;
    void setFavorite(QDeclarativePlace *favorite);

    Q_INVOKABLE void copyFrom(QDeclarativePlace *original);
    Q_INVOKABLE void initializeFavorite(QDeclarativeGeoServiceProvider *plugin);

Q_SIGNALS:
    void pluginChanged();
    void categoriesChanged();
    void locationChanged();
    void ratingsChanged();
    void supplierChanged();
    void iconChanged();
    void nameChanged();
    void placeIdChanged();
    void attributionChanged();
    void detailsFetchedChanged();
    void reviewModelChanged();
    void imageModelChanged();
    void editorialModelChanged();

    void primaryPhoneChanged();
    void primaryFaxChanged();
    void primaryEmailChanged();
    void primaryWebsiteChanged();

    void extendedAttributesChanged();
    void contactDetailsChanged();
    void statusChanged();
    void visibilityChanged();
    void favoriteChanged();

private Q_SLOTS:
    void finished();
    void contactsModified(const QString &, const QVariant &);
    void pluginReady();
    void cleanupDeletedCategories();
private:
    void synchronizeCategories();
    void pullExtendedAttributes();
    void synchronizeContacts();
    void primarySignalsEmission(const QString &type = QString());
    QString primaryValue(const QString &contactType) const;

private:
    QPlaceManager *manager();

    QList<QDeclarativeCategory *> m_categories;
    QDeclarativeGeoLocation *m_location;
    QDeclarativeRatings *m_ratings;
    QDeclarativeSupplier *m_supplier;
    QDeclarativePlaceIcon *m_icon;
    QDeclarativeReviewModel *m_reviewModel;
    QDeclarativePlaceImageModel *m_imageModel;
    QDeclarativePlaceEditorialModel *m_editorialModel;
    QQmlPropertyMap *m_extendedAttributes;
    QDeclarativeContactDetails *m_contactDetails;

    QPlace m_src;

    QPlaceReply *m_reply;

    QDeclarativeGeoServiceProvider *m_plugin;
    bool m_complete;

    QString m_prevPrimaryPhone;
    QString m_prevPrimaryEmail;
    QString m_prevPrimaryFax;
    QUrl m_prevPrimaryWebsite;

    QDeclarativePlace *m_favorite;

    Status m_status;
    QString m_errorString;

    QList<QDeclarativeCategory *>m_categoriesToBeDeleted;
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QDeclarativePlace)

#endif // QDECLARATIVEPLACE_P_H
