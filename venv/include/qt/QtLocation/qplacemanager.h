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

#ifndef QPLACEMANAGER_H
#define QPLACEMANAGER_H

#include <QtLocation/QPlaceContentReply>
#include <QtLocation/QPlaceContentRequest>
#include <QtLocation/QPlaceIdReply>
#include <QtLocation/QPlaceReply>
#include <QtLocation/QPlaceDetailsReply>
#include <QtLocation/QPlaceMatchReply>
#include <QtLocation/QPlaceMatchRequest>
#include <QtLocation/QPlaceSearchSuggestionReply>
#include <QtLocation/QPlaceSearchRequest>
#include <QtLocation/QPlaceSearchResult>

#include <QtCore/QLocale>
#include <QtCore/QVector>
#include <QtCore/QString>
#include <QtCore/QObject>
#include <QtLocation/QPlaceIcon>

QT_BEGIN_NAMESPACE

class QPlaceManagerEngine;
class QPlaceSearchRequest;
class QPlaceSearchReply;

class Q_LOCATION_EXPORT QPlaceManager : public QObject
{
    Q_OBJECT
public:
    ~QPlaceManager();

    QString managerName() const;
    int managerVersion() const;

    QPlaceDetailsReply *getPlaceDetails(const QString &placeId) const;

    QPlaceContentReply *getPlaceContent(const QPlaceContentRequest &request) const;

    QPlaceSearchReply *search(const QPlaceSearchRequest &query) const;

    QPlaceSearchSuggestionReply *searchSuggestions(const QPlaceSearchRequest &request) const;

    QPlaceIdReply *savePlace(const QPlace &place);
    QPlaceIdReply *removePlace(const QString &placeId);

    QPlaceIdReply *saveCategory(const QPlaceCategory &category, const QString &parentId = QString());
    QPlaceIdReply *removeCategory(const QString &categoryId);

    QPlaceReply *initializeCategories();
    QString parentCategoryId(const QString &categoryId) const;
    QStringList childCategoryIds(const QString &parentId = QString()) const;

    QPlaceCategory category(const QString &categoryId) const;
    QList<QPlaceCategory> childCategories(const QString &parentId = QString()) const;

    QList<QLocale> locales() const;
    void setLocale(const QLocale &locale);
    void setLocales(const QList<QLocale> &locale);

    QPlace compatiblePlace(const QPlace &place);

    QPlaceMatchReply *matchingPlaces(const QPlaceMatchRequest &request) const;

Q_SIGNALS:
    void finished(QPlaceReply *reply);
    void error(QPlaceReply *, QPlaceReply::Error error, const QString &errorString = QString());

    void placeAdded(const QString &placeId);
    void placeUpdated(const QString &placeId);
    void placeRemoved(const QString &placeId);

    void categoryAdded(const QPlaceCategory &category, const QString &parentId);
    void categoryUpdated(const QPlaceCategory &category, const QString &parentId);
    void categoryRemoved(const QString &categoryId, const QString &parentId);
    void dataChanged();

private:
    explicit QPlaceManager(QPlaceManagerEngine *engine, QObject *parent = nullptr);
    Q_DISABLE_COPY(QPlaceManager)

    QPlaceManagerEngine *d;

    friend class QGeoServiceProvider;
    friend class QGeoServiceProviderPrivate;
    friend class QPlaceIcon;
};

QT_END_NAMESPACE

#endif // QPLACEMANAGER_H
