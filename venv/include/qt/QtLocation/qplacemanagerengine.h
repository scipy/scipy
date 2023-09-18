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

#ifndef QPLACEMANAGERENGINE_H
#define QPLACEMANAGERENGINE_H

#include <QtLocation/QPlaceManager>

QT_BEGIN_NAMESPACE

class QPlaceManagerEnginePrivate;
class QPlaceMatchReply;
class QPlaceMatchRequest;
class QPlaceSearchReply;
class QPlaceSearchRequest;
class QPlaceSearchSuggestionReply;

class Q_LOCATION_EXPORT QPlaceManagerEngine : public QObject
{
    Q_OBJECT

public:
    explicit QPlaceManagerEngine(const QVariantMap &parameters, QObject *parent = nullptr);
    virtual ~QPlaceManagerEngine();

    QString managerName() const;
    int managerVersion() const;

    virtual QPlaceDetailsReply *getPlaceDetails(const QString &placeId);

    virtual QPlaceContentReply *getPlaceContent(const QPlaceContentRequest &request);

    virtual QPlaceSearchReply *search(const QPlaceSearchRequest &request);

    virtual QPlaceSearchSuggestionReply *searchSuggestions(const QPlaceSearchRequest &request);

    virtual QPlaceIdReply *savePlace(const QPlace &place);
    virtual QPlaceIdReply *removePlace(const QString &placeId);

    virtual QPlaceIdReply *saveCategory(const QPlaceCategory &category, const QString &parentId);
    virtual QPlaceIdReply *removeCategory(const QString &categoryId);

    virtual QPlaceReply *initializeCategories();
    virtual QString parentCategoryId(const QString &categoryId) const;
    virtual QStringList childCategoryIds(const QString &categoryId) const;
    virtual QPlaceCategory category(const QString &categoryId) const;

    virtual QList<QPlaceCategory> childCategories(const QString &parentId) const;

    virtual QList<QLocale> locales() const;
    virtual void setLocales(const QList<QLocale> &locales);

    virtual QUrl constructIconUrl(const QPlaceIcon &icon, const QSize &size) const;

    virtual QPlace compatiblePlace(const QPlace &original) const;

    virtual QPlaceMatchReply *matchingPlaces(const QPlaceMatchRequest &request);

Q_SIGNALS:
    void finished(QPlaceReply *reply);
    void error(QPlaceReply *, QPlaceReply::Error error, const QString &errorString = QString());

    void placeAdded(const QString &placeId);
    void placeUpdated(const QString &placeId);
    void placeRemoved(const QString &placeId);

    void categoryAdded(const QPlaceCategory &category, const QString &parentCategoryId);
    void categoryUpdated(const QPlaceCategory &category, const QString &parentCategoryId);
    void categoryRemoved(const QString &categoryId, const QString &parentCategoryId);
    void dataChanged();

protected:
    QPlaceManager *manager() const;

private:
    void setManagerName(const QString &managerName);
    void setManagerVersion(int managerVersion);

    QPlaceManagerEnginePrivate *d_ptr;
    Q_DISABLE_COPY(QPlaceManagerEngine)

    friend class QGeoServiceProviderPrivate;
    friend class QPlaceManager;
};

QT_END_NAMESPACE

#endif
