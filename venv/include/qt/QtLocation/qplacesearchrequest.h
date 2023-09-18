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

#ifndef QPLACESEARCHREQUEST_H
#define QPLACESEARCHREQUEST_H

#include <QtCore/QSharedDataPointer>
#include <QtLocation/placemacro.h>
#include <QtLocation/QPlaceCategory>
#include <QtLocation/QLocation>

QT_BEGIN_NAMESPACE

class QGeoShape;
class QPlaceSearchRequestPrivate;

class Q_LOCATION_EXPORT QPlaceSearchRequest
{
public:
    enum RelevanceHint {
        UnspecifiedHint,
        DistanceHint,
        LexicalPlaceNameHint
    };

    QPlaceSearchRequest();
    QPlaceSearchRequest(const QPlaceSearchRequest &other);


    QPlaceSearchRequest &operator=(const QPlaceSearchRequest &other);

    bool operator==(const QPlaceSearchRequest &other) const;
    bool operator!=(const QPlaceSearchRequest &other) const;

    ~QPlaceSearchRequest();

    QString searchTerm() const;
    void setSearchTerm(const QString &term);

    QList<QPlaceCategory> categories() const;
    void setCategory(const QPlaceCategory &category);
    void setCategories(const QList<QPlaceCategory> &categories);

    QGeoShape searchArea() const;
    void setSearchArea(const QGeoShape &area);

    QString recommendationId() const;
    void setRecommendationId(const QString &recommendationId);

    QVariant searchContext() const;
    void setSearchContext(const QVariant &context);

    QLocation::VisibilityScope visibilityScope() const;
    void setVisibilityScope(QLocation::VisibilityScope visibilityScopes);

    RelevanceHint relevanceHint() const;
    void setRelevanceHint(RelevanceHint hint);

    int limit() const;
    void setLimit(int limit);

    void clear();

private:
    QSharedDataPointer<QPlaceSearchRequestPrivate> d_ptr;
    inline QPlaceSearchRequestPrivate *d_func();
    inline const QPlaceSearchRequestPrivate *d_func() const;

    friend class QPlaceSearchRequestPrivate;
};

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QPlaceSearchRequest::RelevanceHint)

#endif // QPLACESEARCHQUERY_H
