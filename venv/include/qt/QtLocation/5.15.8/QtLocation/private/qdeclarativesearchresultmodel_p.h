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

#ifndef QDECLARATIVESEARCHRESULTMODEL_P_H
#define QDECLARATIVESEARCHRESULTMODEL_P_H

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
#include <QtLocation/private/qdeclarativesearchmodelbase_p.h>
#include <QtLocation/private/qdeclarativecategory_p.h>
#include <QtLocation/private/qdeclarativeplace_p.h>
#include <QtLocation/private/qdeclarativeplaceicon_p.h>

QT_BEGIN_NAMESPACE

class QDeclarativeGeoServiceProvider;

class Q_LOCATION_PRIVATE_EXPORT QDeclarativeSearchResultModel : public QDeclarativeSearchModelBase
{
    Q_OBJECT

    Q_PROPERTY(QString searchTerm READ searchTerm WRITE setSearchTerm NOTIFY searchTermChanged)
    Q_PROPERTY(QQmlListProperty<QDeclarativeCategory> categories READ categories NOTIFY categoriesChanged)
    Q_PROPERTY(QString recommendationId READ recommendationId WRITE setRecommendationId NOTIFY recommendationIdChanged)
    Q_PROPERTY(RelevanceHint relevanceHint READ relevanceHint WRITE setRelevanceHint NOTIFY relevanceHintChanged)
    Q_PROPERTY(QDeclarativePlace::Visibility visibilityScope READ visibilityScope WRITE setVisibilityScope NOTIFY visibilityScopeChanged)

    Q_PROPERTY(int count READ rowCount NOTIFY rowCountChanged)
    Q_PROPERTY(QDeclarativeGeoServiceProvider *favoritesPlugin READ favoritesPlugin WRITE setFavoritesPlugin NOTIFY favoritesPluginChanged)
    Q_PROPERTY(QVariantMap favoritesMatchParameters READ favoritesMatchParameters WRITE setFavoritesMatchParameters NOTIFY favoritesMatchParametersChanged)

    Q_PROPERTY(bool incremental MEMBER m_incremental NOTIFY incrementalChanged REVISION 12)

    Q_ENUMS(SearchResultType RelevanceHint)

public:
    enum SearchResultType {
        UnknownSearchResult = QPlaceSearchResult::UnknownSearchResult,
        PlaceResult = QPlaceSearchResult::PlaceResult,
        ProposedSearchResult = QPlaceSearchResult::ProposedSearchResult
    };

    enum RelevanceHint {
        UnspecifiedHint = QPlaceSearchRequest::UnspecifiedHint,
        DistanceHint = QPlaceSearchRequest::DistanceHint,
        LexicalPlaceNameHint = QPlaceSearchRequest::LexicalPlaceNameHint
    };

    explicit QDeclarativeSearchResultModel(QObject *parent = 0);
    ~QDeclarativeSearchResultModel();

    QString searchTerm() const;
    void setSearchTerm(const QString &searchTerm);

    QQmlListProperty<QDeclarativeCategory> categories();
    static void categories_append(QQmlListProperty<QDeclarativeCategory> *list,
                                  QDeclarativeCategory *category);
    static int categories_count(QQmlListProperty<QDeclarativeCategory> *list);
    static QDeclarativeCategory *category_at(QQmlListProperty<QDeclarativeCategory> *list, int index);
    static void categories_clear(QQmlListProperty<QDeclarativeCategory> *list);

    QString recommendationId() const;
    void setRecommendationId(const QString &recommendationId);

    QDeclarativeSearchResultModel::RelevanceHint relevanceHint() const;
    void setRelevanceHint(QDeclarativeSearchResultModel::RelevanceHint hint);

    QDeclarativePlace::Visibility visibilityScope() const;
    void setVisibilityScope(QDeclarativePlace::Visibility visibilityScope);

    QDeclarativeGeoServiceProvider *favoritesPlugin() const;
    void setFavoritesPlugin(QDeclarativeGeoServiceProvider *plugin);

    QVariantMap favoritesMatchParameters() const;
    void setFavoritesMatchParameters(const QVariantMap &parameters);

    int rowCount(const QModelIndex &parent = QModelIndex()) const override;

    virtual void clearData(bool suppressSignal = false) override;
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const override;
    Q_INVOKABLE QVariant data(int index, const QString &roleName) const;
    QHash<int, QByteArray> roleNames() const override;

    Q_INVOKABLE void updateWith(int proposedSearchIndex);

    void updateSearchRequest();

Q_SIGNALS:
    void searchTermChanged();
    void categoriesChanged();
    void recommendationIdChanged();
    void relevanceHintChanged();
    void visibilityScopeChanged();

    void rowCountChanged();
    void favoritesPluginChanged();
    void favoritesMatchParametersChanged();
    void dataChanged();
    void incrementalChanged();

protected:
    QPlaceReply *sendQuery(QPlaceManager *manager, const QPlaceSearchRequest &request) override;
    virtual void initializePlugin(QDeclarativeGeoServiceProvider *plugin) override;

protected Q_SLOTS:
    virtual void queryFinished() override;
    virtual void onContentUpdated() override;

private Q_SLOTS:
    void updateLayout(const QList<QPlace> &favoritePlaces = QList<QPlace>());

    void placeUpdated(const QString &placeId);
    void placeRemoved(const QString &placeId);

private:
    enum Roles {
        SearchResultTypeRole = Qt::UserRole,
        TitleRole,
        IconRole,
        DistanceRole,
        PlaceRole,
        SponsoredRole
    };

    int getRow(const QString &placeId) const;
    QList<QPlaceSearchResult> resultsFromPages() const;
    void removePageRow(int row);

    QList<QDeclarativeCategory *> m_categories;
    QLocation::VisibilityScope m_visibilityScope;

    QMap<int, QList<QPlaceSearchResult>> m_pages;
    QList<QPlaceSearchResult> m_results;
    QList<QPlaceSearchResult> m_resultsBuffer;
    QList<QDeclarativePlace *> m_places;
    QList<QDeclarativePlaceIcon *> m_icons;

    QDeclarativeGeoServiceProvider *m_favoritesPlugin;
    QVariantMap m_matchParameters;
    bool m_incremental = false;
};

QT_END_NAMESPACE

#endif // QDECLARATIVESEARCHRESULTMODEL_P_H
