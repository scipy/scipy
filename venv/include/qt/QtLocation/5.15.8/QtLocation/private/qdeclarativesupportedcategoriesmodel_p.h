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

#ifndef QDECLARATIVESUPPORTEDCATEGORIESMODEL_H
#define QDECLARATIVESUPPORTEDCATEGORIESMODEL_H

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
#include <QtLocation/private/qdeclarativegeoserviceprovider_p.h>

#include <QObject>
#include <QtCore/QStringList>
#include <QtCore/QSharedPointer>
#include <QAbstractListModel>
#include <QQmlListProperty>
#include <QtQml/QQmlParserStatus>

#include <QtLocation/QPlaceCategory>

#include <QtLocation/private/qdeclarativecategory_p.h>

QT_BEGIN_NAMESPACE

class QGeoServiceProvider;
class QPlaceManager;
class QPlaceReply;

class PlaceCategoryNode
{
public:
    QString parentId;
    QStringList childIds;
    QSharedPointer<QDeclarativeCategory> declCategory;
};

class Q_LOCATION_PRIVATE_EXPORT QDeclarativeSupportedCategoriesModel : public QAbstractItemModel, public QQmlParserStatus
{
    Q_OBJECT

    Q_ENUMS(Status)

    Q_PROPERTY(QDeclarativeGeoServiceProvider *plugin READ plugin WRITE setPlugin NOTIFY pluginChanged)
    Q_PROPERTY(bool hierarchical READ hierarchical WRITE setHierarchical NOTIFY hierarchicalChanged)
    Q_PROPERTY(Status status READ status NOTIFY statusChanged)

    Q_INTERFACES(QQmlParserStatus)
    Q_ENUMS(Roles) //The Roles enum is for internal usage only.

public:
    explicit QDeclarativeSupportedCategoriesModel(QObject *parent = 0);
    virtual ~QDeclarativeSupportedCategoriesModel();

    // From QQmlParserStatus
    virtual void classBegin() {}
    virtual void componentComplete();

    // From QAbstractItemModel
    int rowCount(const QModelIndex &parent) const;
    int columnCount(const QModelIndex &parent) const;

    QModelIndex index(int row, int column, const QModelIndex &parent) const;
    QModelIndex parent(const QModelIndex &child) const;

    Q_INVOKABLE QVariant data(const QModelIndex &index, int role) const;
    QHash<int, QByteArray> roleNames() const;

    enum Roles {
        CategoryRole = Qt::UserRole,
        ParentCategoryRole
    };  //for internal usage only

    enum Status {Null, Ready, Loading, Error};

    void setPlugin(QDeclarativeGeoServiceProvider *plugin);
    QDeclarativeGeoServiceProvider *plugin() const;

    void setHierarchical(bool hierarchical);
    bool hierarchical() const;

    Q_INVOKABLE QString errorString() const;

    Status status() const;
    void setStatus(Status status, const QString &errorString = QString());

    using QAbstractItemModel::dataChanged;
Q_SIGNALS:
    void pluginChanged();
    void hierarchicalChanged();
    void statusChanged();
    void dataChanged();

public Q_SLOTS:
    void update();

private Q_SLOTS:
    void replyFinished();
    void addedCategory(const QPlaceCategory &category, const QString &parentId);
    void updatedCategory(const QPlaceCategory &category, const QString &parentId);
    void removedCategory(const QString &categoryId, const QString &parentId);
    void connectNotificationSignals();

private:
    QStringList populateCategories(QPlaceManager *, const QPlaceCategory &parent);
    QModelIndex index(const QString &categoryId) const;
    int rowToAddChild(PlaceCategoryNode *, const QPlaceCategory &category);
    void updateLayout();

    QPlaceReply *m_response;

    QDeclarativeGeoServiceProvider *m_plugin;
    bool m_hierarchical;
    bool m_complete;
    Status m_status;
    QString m_errorString;

    QHash<QString, PlaceCategoryNode *> m_categoriesTree;
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QDeclarativeSupportedCategoriesModel)

#endif // QDECLARATIVESUPPORTEDCATEGORIESMODEL_H
