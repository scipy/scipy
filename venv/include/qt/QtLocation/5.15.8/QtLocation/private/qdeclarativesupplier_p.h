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

#ifndef QDECLARATIVESUPPLIER_P_H
#define QDECLARATIVESUPPLIER_P_H

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
#include <QObject>
#include <QtCore/QUrl>
#include <QtQml/qqml.h>
#include <QtQml/QQmlParserStatus>
#include <QtLocation/qplacesupplier.h>

#include <QtLocation/private/qdeclarativeplaceicon_p.h>

QT_BEGIN_NAMESPACE

class Q_LOCATION_PRIVATE_EXPORT QDeclarativeSupplier : public QObject, public QQmlParserStatus
{
    Q_OBJECT

    Q_PROPERTY(QPlaceSupplier supplier READ supplier WRITE setSupplier)
    Q_PROPERTY(QString name READ name WRITE setName NOTIFY nameChanged)
    Q_PROPERTY(QString supplierId READ supplierId WRITE setSupplierId NOTIFY supplierIdChanged)
    Q_PROPERTY(QUrl url READ url WRITE setUrl NOTIFY urlChanged)
    Q_PROPERTY(QDeclarativePlaceIcon *icon READ icon WRITE setIcon NOTIFY iconChanged)

    Q_INTERFACES(QQmlParserStatus)

public:
    explicit QDeclarativeSupplier(QObject *parent = 0);
    explicit QDeclarativeSupplier(const QPlaceSupplier &src, QDeclarativeGeoServiceProvider *plugin, QObject *parent = 0);
    ~QDeclarativeSupplier();

    // From QQmlParserStatus
    void classBegin() { }
    void componentComplete();

    QPlaceSupplier supplier();
    void setSupplier(const QPlaceSupplier &src, QDeclarativeGeoServiceProvider *plugin = 0);

    QString name() const;
    void setName(const QString &data);
    QString supplierId() const;
    void setSupplierId(const QString &data);
    QUrl url() const;
    void setUrl(const QUrl &data);

    QDeclarativePlaceIcon *icon() const;
    void setIcon(QDeclarativePlaceIcon *icon);

Q_SIGNALS:
    void nameChanged();
    void supplierIdChanged();
    void urlChanged();
    void iconChanged();

private:
    QPlaceSupplier m_src;
    QDeclarativePlaceIcon *m_icon;
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QDeclarativeSupplier)

#endif // QDECLARATIVESUPPLIER_P_H
