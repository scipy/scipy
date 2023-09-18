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

#ifndef QPLACECATEGORY_H
#define QPLACECATEGORY_H

#include <QtLocation/qlocationglobal.h>

#include <QtCore/QSharedDataPointer>
#include <QtCore/QString>
#include <QtCore/QMetaType>
#include <QtLocation/QLocation>

QT_BEGIN_NAMESPACE

class QPlaceIcon;

class QPlaceCategoryPrivate;
class Q_LOCATION_EXPORT QPlaceCategory
{
public:
    QPlaceCategory();
    QPlaceCategory(const QPlaceCategory &other);

    virtual ~QPlaceCategory();

    QPlaceCategory &operator=(const QPlaceCategory &other);

    bool operator==(const QPlaceCategory &other) const;
    bool operator!=(const QPlaceCategory &other) const {
        return !(other == *this);
    }

    QString categoryId() const;
    void setCategoryId(const QString &identifier);

    QString name() const;
    void setName(const QString &name);

    QLocation::Visibility visibility() const;
    void setVisibility(QLocation::Visibility visibility);

    QPlaceIcon icon() const;
    void setIcon(const QPlaceIcon &icon);

    bool isEmpty() const;

private:
    QSharedDataPointer<QPlaceCategoryPrivate> d;
};

Q_DECLARE_TYPEINFO(QPlaceCategory, Q_MOVABLE_TYPE);

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QPlaceCategory)

#endif // QPLACECATEGORY_H
