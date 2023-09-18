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

#ifndef QPLACECONTACTDETAIL_H
#define QPLACECONTACTDETAIL_H

#include <QtCore/QString>
#include <QtCore/QVariant>
#include <QtCore/QSharedDataPointer>

#include <QtLocation/qlocationglobal.h>

QT_BEGIN_NAMESPACE

class QPlaceContactDetailPrivate;
class Q_LOCATION_EXPORT QPlaceContactDetail
{
public:
    static const QString Phone;
    static const QString Email;
    static const QString Website;
    static const QString Fax;

    QPlaceContactDetail();
    QPlaceContactDetail(const QPlaceContactDetail &other);
    virtual ~QPlaceContactDetail();

    QPlaceContactDetail &operator=(const QPlaceContactDetail &other);

    bool operator==(const QPlaceContactDetail &other) const;
    bool operator!=(const QPlaceContactDetail &other) const;

    QString label() const;
    void setLabel(const QString &label);

    QString value() const;
    void setValue(const QString &value);

    void clear();

private:
    QSharedDataPointer<QPlaceContactDetailPrivate> d_ptr;

};

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QPlaceContactDetail)

#endif
