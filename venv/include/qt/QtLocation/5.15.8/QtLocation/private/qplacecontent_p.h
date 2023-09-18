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

#ifndef QPLACECONTENT_P_H
#define QPLACECONTENT_P_H

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

#include "qplacecontent.h"
#include "qplacesupplier.h"
#include "qplaceuser.h"

#include <QtCore/QSharedData>
#include <QtCore/QString>
#include <QtCore/QUrl>

QT_BEGIN_NAMESPACE


#define Q_IMPLEMENT_CONTENT_D_FUNC(Class) \
    Class##Private *Class::d_func() { return reinterpret_cast<Class##Private *>(d_ptr.data()); } \
    const Class##Private *Class::d_func() const { return reinterpret_cast<const Class##Private *>(d_ptr.constData()); } \

#define Q_IMPLEMENT_CONTENT_COPY_CTOR(Class) \
    Class::Class(const QPlaceContent &other) : QPlaceContent() { Class##Private::copyIfPossible(d_ptr, other); }

#define Q_DEFINE_CONTENT_PRIVATE_HELPER(Class, ContentType) \
    virtual QPlaceContentPrivate *clone() const { return new Class##Private(*this); } \
    virtual QPlaceContent::Type type() const {return ContentType;} \
    static void copyIfPossible(QSharedDataPointer<QPlaceContentPrivate> &d_ptr, const QPlaceContent &other) \
    { \
        if (other.type() == ContentType) \
            d_ptr = extract_d(other); \
        else \
            d_ptr = new Class##Private; \
    }

class QPlaceContentPrivate : public QSharedData
{
public:
    QPlaceContentPrivate(){}
    virtual ~QPlaceContentPrivate(){}

    virtual bool compare(const QPlaceContentPrivate *other) const;
    virtual QPlaceContentPrivate *clone() const = 0;
    virtual QPlaceContent::Type type() const = 0;

    /* Helper functions for C++ protection rules */
    static const QSharedDataPointer<QPlaceContentPrivate> &extract_d(const QPlaceContent &other) {return other.d_ptr;}

    QPlaceSupplier supplier;
    QPlaceUser user;
    QString attribution;
};

template<> QPlaceContentPrivate *QSharedDataPointer<QPlaceContentPrivate>::clone();

QT_END_NAMESPACE

#endif

