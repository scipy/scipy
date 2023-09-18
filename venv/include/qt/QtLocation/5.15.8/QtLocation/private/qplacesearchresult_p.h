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

#ifndef QPLACESEARCHRESULT_P_H
#define QPLACESEARCHRESULT_P_H

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

#include "qplacesearchresult.h"
#include "qplacesearchrequest.h"

#include <QSharedData>
#include <QtLocation/QPlaceIcon>

QT_BEGIN_NAMESPACE

// defines must be in sync with class below
#define Q_IMPLEMENT_SEARCHRESULT_D_FUNC(Class) \
    Class##Private *Class::d_func() { return reinterpret_cast<Class##Private *>(d_ptr.data()); } \
    const Class##Private *Class::d_func() const { return reinterpret_cast<const Class##Private *>(d_ptr.constData()); } \

#define Q_IMPLEMENT_SEARCHRESULT_COPY_CTOR(Class) \
    Class::Class(const QPlaceSearchResult &other) : QPlaceSearchResult() { Class##Private::copyIfPossible(d_ptr, other); }

#define Q_DEFINE_SEARCHRESULT_PRIVATE_HELPER(Class, ResultType) \
    virtual QPlaceSearchResultPrivate *clone() const override { return new Class##Private(*this); } \
    virtual QPlaceSearchResult::SearchResultType type() const override {return ResultType;} \
    static void copyIfPossible(QSharedDataPointer<QPlaceSearchResultPrivate> &d_ptr, const QPlaceSearchResult &other) \
    { \
        if (other.type() == ResultType) \
            d_ptr = extract_d(other); \
        else \
            d_ptr = new Class##Private; \
    }

class QPlaceSearchResultPrivate : public QSharedData
{
public:
    QPlaceSearchResultPrivate() {}
    virtual ~QPlaceSearchResultPrivate() {}

    virtual bool compare(const QPlaceSearchResultPrivate *other) const;

    static const QSharedDataPointer<QPlaceSearchResultPrivate>
            &extract_d(const QPlaceSearchResult &other) { return other.d_ptr; }

    virtual QPlaceSearchResultPrivate *clone() const { return new QPlaceSearchResultPrivate(*this); }
    virtual QPlaceSearchResult::SearchResultType type() const  { return QPlaceSearchResult::UnknownSearchResult; }
    static void copyIfPossible(QSharedDataPointer<QPlaceSearchResultPrivate> &d_ptr, const QPlaceSearchResult &other)
    {
        if (other.type() == QPlaceSearchResult::UnknownSearchResult)
            d_ptr = extract_d(other);
        else
            d_ptr = new QPlaceSearchResultPrivate;
    }

    QString title;
    QPlaceIcon icon;
};

template<> QPlaceSearchResultPrivate *QSharedDataPointer<QPlaceSearchResultPrivate>::clone();

QT_END_NAMESPACE

#endif // QPLACESEARCHRESULT_P_H
