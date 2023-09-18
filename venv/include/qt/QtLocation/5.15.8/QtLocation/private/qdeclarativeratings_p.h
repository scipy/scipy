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

#ifndef QDECLARATIVERATINGS_P_H
#define QDECLARATIVERATINGS_P_H

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
#include <QtLocation/qplaceratings.h>
#include <QtQml/qqml.h>

#include <QObject>

QT_BEGIN_NAMESPACE

class Q_LOCATION_PRIVATE_EXPORT QDeclarativeRatings : public QObject
{
    Q_OBJECT

    Q_PROPERTY(QPlaceRatings ratings READ ratings WRITE setRatings)
    Q_PROPERTY(qreal average READ average WRITE setAverage NOTIFY averageChanged)
    Q_PROPERTY(qreal maximum READ maximum WRITE setMaximum NOTIFY maximumChanged)
    Q_PROPERTY(int count READ count WRITE setCount NOTIFY countChanged)

public:
    explicit QDeclarativeRatings(QObject *parent = 0);
    explicit QDeclarativeRatings(const QPlaceRatings &src, QObject *parent = 0);
    ~QDeclarativeRatings();

    QPlaceRatings ratings() const;
    void setRatings(const QPlaceRatings &src);

    qreal average() const;
    void setAverage(qreal average);

    qreal maximum() const;
    void setMaximum(qreal max);

    int count() const;
    void setCount(int count);

Q_SIGNALS:
    void averageChanged();
    void maximumChanged();
    void countChanged();

private:
    QPlaceRatings m_ratings;
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QDeclarativeRatings)

#endif // QDECLARATIVERATING_P_H
