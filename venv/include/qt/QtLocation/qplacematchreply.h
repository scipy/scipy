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

#ifndef QPLACEMATCHREPLY_H
#define QPLACEMATCHREPLY_H

#include <QtLocation/QPlaceReply>
#include <QtLocation/QPlaceMatchRequest>

QT_BEGIN_NAMESPACE

class QPlaceMatchReplyPrivate;
class Q_LOCATION_EXPORT QPlaceMatchReply : public QPlaceReply
{
    Q_OBJECT
public:
    explicit QPlaceMatchReply(QObject *parent = nullptr);
    ~QPlaceMatchReply();

    QPlaceReply::Type type() const;

    QList<QPlace> places() const;
    QPlaceMatchRequest request() const;

protected:
    void setPlaces(const QList<QPlace> &results);
    void setRequest(const QPlaceMatchRequest &request);
private:
    Q_DISABLE_COPY(QPlaceMatchReply)
    Q_DECLARE_PRIVATE(QPlaceMatchReply)
};

QT_END_NAMESPACE

#endif
