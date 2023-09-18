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

#ifndef QGEOCODEREPLY_H
#define QGEOCODEREPLY_H

#include <QtCore/QObject>
#include <QtCore/QList>
#include <QtPositioning/QGeoLocation>

#include <QtLocation/qlocationglobal.h>

QT_BEGIN_NAMESPACE

class QGeoShape;
class QGeoCodeReplyPrivate;

class Q_LOCATION_EXPORT QGeoCodeReply : public QObject
{
    Q_OBJECT

public:
    enum Error {
        NoError,
        EngineNotSetError,
        CommunicationError,
        ParseError,
        UnsupportedOptionError,
        CombinationError,
        UnknownError
    };

    explicit QGeoCodeReply(Error error, const QString &errorString, QObject *parent = nullptr);
    virtual ~QGeoCodeReply();

    bool isFinished() const;
    Error error() const;
    QString errorString() const;

    QGeoShape viewport() const;
    QList<QGeoLocation> locations() const;

    int limit() const;
    int offset() const;

    virtual void abort();

Q_SIGNALS:
    void finished();
    void aborted();
    void error(QGeoCodeReply::Error error, const QString &errorString = QString());

protected:
    explicit QGeoCodeReply(QObject *parent = nullptr);
    explicit QGeoCodeReply(QGeoCodeReplyPrivate &dd, QObject *parent = nullptr);

    void setError(Error error, const QString &errorString);
    void setFinished(bool finished);

    void setViewport(const QGeoShape &viewport);
    void addLocation(const QGeoLocation &location);
    void setLocations(const QList<QGeoLocation> &locations);

    void setLimit(int limit);
    void setOffset(int offset);

private:
    QGeoCodeReplyPrivate *d_ptr;
    Q_DISABLE_COPY(QGeoCodeReply)
    friend class QGeoCodeReplyPrivate;
};

QT_END_NAMESPACE

#endif
