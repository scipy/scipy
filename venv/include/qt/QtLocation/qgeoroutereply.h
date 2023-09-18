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

#ifndef QGEOROUTEREPLY_H
#define QGEOROUTEREPLY_H

#include <QtLocation/QGeoRoute>

#include <QtCore/QList>
#include <QtCore/QObject>

QT_BEGIN_NAMESPACE

class QGeoRouteRequest;
class QGeoRouteReplyPrivate;

class Q_LOCATION_EXPORT QGeoRouteReply : public QObject
{
    Q_OBJECT
public:
    enum Error {
        NoError,
        EngineNotSetError,
        CommunicationError,
        ParseError,
        UnsupportedOptionError,
        UnknownError
    };

    explicit QGeoRouteReply(Error error, const QString &errorString, QObject *parent = nullptr);
    virtual ~QGeoRouteReply();

    bool isFinished() const;
    Error error() const;
    QString errorString() const;

    QGeoRouteRequest request() const;
    QList<QGeoRoute> routes() const;

    virtual void abort();

Q_SIGNALS:
    void finished();
    void aborted();
    void error(QGeoRouteReply::Error error, const QString &errorString = QString());

protected:
    explicit QGeoRouteReply(const QGeoRouteRequest &request, QObject *parent = nullptr);

    void setError(Error error, const QString &errorString);
    void setFinished(bool finished);

    void setRoutes(const QList<QGeoRoute> &routes);
    void addRoutes(const QList<QGeoRoute> &routes);

private:
    QGeoRouteReplyPrivate *d_ptr;
    Q_DISABLE_COPY(QGeoRouteReply)
};

QT_END_NAMESPACE

#endif
