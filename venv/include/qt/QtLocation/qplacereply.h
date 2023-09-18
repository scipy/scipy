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

#ifndef QPLACEREPLY_H
#define QPLACEREPLY_H

#include <QtCore/QMetaType>
#include <QtCore/QObject>
#include <QtLocation/qlocationglobal.h>

QT_BEGIN_NAMESPACE

class QPlaceReplyPrivate;
class Q_LOCATION_EXPORT QPlaceReply : public QObject
{
    Q_OBJECT
public:
    enum Error {
        NoError,
        PlaceDoesNotExistError,
        CategoryDoesNotExistError,
        CommunicationError,
        ParseError,
        PermissionsError,
        UnsupportedError,
        BadArgumentError,
        CancelError,
        UnknownError
    };

    enum Type {
        Reply,
        DetailsReply,
        SearchReply,
        SearchSuggestionReply,
        ContentReply,
        IdReply,
        MatchReply
    };

    explicit QPlaceReply(QObject *parent = nullptr);
    ~QPlaceReply();

    bool isFinished() const;

    virtual Type type() const;

    QString errorString() const;
    QPlaceReply::Error error() const;

public Q_SLOTS:
    virtual void abort();

Q_SIGNALS:
    void finished();
    void contentUpdated();
    void aborted();
    void error(QPlaceReply::Error error, const QString &errorString = QString());

protected:
    explicit QPlaceReply(QPlaceReplyPrivate *, QObject *parent = nullptr);
    void setFinished(bool finished);
    void setError(QPlaceReply::Error error, const QString &errorString);
    QPlaceReplyPrivate *d_ptr;

private:
    Q_DISABLE_COPY(QPlaceReply)
};

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QPlaceReply::Error)
Q_DECLARE_METATYPE(QPlaceReply *)

#endif // QPLACEREPLY_H
