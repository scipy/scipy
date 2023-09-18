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

#ifndef UNSUPPORTEDREPLIES_P_H
#define UNSUPPORTEDREPLIES_P_H

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
#include "qplacedetailsreply.h"
#include "qplacecontentreply.h"
#include "qplacesearchreply.h"
#include "qplacesearchsuggestionreply.h"
#include "qplaceidreply.h"

#include "qplacematchreply.h"
#include "qplacemanagerengine.h"

class Q_LOCATION_PRIVATE_EXPORT QPlaceDetailsReplyUnsupported : public QPlaceDetailsReply
{
    Q_OBJECT

public:
    QPlaceDetailsReplyUnsupported(QPlaceManagerEngine *parent)
    :   QPlaceDetailsReply(parent)
    {
        setError(QPlaceReply::UnsupportedError,
                 QStringLiteral("Getting place details is not supported."));
        setFinished(true);
        QMetaObject::invokeMethod(this, "error", Qt::QueuedConnection,
                                  Q_ARG(QPlaceReply::Error, error()),
                                  Q_ARG(QString, errorString()));
        QMetaObject::invokeMethod(parent, "error", Qt::QueuedConnection,
                                  Q_ARG(QPlaceReply *, this),
                                  Q_ARG(QPlaceReply::Error, error()),
                                  Q_ARG(QString, errorString()));
        QMetaObject::invokeMethod(this, "finished", Qt::QueuedConnection);
        QMetaObject::invokeMethod(parent, "finished", Qt::QueuedConnection,
                                  Q_ARG(QPlaceReply *, this));
    }
};

class Q_LOCATION_PRIVATE_EXPORT QPlaceContentReplyUnsupported : public QPlaceContentReply
{
    Q_OBJECT

public:
    QPlaceContentReplyUnsupported(QPlaceManagerEngine *parent)
    :   QPlaceContentReply(parent)
    {
        setError(QPlaceReply::UnsupportedError,
                 QStringLiteral("Place content is not supported."));
        setFinished(true);
        QMetaObject::invokeMethod(this, "error", Qt::QueuedConnection,
                                  Q_ARG(QPlaceReply::Error, error()),
                                  Q_ARG(QString, errorString()));
        QMetaObject::invokeMethod(parent, "error", Qt::QueuedConnection,
                                  Q_ARG(QPlaceReply *, this),
                                  Q_ARG(QPlaceReply::Error, error()),
                                  Q_ARG(QString, errorString()));
        QMetaObject::invokeMethod(this, "finished", Qt::QueuedConnection);
        QMetaObject::invokeMethod(parent, "finished", Qt::QueuedConnection,
                                  Q_ARG(QPlaceReply *, this));
    }
};

class Q_LOCATION_PRIVATE_EXPORT QPlaceSearchReplyUnsupported : public QPlaceSearchReply
{
    Q_OBJECT

public:
    QPlaceSearchReplyUnsupported(QPlaceReply::Error errorCode, const QString &message,
                                 QPlaceManagerEngine *parent)
    :   QPlaceSearchReply(parent)
    {
        setError(errorCode, message);
        setFinished(true);
        QMetaObject::invokeMethod(this, "error", Qt::QueuedConnection,
                                  Q_ARG(QPlaceReply::Error, error()),
                                  Q_ARG(QString, errorString()));
        QMetaObject::invokeMethod(parent, "error", Qt::QueuedConnection,
                                  Q_ARG(QPlaceReply *, this),
                                  Q_ARG(QPlaceReply::Error, error()),
                                  Q_ARG(QString, errorString()));
        QMetaObject::invokeMethod(this, "finished", Qt::QueuedConnection);
        QMetaObject::invokeMethod(parent, "finished", Qt::QueuedConnection,
                                  Q_ARG(QPlaceReply *, this));
    }
};

class Q_LOCATION_PRIVATE_EXPORT QPlaceSearchSuggestionReplyUnsupported : public QPlaceSearchSuggestionReply
{
    Q_OBJECT

public:
    QPlaceSearchSuggestionReplyUnsupported(QPlaceManagerEngine *parent)
    :   QPlaceSearchSuggestionReply(parent)
    {
        setError(QPlaceReply::UnsupportedError,
                 QStringLiteral("Place search suggestions are not supported."));
        setFinished(true);
        QMetaObject::invokeMethod(this, "error", Qt::QueuedConnection,
                                  Q_ARG(QPlaceReply::Error, error()),
                                  Q_ARG(QString, errorString()));
        QMetaObject::invokeMethod(parent, "error", Qt::QueuedConnection,
                                  Q_ARG(QPlaceReply *, this),
                                  Q_ARG(QPlaceReply::Error, error()),
                                  Q_ARG(QString, errorString()));
        QMetaObject::invokeMethod(this, "finished", Qt::QueuedConnection);
        QMetaObject::invokeMethod(parent, "finished", Qt::QueuedConnection,
                                  Q_ARG(QPlaceReply *, this));
    }
};

class Q_LOCATION_PRIVATE_EXPORT QPlaceIdReplyUnsupported : public QPlaceIdReply
{
    Q_OBJECT

public:
    QPlaceIdReplyUnsupported(const QString &message, QPlaceIdReply::OperationType type,
                             QPlaceManagerEngine *parent)
    :   QPlaceIdReply(type, parent)
    {
        setError(QPlaceReply::UnsupportedError, message);
        setFinished(true);
        QMetaObject::invokeMethod(this, "error", Qt::QueuedConnection,
                                  Q_ARG(QPlaceReply::Error, error()),
                                  Q_ARG(QString, errorString()));
        QMetaObject::invokeMethod(parent, "error", Qt::QueuedConnection,
                                  Q_ARG(QPlaceReply *, this),
                                  Q_ARG(QPlaceReply::Error, error()),
                                  Q_ARG(QString, errorString()));
        QMetaObject::invokeMethod(this, "finished", Qt::QueuedConnection);
        QMetaObject::invokeMethod(parent, "finished", Qt::QueuedConnection,
                                  Q_ARG(QPlaceReply *, this));
    }
};

class Q_LOCATION_PRIVATE_EXPORT QPlaceReplyUnsupported : public QPlaceReply
{
    Q_OBJECT

public:
    QPlaceReplyUnsupported(const QString &message, QPlaceManagerEngine *parent)
    :   QPlaceReply(parent)
    {
        setError(QPlaceReply::UnsupportedError, message);
        setFinished(true);
        QMetaObject::invokeMethod(this, "error", Qt::QueuedConnection,
                                  Q_ARG(QPlaceReply::Error, error()),
                                  Q_ARG(QString, errorString()));
        QMetaObject::invokeMethod(parent, "error", Qt::QueuedConnection,
                                  Q_ARG(QPlaceReply *, this),
                                  Q_ARG(QPlaceReply::Error, error()),
                                  Q_ARG(QString, errorString()));
        QMetaObject::invokeMethod(this, "finished", Qt::QueuedConnection);
        QMetaObject::invokeMethod(parent, "finished", Qt::QueuedConnection,
                                  Q_ARG(QPlaceReply *, this));
    }
};

class Q_LOCATION_PRIVATE_EXPORT QPlaceMatchReplyUnsupported : public QPlaceMatchReply
{
    Q_OBJECT

public:
    QPlaceMatchReplyUnsupported(QPlaceManagerEngine *parent)
    :   QPlaceMatchReply(parent)
    {
        setError(QPlaceReply::UnsupportedError,
                 QStringLiteral("Place matching is not supported."));
        setFinished(true);
        QMetaObject::invokeMethod(this, "error", Qt::QueuedConnection,
                                  Q_ARG(QPlaceReply::Error, error()),
                                  Q_ARG(QString, errorString()));
        QMetaObject::invokeMethod(parent, "error", Qt::QueuedConnection,
                                  Q_ARG(QPlaceReply *, this),
                                  Q_ARG(QPlaceReply::Error, error()),
                                  Q_ARG(QString, errorString()));
        QMetaObject::invokeMethod(this, "finished", Qt::QueuedConnection);
        QMetaObject::invokeMethod(parent, "finished", Qt::QueuedConnection,
                                  Q_ARG(QPlaceReply *, this));
    }
};

#endif
