/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtScxml module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl-3.0.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or (at your option) the GNU General
** Public license version 3 or any later version approved by the KDE Free
** Qt Foundation. The licenses are as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL2 and LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-2.0.html and
** https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QSCXMLEVENT_H
#define QSCXMLEVENT_H

#include <QtScxml/qscxmlglobals.h>

#include <QtCore/qstringlist.h>
#include <QtCore/qvariant.h>

QT_BEGIN_NAMESPACE

class QScxmlEventPrivate;

class Q_SCXML_EXPORT QScxmlEvent
{
    Q_GADGET
    Q_PROPERTY(QString name READ name WRITE setName)
    Q_PROPERTY(EventType eventType READ eventType WRITE setEventType)
    Q_PROPERTY(QString scxmlType READ scxmlType)
    Q_PROPERTY(QString sendId READ sendId WRITE setSendId)
    Q_PROPERTY(QString origin READ origin WRITE setOrigin)
    Q_PROPERTY(QString originType READ originType WRITE setOriginType)
    Q_PROPERTY(QString invokeId READ invokeId WRITE setInvokeId)
    Q_PROPERTY(int delay READ delay WRITE setDelay)
    Q_PROPERTY(QVariant data READ data WRITE setData)
    Q_PROPERTY(bool errorEvent READ isErrorEvent)
    Q_PROPERTY(QString errorMessage READ errorMessage WRITE setErrorMessage)

public:
    QScxmlEvent();
    ~QScxmlEvent();

    QScxmlEvent &operator=(const QScxmlEvent &other);
    QScxmlEvent(const QScxmlEvent &other);

    enum EventType {
        PlatformEvent,
        InternalEvent,
        ExternalEvent
    };
    Q_ENUM(EventType)

    QString name() const;
    void setName(const QString &name);

    EventType eventType() const;
    void setEventType(const EventType &type);

    QString scxmlType() const;

    QString sendId() const;
    void setSendId(const QString &sendId);

    QString origin() const;
    void setOrigin(const QString &origin);

    QString originType() const;
    void setOriginType(const QString &originType);

    QString invokeId() const;
    void setInvokeId(const QString &invokeId);

    int delay() const;
    void setDelay(int delayInMiliSecs);

    Q_INVOKABLE void clear();

    QVariant data() const;
    void setData(const QVariant &data);

    bool isErrorEvent() const;
    QString errorMessage() const;
    void setErrorMessage(const QString &message);

private:
    QScxmlEventPrivate *d;

};

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QScxmlEvent)

#endif // QSCXMLEVENT_H
