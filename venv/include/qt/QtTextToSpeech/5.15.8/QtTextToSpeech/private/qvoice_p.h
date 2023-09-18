/****************************************************************************
**
** Copyright (C) 2015 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the Qt Speech module of the Qt Toolkit.
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



#ifndef QVOICE_P_H
#define QVOICE_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of other Qt classes.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <qvoice.h>

#include <QString>
#include <QCoreApplication>
#include <QVariant>

QT_BEGIN_NAMESPACE

class QVoicePrivate : public QSharedData
{
public:
    QVoicePrivate();
    QVoicePrivate(const QVoicePrivate &other);
    QVoicePrivate(const QString &n, QVoice::Gender g, QVoice::Age a, const QVariant &d);
    ~QVoicePrivate() {}

    QString name;
    QVoice::Gender gender;
    QVoice::Age age;
    // Various data depending on the platform:
    // On OS X the VoiceIdentifier is stored.
    // On unix the synthesizer (output module) is stored.
    QVariant data;
};

QVoicePrivate::QVoicePrivate()
{
    name = qApp->translate("QVoice", "Default");
    age = QVoice::Other;
    gender = QVoice::Unknown;
    data = QVariant();
}

QVoicePrivate::QVoicePrivate(const QVoicePrivate &other)
 :QSharedData(other), name(other.name), gender(other.gender), age(other.age), data(other.data)
{
}

QVoicePrivate::QVoicePrivate(const QString &n, QVoice::Gender g, QVoice::Age a, const QVariant &d)
    :name(n), gender(g), age(a), data(d)
{
}

QT_END_NAMESPACE

#endif

