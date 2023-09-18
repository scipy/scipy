/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtNfc module of the Qt Toolkit.
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

#ifndef QNEARFIELDMANAGER_H
#define QNEARFIELDMANAGER_H

#include <QtCore/QObject>
#include <QtNfc/qtnfcglobal.h>
#include <QtNfc/QNearFieldTarget>
#include <QtNfc/QNdefRecord>
#include <QtNfc/QNdefFilter>

QT_BEGIN_NAMESPACE

class QNearFieldManagerPrivate;
class Q_NFC_EXPORT QNearFieldManager : public QObject
{
    Q_OBJECT

    Q_DECLARE_PRIVATE(QNearFieldManager)

public:
    enum class AdapterState {
        Offline = 1,
        TurningOn = 2,
        Online = 3,
        TurningOff = 4
    };
    Q_ENUM(AdapterState)
    enum TargetAccessMode {
        NoTargetAccess = 0x00,
        NdefReadTargetAccess = 0x01,
        NdefWriteTargetAccess = 0x02,
        TagTypeSpecificTargetAccess = 0x04
    };
    Q_ENUM(TargetAccessMode)
    Q_DECLARE_FLAGS(TargetAccessModes, TargetAccessMode)

    explicit QNearFieldManager(QObject *parent = nullptr);
    explicit QNearFieldManager(QNearFieldManagerPrivate *backend, QObject *parent = nullptr);
    ~QNearFieldManager();

    bool isAvailable() const;
    bool isSupported() const;

    void setTargetAccessModes(TargetAccessModes accessModes);
    TargetAccessModes targetAccessModes() const;

    bool startTargetDetection();
    void stopTargetDetection();

    //TODO Qt 6 Consider removal of this registration mechanism
    //None of the currently supported platforms supports the feature
    //or in fact the implementation (on Android) is not what the
    //function is supposed to do.
    int registerNdefMessageHandler(QObject *object, const char *method);
    int registerNdefMessageHandler(QNdefRecord::TypeNameFormat typeNameFormat,
                                   const QByteArray &type,
                                   QObject *object, const char *method);
    int registerNdefMessageHandler(const QNdefFilter &filter,
                                   QObject *object, const char *method);

    bool unregisterNdefMessageHandler(int handlerId);

Q_SIGNALS:
    void adapterStateChanged(QNearFieldManager::AdapterState state);
    void targetDetected(QNearFieldTarget *target);
    void targetLost(QNearFieldTarget *target);

private:
    QNearFieldManagerPrivate *d_ptr;
};

Q_DECLARE_OPERATORS_FOR_FLAGS(QNearFieldManager::TargetAccessModes)

QT_END_NAMESPACE

#endif // QNEARFIELDMANAGER_H
