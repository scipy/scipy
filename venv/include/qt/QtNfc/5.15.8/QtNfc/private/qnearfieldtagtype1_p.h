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

#ifndef QNEARFIELDTAGTYPE1_H
#define QNEARFIELDTAGTYPE1_H

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

#include <QtNfc/QNearFieldTarget>

QT_BEGIN_NAMESPACE

class QNearFieldTagType1Private;

class Q_AUTOTEST_EXPORT QNearFieldTagType1 : public QNearFieldTarget
{
    Q_OBJECT

    Q_DECLARE_PRIVATE(QNearFieldTagType1)

public:
    enum WriteMode {
        EraseAndWrite,
        WriteOnly
    };
    Q_ENUM(WriteMode)

    explicit QNearFieldTagType1(QObject *parent = 0);
    ~QNearFieldTagType1();

    Type type() const { return NfcTagType1; }

    bool hasNdefMessage();
    RequestId readNdefMessages();
    RequestId writeNdefMessages(const QList<QNdefMessage> &messages);

    quint8 version();
    virtual int memorySize();

    // DIGPROTO
    virtual RequestId readIdentification();

    // static memory functions
    virtual RequestId readAll();
    virtual RequestId readByte(quint8 address);
    virtual RequestId writeByte(quint8 address, quint8 data, WriteMode mode = EraseAndWrite);

    // dynamic memory functions
    virtual RequestId readSegment(quint8 segmentAddress);
    virtual RequestId readBlock(quint8 blockAddress);
    virtual RequestId writeBlock(quint8 blockAddress, const QByteArray &data,
                                 WriteMode mode = EraseAndWrite);

protected:
    bool handleResponse(const QNearFieldTarget::RequestId &id, const QByteArray &response);

private:
    QNearFieldTagType1Private *d_ptr;
};

QT_END_NAMESPACE

#endif // QNEARFIELDTAGTYPE1_H
