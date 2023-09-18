/****************************************************************************
**
** Copyright (C) 2016 Kurt Pattyn <pattyn.kurt@gmail.com>.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtWebSockets module of the Qt Toolkit.
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

#ifndef QWEBSOCKETDATAPROCESSOR_P_H
#define QWEBSOCKETDATAPROCESSOR_P_H

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

#include <QtCore/QObject>
#include <QtCore/QByteArray>
#include <QtCore/QString>
#include <QtCore/QTextCodec>
#include <QTimer>
#include "qwebsocketframe_p.h"
#include "qwebsocketprotocol.h"
#include "qwebsocketprotocol_p.h"

QT_BEGIN_NAMESPACE

class QIODevice;
class QWebSocketFrame;

const quint64 MAX_MESSAGE_SIZE_IN_BYTES = std::numeric_limits<int>::max() - 1;

class Q_AUTOTEST_EXPORT QWebSocketDataProcessor : public QObject
{
    Q_OBJECT
    Q_DISABLE_COPY(QWebSocketDataProcessor)

public:
    explicit QWebSocketDataProcessor(QObject *parent = nullptr);
    ~QWebSocketDataProcessor() override;

    void setMaxAllowedFrameSize(quint64 maxAllowedFrameSize);
    quint64 maxAllowedFrameSize() const;
    void setMaxAllowedMessageSize(quint64 maxAllowedMessageSize);
    quint64 maxAllowedMessageSize() const;
    static quint64 maxMessageSize();
    static quint64 maxFrameSize();

Q_SIGNALS:
    void pingReceived(const QByteArray &data);
    void pongReceived(const QByteArray &data);
    void closeReceived(QWebSocketProtocol::CloseCode closeCode, const QString &closeReason);
    void textFrameReceived(const QString &frame, bool lastFrame);
    void binaryFrameReceived(const QByteArray &frame, bool lastFrame);
    void textMessageReceived(const QString &message);
    void binaryMessageReceived(const QByteArray &message);
    void errorEncountered(QWebSocketProtocol::CloseCode code, const QString &description);

public Q_SLOTS:
    bool process(QIODevice *pIoDevice);
    void clear();

private:
    enum
    {
        PS_READ_HEADER,
        PS_READ_PAYLOAD_LENGTH,
        PS_READ_BIG_PAYLOAD_LENGTH,
        PS_READ_MASK,
        PS_READ_PAYLOAD,
        PS_DISPATCH_RESULT
    } m_processingState;

    bool m_isFinalFrame;
    bool m_isFragmented;
    QWebSocketProtocol::OpCode m_opCode;
    bool m_isControlFrame;
    bool m_hasMask;
    quint32 m_mask;
    QByteArray m_binaryMessage;
    QString m_textMessage;
    quint64 m_payloadLength;
    QTextCodec::ConverterState *m_pConverterState;
    QTextCodec *m_pTextCodec;
    QWebSocketFrame frame;
    QTimer *m_waitTimer;
    quint64 m_maxAllowedMessageSize = MAX_MESSAGE_SIZE_IN_BYTES;

    bool processControlFrame(const QWebSocketFrame &frame);
    void timeout();
};

QT_END_NAMESPACE

#endif // QWEBSOCKETDATAPROCESSOR_P_H
