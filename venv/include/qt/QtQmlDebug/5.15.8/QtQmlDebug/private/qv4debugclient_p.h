/****************************************************************************
**
** Copyright (C) 2018 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtQml module of the Qt Toolkit.
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

#ifndef QV4DEBUGCLIENT_P_H
#define QV4DEBUGCLIENT_P_H

#include <QtQmlDebug/private/qqmldebugclient_p.h>
#include <QtCore/qjsonvalue.h>

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

QT_BEGIN_NAMESPACE

class QV4DebugClientPrivate;
class QV4DebugClient : public QQmlDebugClient
{
    Q_OBJECT
    Q_DECLARE_PRIVATE(QV4DebugClient)

public:
    enum StepAction
    {
        Continue,
        In,
        Out,
        Next
    };

    enum Exception
    {
        All,
        Uncaught
    };

    struct Response
    {
        QString command;
        QJsonValue body;
    };

    QV4DebugClient(QQmlDebugConnection *connection);

    void connect();
    void disconnect();

    void interrupt();
    void continueDebugging(StepAction stepAction);
    void evaluate(const QString &expr, int frame = -1, int context = -1);
    void lookup(const QList<int> &handles, bool includeSource = false);
    void backtrace(int fromFrame = -1, int toFrame = -1, bool bottom = false);
    void frame(int number = -1);
    void scope(int number = -1, int frameNumber = -1);
    void scripts(int types = 4, const QList<int> &ids = QList<int>(), bool includeSource = false);
    void setBreakpoint(const QString &target, int line = -1, int column = -1, bool enabled = true,
                       const QString &condition = QString(), int ignoreCount = -1);
    void clearBreakpoint(int breakpoint);
    void changeBreakpoint(int breakpoint, bool enabled);
    void setExceptionBreak(Exception type, bool enabled = false);
    void version();

    Response response() const;

protected:
    void messageReceived(const QByteArray &data) override;

signals:
    void connected();
    void interrupted();
    void result();
    void failure();
    void stopped();
};

QT_END_NAMESPACE

#endif // QV4DEBUGCLIENT_P_H
