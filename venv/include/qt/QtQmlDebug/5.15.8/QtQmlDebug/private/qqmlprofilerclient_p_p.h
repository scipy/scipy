/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
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

#ifndef QQMLPROFILERCLIENT_P_P_H
#define QQMLPROFILERCLIENT_P_P_H

#include "qqmldebugclient_p_p.h"
#include "qqmldebugmessageclient_p.h"
#include "qqmlenginecontrolclient_p.h"
#include "qqmlprofilerclient_p.h"
#include "qqmlprofilertypedevent_p.h"
#include "qqmlprofilerclientdefinitions_p.h"

#include <QtCore/qqueue.h>
#include <QtCore/qstack.h>

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

class QQmlProfilerClientPrivate : public QQmlDebugClientPrivate {
    Q_DECLARE_PUBLIC(QQmlProfilerClient)
public:
    QQmlProfilerClientPrivate(QQmlDebugConnection *connection,
                              QQmlProfilerEventReceiver *eventReceiver)
        : QQmlDebugClientPrivate(QLatin1String("CanvasFrameRate"), connection)
        , eventReceiver(eventReceiver)
        , engineControl(new QQmlEngineControlClient(connection))
        , maximumTime(0)
        , recording(false)
        , requestedFeatures(0)
        , recordedFeatures(0)
        , flushInterval(0)
    {
    }

    virtual ~QQmlProfilerClientPrivate() override;

    void sendRecordingStatus(int engineId);
    bool updateFeatures(ProfileFeature feature);
    int resolveType(const QQmlProfilerTypedEvent &type);
    int resolveStackTop();
    void forwardEvents(const QQmlProfilerEvent &last);
    void forwardDebugMessages(qint64 untilTimestamp);
    void processCurrentEvent();
    void finalize();

    QQmlProfilerEventReceiver *eventReceiver;
    QScopedPointer<QQmlEngineControlClient> engineControl;
    QScopedPointer<QQmlDebugMessageClient> messageClient;
    qint64 maximumTime;
    bool recording;
    quint64 requestedFeatures;
    quint64 recordedFeatures;
    quint32 flushInterval;

    // Reuse the same event, so that we don't have to constantly reallocate all the data.
    QQmlProfilerTypedEvent currentEvent;
    QHash<QQmlProfilerEventType, int> eventTypeIds;
    QHash<qint64, int> serverTypeIds;
    QStack<QQmlProfilerTypedEvent> rangesInProgress;
    QQueue<QQmlProfilerEvent> pendingMessages;
    QQueue<QQmlProfilerEvent> pendingDebugMessages;

    QList<int> trackedEngines;
};

QT_END_NAMESPACE

#endif // QQMLPROFILERCLIENT_P_P_H
