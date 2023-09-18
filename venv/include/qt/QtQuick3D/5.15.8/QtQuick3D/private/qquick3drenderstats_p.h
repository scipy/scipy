/****************************************************************************
**
** Copyright (C) 2019 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of Qt Quick 3D.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QQUICK3DRENDERSTATS_H
#define QQUICK3DRENDERSTATS_H

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

#include <QtQuick3D/qtquick3dglobal.h>
#include <QtCore/qobject.h>
#include <QtQuick3DRuntimeRender/private/qssgrendercontextcore_p.h>

QT_BEGIN_NAMESPACE

class Q_QUICK3D_EXPORT QQuick3DRenderStats : public QObject
{
    Q_OBJECT
    Q_PROPERTY(int fps READ fps NOTIFY fpsChanged)
    Q_PROPERTY(float frameTime READ frameTime NOTIFY frameTimeChanged)
    Q_PROPERTY(float renderTime READ renderTime NOTIFY renderTimeChanged)
    Q_PROPERTY(float syncTime READ syncTime NOTIFY syncTimeChanged)
    Q_PROPERTY(float maxFrameTime READ maxFrameTime NOTIFY maxFrameTimeChanged)

public:
    QQuick3DRenderStats(QObject *parent = nullptr);

    int fps() const;
    float frameTime() const;
    float renderTime() const;
    float syncTime() const;
    float maxFrameTime() const;

    void startSync();
    void endSync(bool dump = false);
    void startRender();
    void endRender(bool dump = false);

Q_SIGNALS:
    void fpsChanged();
    void frameTimeChanged();
    void renderTimeChanged();
    void syncTimeChanged();
    void maxFrameTimeChanged();

private:
    float timestamp() const;

    QElapsedTimer m_frameTimer;
    int m_frameCount = 0;
    float m_secTimer = 0;
    float m_notifyTimer = 0;
    float m_renderStartTime = 0;
    float m_syncStartTime = 0;

    float m_internalMaxFrameTime = 0;
    float m_notifiedFrameTime = 0;
    float m_notifiedRenderTime = 0;
    float m_notifiedSyncTime = 0;

    int m_fps = 0;
    float m_frameTime = 0;
    float m_renderTime = 0;
    float m_syncTime = 0;
    float m_maxFrameTime = 0;
};

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QQuick3DRenderStats *)

#endif // QQUICK3DRENDERSTATS_H
