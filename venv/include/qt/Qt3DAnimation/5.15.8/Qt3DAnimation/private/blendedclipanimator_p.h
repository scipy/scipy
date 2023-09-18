/****************************************************************************
**
** Copyright (C) 2017 Klaralvdalens Datakonsult AB (KDAB).
** Contact: http://www.qt-project.org/legal
**
** This file is part of the Qt3D module of the Qt Toolkit.
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

#ifndef QT3DANIMATION_ANIMATION_BLENDEDCLIPANIMATOR_P_H
#define QT3DANIMATION_ANIMATION_BLENDEDCLIPANIMATOR_P_H

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

#include <Qt3DAnimation/private/backendnode_p.h>
#include <Qt3DAnimation/private/animationutils_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DAnimation {
namespace Animation {

class Handler;

class Q_AUTOTEST_EXPORT BlendedClipAnimator : public BackendNode
{
public:
    BlendedClipAnimator();

    void cleanup();
    void syncFromFrontEnd(const Qt3DCore::QNode *frontEnd, bool firstTime) override;

    Qt3DCore::QNodeId blendTreeRootId() const;
    Qt3DCore::QNodeId mapperId() const { return m_mapperId; }
    Qt3DCore::QNodeId clockId() const { return m_clockId; }
    bool isRunning() const { return m_running; }
    void setNormalizedLocalTime(float normalizedTime, bool allowMarkDirty = true);
    float normalizedLocalTime() const { return m_normalizedLocalTime; }

    //  Called by BuildBlendTreeJob
    bool canRun() const { return !m_mapperId.isNull() && !m_blendTreeRootId.isNull(); }

    void setBlendTreeRootId(Qt3DCore::QNodeId blendTreeRootId);
    void setMapperId(Qt3DCore::QNodeId mapperId);
    void setClockId(Qt3DCore::QNodeId clockId);
    void setRunning(bool running);

    void setStartTime(qint64 globalTime) { m_lastGlobalTimeNS = globalTime; }

    void setLoops(int loops) { m_loops = loops; }
    int loops() const { return m_loops; }

    int currentLoop() const { return m_currentLoop; }
    void setCurrentLoop(int currentLoop) { m_currentLoop = currentLoop; }

    void setMappingData(const QVector<MappingData> &mappingData) { m_mappingData = mappingData; }
    QVector<MappingData> mappingData() const { return m_mappingData; }

    void animationClipMarkedDirty() { setDirty(Handler::BlendedClipAnimatorDirty); }

    qint64 nsSincePreviousFrame(qint64 currentGlobalTimeNS);
    void setLastGlobalTimeNS(const qint64 &lastGlobalTimeNS);

    double lastLocalTime() const;
    void setLastLocalTime(double lastLocalTime);

    float lastNormalizedLocalTime() { return m_lastNormalizedLocalTime; }
    void setLastNormalizedLocalTime(float normalizedTime);
    bool isSeeking() const
    {
        return isValidNormalizedTime(m_normalizedLocalTime)
                && !qFuzzyCompare(m_lastNormalizedLocalTime, m_normalizedLocalTime);
    }

private:
    Qt3DCore::QNodeId m_blendTreeRootId;
    Qt3DCore::QNodeId m_mapperId;
    Qt3DCore::QNodeId m_clockId;
    bool m_running;

    qint64 m_lastGlobalTimeNS;
    double m_lastLocalTime;

    int m_currentLoop;
    int m_loops;

    float m_normalizedLocalTime;
    float m_lastNormalizedLocalTime;

    QVector<MappingData> m_mappingData;
};

} // namespace Animation
} // namespace Qt3DAnimation


QT_END_NAMESPACE

#endif // QT3DANIMATION_ANIMATION_BLENDEDCLIPANIMATOR_P_H
