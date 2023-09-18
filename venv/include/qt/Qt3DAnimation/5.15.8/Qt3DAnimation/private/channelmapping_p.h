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

#ifndef QT3DANIMATION_ANIMATION_CHANNELMAPPING_P_H
#define QT3DANIMATION_ANIMATION_CHANNELMAPPING_P_H

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
#include <Qt3DAnimation/private/fcurve_p.h>
#include <Qt3DAnimation/qanimationcallback.h>

#include <Qt3DCore/qnodeid.h>
#include <QtCore/QMetaProperty>

QT_BEGIN_NAMESPACE

namespace Qt3DAnimation {

namespace Animation {

class Handler;

class Q_AUTOTEST_EXPORT ChannelMapping : public BackendNode
{
public:
    enum MappingType {
        ChannelMappingType = 0,
        SkeletonMappingType,
        CallbackMappingType
    };

    ChannelMapping();

    void cleanup();

    void syncFromFrontEnd(const Qt3DCore::QNode *frontEnd, bool firstTime) override;

    void setChannelName(const QString &channelName) { m_channelName = channelName; }
    QString channelName() const { return m_channelName; }

    void setTargetId(Qt3DCore::QNodeId targetId) { m_targetId = targetId; }
    Qt3DCore::QNodeId targetId() const { return m_targetId; }

    void setType(int type) { m_type = type; }
    int type() const { return m_type; }

    void setPropertyName(const char *propertyName) { m_propertyName = propertyName; }
    const char *propertyName() const { return m_propertyName; }

    void setComponentCount(int componentCount) { m_componentCount = componentCount; }
    int componentCount() const { return m_componentCount; }

    void setCallback(QAnimationCallback *callback) { m_callback = callback; }
    QAnimationCallback *callback() const { return m_callback; }

    void setCallbackFlags(QAnimationCallback::Flags flags) { m_callbackFlags = flags; }
    QAnimationCallback::Flags callbackFlags() const { return m_callbackFlags; }

    void setSkeletonId(Qt3DCore::QNodeId skeletonId) { m_skeletonId = skeletonId; }
    Qt3DCore::QNodeId skeletonId() const { return m_skeletonId; }
    Skeleton *skeleton() const;

    void setMappingType(MappingType mappingType) { m_mappingType = mappingType; }
    MappingType mappingType() const { return m_mappingType; }

private:
    // Properties from QChannelMapping
    QString m_channelName;
    Qt3DCore::QNodeId m_targetId;
    int m_type;
    int m_componentCount;
    const char *m_propertyName;

    // TODO: Properties from QCallbackMapping
    QAnimationCallback *m_callback;
    QAnimationCallback::Flags m_callbackFlags;

    // Properties from QSkeletonMapping
    Qt3DCore::QNodeId m_skeletonId;

    MappingType m_mappingType;
};

} // namespace Animation
} // namespace Qt3DAnimation


QT_END_NAMESPACE

#endif // QT3DANIMATION_ANIMATION_CHANNELMAPPING_P_H
