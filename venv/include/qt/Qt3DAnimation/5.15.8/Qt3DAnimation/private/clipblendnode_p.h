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

#ifndef QT3DANIMATION_ANIMATION_CLIPBLENDNODE_P_H
#define QT3DANIMATION_ANIMATION_CLIPBLENDNODE_P_H

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

#include <Qt3DCore/qnodeid.h>
#include <Qt3DAnimation/private/backendnode_p.h>
#include <Qt3DAnimation/private/managers_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DAnimation {

namespace Animation {

class ClipBlendNodeManager;

class Q_AUTOTEST_EXPORT ClipBlendNode : public BackendNode
{
public:
    ~ClipBlendNode();

    enum BlendType {
        NoneBlendType,
        LerpBlendType,
        AdditiveBlendType,
        ValueType
    };

    void setClipBlendNodeManager(ClipBlendNodeManager *manager);
    inline ClipBlendNodeManager *clipBlendNodeManager() const { return m_manager; }

    BlendType blendType() const;

    void blend(Qt3DCore::QNodeId animatorId);

    void setClipResults(Qt3DCore::QNodeId animatorId, const ClipResults &clipResults);
    ClipResults clipResults(Qt3DCore::QNodeId animatorId) const;

    virtual QVector<Qt3DCore::QNodeId> allDependencyIds() const = 0;
    virtual QVector<Qt3DCore::QNodeId> currentDependencyIds() const = 0;
    virtual double duration() const = 0;

protected:
    explicit ClipBlendNode(BlendType blendType);
    void initializeFromPeer(const Qt3DCore::QNodeCreatedChangeBasePtr &change) override;
    virtual ClipResults doBlend(const QVector<ClipResults> &blendData) const = 0;

private:
    ClipBlendNodeManager *m_manager;
    BlendType m_blendType;

    // Store the results of evaluations indexed by animator id
    QVector<Qt3DCore::QNodeId> m_animatorIds;
    QVector<ClipResults> m_clipResults;
};

template<typename Backend, typename Frontend>
class ClipBlendNodeFunctor : public Qt3DCore::QBackendNodeMapper
{
public:
    explicit ClipBlendNodeFunctor(Handler *handler, ClipBlendNodeManager *manager)
        : m_handler(handler)
        , m_manager(manager)
    {
    }

    Qt3DCore::QBackendNode *create(const Qt3DCore::QNodeCreatedChangeBasePtr &change) const final
    {
        if (m_manager->containsNode(change->subjectId()))
            return static_cast<Backend *>(m_manager->lookupNode(change->subjectId()));
        Backend *backend = new Backend();
        backend->setClipBlendNodeManager(m_manager);
        backend->setHandler(m_handler);
        m_manager->appendNode(change->subjectId(), backend);
        return backend;
    }

    Qt3DCore::QBackendNode *get(Qt3DCore::QNodeId id) const final
    {
        return m_manager->lookupNode(id);
    }

    void destroy(Qt3DCore::QNodeId id) const final
    {
        m_manager->releaseNode(id);
    }

private:
    Handler *m_handler;
    ClipBlendNodeManager *m_manager;
};

} // Animation

} // Qt3DAnimation

QT_END_NAMESPACE

#endif // QT3DANIMATION_ANIMATION_CLIPBLENDNODE_P_H
