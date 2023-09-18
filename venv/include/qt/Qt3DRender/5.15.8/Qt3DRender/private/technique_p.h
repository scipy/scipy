/****************************************************************************
**
** Copyright (C) 2014 Klaralvdalens Datakonsult AB (KDAB).
** Copyright (C) 2016 The Qt Company Ltd and/or its subsidiary(-ies).
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
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

#ifndef QT3DRENDER_RENDER_TECHNIQUE_H
#define QT3DRENDER_RENDER_TECHNIQUE_H

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

#include <Qt3DRender/private/backendnode_p.h>
#include <Qt3DRender/private/parameterpack_p.h>
#include <Qt3DRender/private/qgraphicsapifilter_p.h>
#include <Qt3DRender/qfilterkey.h>
#include <QVector>
#include <QStringList>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QTechnique;
class QParameter;
class QGraphicsApiFilter;
class QFilterKey;
class QRenderPass;

namespace Render {

class TechniqueManager;

class Q_3DRENDERSHARED_PRIVATE_EXPORT Technique : public BackendNode
{
public:
    Technique();
    ~Technique();
    void cleanup();

    void syncFromFrontEnd(const Qt3DCore::QNode *frontEnd, bool firstTime) override;

    QVector<Qt3DCore::QNodeId> parameters() const;

    void appendRenderPass(Qt3DCore::QNodeId renderPassId);
    void removeRenderPass(Qt3DCore::QNodeId renderPassId);

    void appendFilterKey(Qt3DCore::QNodeId criterionId);
    void removeFilterKey(Qt3DCore::QNodeId criterionId);

    QVector<Qt3DCore::QNodeId> filterKeys() const;
    QVector<Qt3DCore::QNodeId> renderPasses() const;
    const GraphicsApiFilterData *graphicsApiFilter() const;

    bool isCompatibleWithRenderer() const;
    void setCompatibleWithRenderer(bool compatible);

    bool isCompatibleWithFilters(const Qt3DCore::QNodeIdVector &filterKeyIds);

    void setNodeManager(NodeManagers *nodeManager);
    NodeManagers *nodeManager() const;

private:

    GraphicsApiFilterData m_graphicsApiFilterData;
    ParameterPack m_parameterPack;
    QVector<Qt3DCore::QNodeId> m_filterKeyList;
    QVector<Qt3DCore::QNodeId> m_renderPasses;
    bool m_isCompatibleWithRenderer;
    NodeManagers *m_nodeManager;
};

class TechniqueFunctor : public Qt3DCore::QBackendNodeMapper
{
public:
    explicit TechniqueFunctor(AbstractRenderer *renderer, NodeManagers *manager);
    Qt3DCore::QBackendNode *create(const Qt3DCore::QNodeCreatedChangeBasePtr &change) const override;
    Qt3DCore::QBackendNode *get(Qt3DCore::QNodeId id) const override;
    void destroy(Qt3DCore::QNodeId id) const override;
private:
    NodeManagers *m_manager;
    AbstractRenderer *m_renderer;
};

} // namespace Render
} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_RENDER_TECHNIQUE_H
