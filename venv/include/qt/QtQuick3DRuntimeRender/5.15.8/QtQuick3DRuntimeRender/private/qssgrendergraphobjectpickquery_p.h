/****************************************************************************
**
** Copyright (C) 2008-2012 NVIDIA Corporation.
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

#ifndef QSSG_RENDER_GRAPH_OBJECT_PICK_QUERY_H
#define QSSG_RENDER_GRAPH_OBJECT_PICK_QUERY_H

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

#include <QtGui/QVector2D>
#include <QtGui/QMatrix4x4>

#include <QtQuick3DUtils/private/qssgdataref_p.h>

#include <QtQuick3DRender/private/qssgrenderbasetypes_p.h>

#include <QtQuick3DRuntimeRender/private/qssgrendergraphobject_p.h>

#include <limits>

QT_BEGIN_NAMESPACE

struct QSSGRenderPickResult
{
    const QSSGRenderGraphObject *m_hitObject = nullptr;
    float m_cameraDistanceSq = std::numeric_limits<float>::max();
    // The local coordinates in X,Y UV space where the hit occured
    QVector2D m_localUVCoords;
    // The position in world coordinates
    QVector3D m_scenePosition;

    QSSGRenderPickResult(const QSSGRenderGraphObject &inHitObject,
                         float inCameraDistance,
                         const QVector2D &inLocalUVCoords,
                         const QVector3D &inScenePosition);
    QSSGRenderPickResult() = default;
};

Q_STATIC_ASSERT(std::is_trivially_destructible<QSSGRenderPickResult>::value);

class QSSGGraphObjectPickQueryInterface
{
protected:
    virtual ~QSSGGraphObjectPickQueryInterface() {}

public:
    // Implementors have the option of batching the results to allow fewer virtual calls
    // or returning one item each pick.
    // Results are guaranteed to be returned nearest to furthest
    // If the return value has size of zero then we assume nothing more can be picked and the
    // pick
    // is finished.
    virtual QSSGRenderPickResult pick(const QVector2D &inMouseCoords,
                                      const QVector2D &inViewportDimensions,
                                      bool inPickEverything) = 0;
};
QT_END_NAMESPACE
#endif
