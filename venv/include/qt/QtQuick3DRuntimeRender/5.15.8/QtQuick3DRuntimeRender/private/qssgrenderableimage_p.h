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

#ifndef QSSG_RENDERABLE_IMAGE_H
#define QSSG_RENDERABLE_IMAGE_H

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

#include <QtQuick3DRuntimeRender/private/qtquick3druntimerenderglobal_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderimage_p.h>

QT_BEGIN_NAMESPACE

enum class QSSGImageMapTypes
{
    Unknown = 0,
    Diffuse = 1,
    Opacity = 2,
    Specular = 3,
    Emissive = 4,
    Bump = 5,
    SpecularAmountMap = 6,
    Normal = 7,
    Displacement = 8,
    Translucency = 9,
    LightmapIndirect = 10,
    LightmapRadiosity = 11,
    LightmapShadow = 12,
    Roughness = 13,
    BaseColor = 14,
    Metalness = 15,
    Occlusion = 16
};

inline uint qHash(QSSGImageMapTypes t, uint) { return qHash(static_cast<uint>(t)); }
/**
 *	Some precomputed information on a given image.  When generating a renderable, the shader
 *	generator goes through all the possible images on a material and for each valid image
 *	computes this renderable image and attaches it to the renderable.
 */
struct QSSGRenderableImage
{
    QSSGImageMapTypes m_mapType;
    QSSGRenderImage &m_image;
    QSSGRenderableImage *m_nextImage;
    QSSGRenderableImage(QSSGImageMapTypes inMapType, QSSGRenderImage &inImage)
        : m_mapType(inMapType), m_image(inImage), m_nextImage(nullptr)
    {
    }
};
QT_END_NAMESPACE
#endif
