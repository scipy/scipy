/****************************************************************************
**
** Copyright (C) 2015 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DEXTRAS_QSKYBOXENTITY_P_H
#define QT3DEXTRAS_QSKYBOXENTITY_P_H

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

#include <QtGui/QVector3D>

#include <Qt3DCore/private/qentity_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QFilterKey;
class QTextureCubeMap;
class QTextureLoader;
class QShaderProgram;
class QSkyboxEntity;
class QTextureImage;
class QRenderPass;
class QTechnique;
class QParameter;
class QMaterial;
class QEffect;

} // namespace Qt3DRender

namespace Qt3DExtras {

class QCuboidMesh;

class QSkyboxEntityPrivate : public Qt3DCore::QEntityPrivate
{
    QSkyboxEntityPrivate();

    void init();
    void reloadTexture();

    Q_DECLARE_PUBLIC(QSkyboxEntity)

    Qt3DRender::QEffect *m_effect;
    Qt3DRender::QMaterial *m_material;
    Qt3DRender::QTextureCubeMap *m_skyboxTexture;
    Qt3DRender::QTextureLoader *m_loadedTexture;
    Qt3DRender::QShaderProgram *m_gl3Shader;
    Qt3DRender::QShaderProgram *m_gl2es2Shader;
    Qt3DRender::QTechnique *m_gl2Technique;
    Qt3DRender::QTechnique *m_es2Technique;
    Qt3DRender::QTechnique *m_gl3Technique;
    Qt3DRender::QFilterKey *m_filterKey;
    Qt3DRender::QRenderPass *m_gl2RenderPass;
    Qt3DRender::QRenderPass *m_es2RenderPass;
    Qt3DRender::QRenderPass *m_gl3RenderPass;
    QCuboidMesh *m_mesh;
    Qt3DRender::QParameter *m_gammaStrengthParameter;
    Qt3DRender::QParameter *m_textureParameter;
    Qt3DRender::QTextureImage *m_posXImage;
    Qt3DRender:: QTextureImage *m_posYImage;
    Qt3DRender::QTextureImage *m_posZImage;
    Qt3DRender::QTextureImage *m_negXImage;
    Qt3DRender::QTextureImage *m_negYImage;
    Qt3DRender::QTextureImage *m_negZImage;
    QString m_extension;
    QString m_baseName;
    QVector3D m_position;
    bool m_hasPendingReloadTextureCall;
};

} // Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DEXTRAS_QSKYBOXENTITY_P_H

