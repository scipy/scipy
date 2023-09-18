/****************************************************************************
**
** Copyright (C) 2017 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DEXTRAS_QTEXTUREMATERIAL_P_H
#define QT3DEXTRAS_QTEXTUREMATERIAL_P_H

#include <QMatrix3x3>

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

#include <Qt3DRender/private/qmaterial_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QFilterKey;
class QEffect;
class QAbstractTexture;
class QTechnique;
class QParameter;
class QShaderProgram;
class QRenderPass;
class QNoDepthMask;
class QBlendEquationArguments;
class QBlendEquation;

} // namespace Qt3DRender

namespace Qt3DExtras {

class QTextureMaterial;

class QTextureMaterialPrivate : public Qt3DRender::QMaterialPrivate
{
    QTextureMaterialPrivate();

    void init();

    void handleTextureChanged(const QVariant &var);
    void handleTextureTransformChanged(const QVariant &var);

    Qt3DRender::QEffect *m_textureEffect;
    Qt3DRender::QParameter *m_textureParameter;
    Qt3DRender::QParameter *m_textureTransformParameter;
    Qt3DRender::QTechnique *m_textureGL3Technique;
    Qt3DRender::QTechnique *m_textureGL2Technique;
    Qt3DRender::QTechnique *m_textureES2Technique;
    Qt3DRender::QRenderPass *m_textureGL3RenderPass;
    Qt3DRender::QRenderPass *m_textureGL2RenderPass;
    Qt3DRender::QRenderPass *m_textureES2RenderPass;
    Qt3DRender::QShaderProgram *m_textureGL3Shader;
    Qt3DRender::QShaderProgram *m_textureGL2ES2Shader;
    Qt3DRender::QNoDepthMask *m_noDepthMask;
    Qt3DRender::QBlendEquationArguments *m_blendState;
    Qt3DRender::QBlendEquation *m_blendEquation;
    Qt3DRender::QFilterKey *m_filterKey;

    Q_DECLARE_PUBLIC(QTextureMaterial)
};

} // Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DEXTRAS_QTEXTUREMATERIAL_P_H

