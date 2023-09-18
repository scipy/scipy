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

#ifndef QT3DRENDER_RENDER_RENDERSTATES_H
#define QT3DRENDER_RENDER_RENDERSTATES_H

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

#include <Qt3DRender/private/genericstate_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {
namespace Render {

class Q_3DRENDERSHARED_PRIVATE_EXPORT BlendEquationArguments : public GenericState<BlendEquationArguments, BlendEquationArgumentsMask, GLenum, GLenum, GLenum, GLenum, bool, int>
{
public:
    void updateProperties(const QRenderState *node) override;
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT BlendEquation : public GenericState<BlendEquation, BlendStateMask, GLenum>
{
public:
    void updateProperties(const QRenderState *node) override;
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT AlphaFunc : public GenericState<AlphaFunc, AlphaTestMask, GLenum, GLclampf>
{
public:
    void updateProperties(const QRenderState *node) override;
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT MSAAEnabled : public GenericState<MSAAEnabled, MSAAEnabledStateMask, GLboolean>
{
public:
    void updateProperties(const QRenderState *node) override;
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT DepthRange : public GenericState<DepthRange, DepthRangeMask, double, double>
{
public:
    void updateProperties(const QRenderState *node) override;
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT DepthTest : public GenericState<DepthTest, DepthTestStateMask, GLenum>
{
public:
    void updateProperties(const QRenderState *node) override;
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT RasterMode : public GenericState<RasterMode, RasterModeMask, GLenum, GLenum>
{
public:
    void updateProperties(const QRenderState *node) override;
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT NoDepthMask : public GenericState<NoDepthMask, DepthWriteStateMask, GLboolean>
{
public:
    void updateProperties(const QRenderState *node) override;
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT CullFace : public GenericState<CullFace, CullFaceStateMask, GLenum>
{
public:
    void updateProperties(const QRenderState *node) override;
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT FrontFace : public GenericState<FrontFace, FrontFaceStateMask, GLenum>
{
public:
    void updateProperties(const QRenderState *node) override;
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT Dithering : public GenericState<Dithering, DitheringStateMask>
{
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT ScissorTest : public GenericState<ScissorTest, ScissorStateMask, int, int, int, int>
{
public:
    void updateProperties(const QRenderState *node) override;
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT StencilTest : public GenericState<StencilTest, StencilTestStateMask, GLenum, int, uint, GLenum, int, uint>
{
public:
    void updateProperties(const QRenderState *node) override;
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT AlphaCoverage : public GenericState<AlphaCoverage, AlphaCoverageStateMask>
{
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT PointSize : public GenericState<PointSize, PointSizeMask, bool, GLfloat>
{
public:
    void updateProperties(const QRenderState *node) override;
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT PolygonOffset : public GenericState<PolygonOffset, PolygonOffsetStateMask, GLfloat, GLfloat>
{
public:

    void updateProperties(const QRenderState *node) override;
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT ColorMask : public GenericState<ColorMask, ColorStateMask, GLboolean, GLboolean, GLboolean, GLboolean>
{
public:
    void updateProperties(const QRenderState *node) override;
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT ClipPlane : public GenericState<ClipPlane, ClipPlaneMask, int, QVector3D, float>
{
public:
    void updateProperties(const QRenderState *node) override;
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT SeamlessCubemap : public GenericState<SeamlessCubemap, SeamlessCubemapMask>
{
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT StencilOp : public GenericState<StencilOp, StencilOpMask, GLenum, GLenum, GLenum, GLenum, GLenum, GLenum>
{
public:
    void updateProperties(const QRenderState *node) override;
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT StencilMask : public GenericState<StencilMask, StencilWriteStateMask, uint, uint>
{
public:
    void updateProperties(const QRenderState *node) override;
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT LineWidth : public GenericState<LineWidth, LineWidthMask, GLfloat, bool>
{
public:
    void updateProperties(const QRenderState *node) override;
};

} // namespace Render
} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_RENDER_RENDERSTATES_H
