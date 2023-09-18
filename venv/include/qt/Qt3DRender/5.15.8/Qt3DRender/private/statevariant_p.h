/****************************************************************************
**
** Copyright (C) 2016 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DRENDER_RENDER_STATEVARIANT_P_H
#define QT3DRENDER_RENDER_STATEVARIANT_P_H

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

#include <Qt3DRender/private/renderstates_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {
namespace Render {

struct Q_3DRENDERSHARED_PRIVATE_EXPORT StateVariant
{
    StateMask type;

    union u_Data {
        BlendEquationArguments blendEquationArguments;
        BlendEquation blendEquation;
        AlphaFunc alphaFunc;
        MSAAEnabled msaaEnabled;
        DepthRange depthRange;
        DepthTest depthTest;
        NoDepthMask noDepthMask;
        CullFace cullFace;
        FrontFace frontFace;
        Dithering dithering;
        ScissorTest scissorTest;
        StencilTest stencilTest;
        AlphaCoverage alphaCoverage;
        PointSize pointSize;
        PolygonOffset polygonOffset;
        ColorMask colorMask;
        ClipPlane clipPlane;
        SeamlessCubemap seamlessCubemap;
        StencilOp stencilOp;
        StencilMask stencilMask;
        LineWidth lineWidth;

        u_Data()
        {
            // Assumes the above types don't need to have their ctor called
            memset(static_cast<void *>(this), 0, sizeof(u_Data));
        }

        u_Data(const u_Data &other)
        {
            memcpy(static_cast<void *>(this), static_cast<const void *>(&other), sizeof(u_Data));
        }

        u_Data& operator=(const u_Data &other)
        {
            memcpy(static_cast<void *>(this), static_cast<const void *>(&other), sizeof(u_Data));
            return *this;
        }

        ~u_Data()
        {
            // Assumes the above types don't need to have their dtor called
        }
    } data;

    template<class State, typename ... Args>
    static StateVariant createState(Args... values)
    {
        State state;
        state.set(values...);
        return StateVariant::fromValue(state);
    }

    template<typename GenericState>
    static StateVariant fromValue(const GenericState &state)
    {
        StateVariant v;
        v.type = GenericState::type();
        // all union members start at the same memory address
        // so we can just write into whichever we want
        memcpy(static_cast<void *>(&v.data), static_cast<const void *>(&state), sizeof(state));
        return v;
    }

    RenderStateImpl *state();
    const RenderStateImpl *constState() const;

    bool operator ==(const StateVariant &other) const;
    bool operator !=(const StateVariant &other) const;
};

} // namespace Render
} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_RENDER_STATEVARIANT_P_H
