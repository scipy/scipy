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

#ifndef QT3DRENDER_RENDER_GENERICSTATE_H
#define QT3DRENDER_RENDER_GENERICSTATE_H

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
#include <Qt3DRender/private/qrenderstatecreatedchange_p.h>
#include <Qt3DCore/private/qresourcemanager_p.h>
#include <Qt3DRender/private/statemask_p.h>
#include <Qt3DRender/private/qt3drender_global_p.h>
//#include <Qt3DRender/private/statevariant_p.h>
#include <QList>
#include <QVector3D>
#include <QOpenGLContext>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QRenderState;

namespace Render {

class GraphicsContext;

typedef quint64 StateMaskSet;

class Q_3DRENDERSHARED_PRIVATE_EXPORT RenderStateImpl
{
public:
    virtual ~RenderStateImpl() {}

    virtual StateMask mask() const = 0;
    virtual bool equalTo(const RenderStateImpl &renderState) const = 0;
    virtual void updateProperties(const QRenderState *);
};

template <class StateSetImpl, StateMask stateMask, typename ... T>
class GenericState : public RenderStateImpl
{
public:
    GenericState *set(T... values)
    {
        m_values = std::tuple<T ...>(values...);
        return this;
    }

    bool equalTo(const RenderStateImpl &renderState) const override
    {
        const GenericState *other = static_cast<const GenericState*>(&renderState);
        return (other != nullptr && other->m_values == m_values);
    }

    StateMask mask() const override
    {
        return GenericState::type();
    }

    static StateMask type()
    {
        return stateMask;
    }

    std::tuple<T ...> values() const
    {
        return m_values;
    }

protected:
    std::tuple<T ...> m_values;
};

} // namespace Render
} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_RENDER_GENERICSTATE_H
