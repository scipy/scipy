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

#ifndef QRENDERSURFACESELECTOR_P_H
#define QRENDERSURFACESELECTOR_P_H

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

#include <Qt3DRender/private/qframegraphnode_p.h>
#include <Qt3DRender/private/platformsurfacefilter_p.h>
#include <Qt3DRender/private/qt3drender_global_p.h>
#include <QtGui/qsurface.h>
#include <QtGui/qwindow.h>
#include <QtCore/qpointer.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QRenderSurfaceSelector;

class Q_3DRENDERSHARED_PRIVATE_EXPORT QRenderSurfaceSelectorPrivate : public Qt3DRender::QFrameGraphNodePrivate
{
public:
    QRenderSurfaceSelectorPrivate();
    ~QRenderSurfaceSelectorPrivate();

    // TODO: Qt 5.8, make it public
    static QRenderSurfaceSelector *find(QObject *rootObject);

    void setExternalRenderTargetSize(const QSize &size);
    QSize externalRenderTargetSize() const { return m_externalRenderTargetSize; }

    QSurface *m_surface;
    QSize m_externalRenderTargetSize;
    QScopedPointer<Qt3DRender::Render::PlatformSurfaceFilter> m_surfaceEventFilter;
    float m_surfacePixelRatio;
    QMetaObject::Connection m_heightConn;
    QMetaObject::Connection m_widthConn;
    QMetaObject::Connection m_screenConn;
};

struct QRenderSurfaceSelectorData
{
    QPointer<QObject> surface;
    QSize externalRenderTargetSize;
    float surfacePixelRatio;
};

} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QRENDERSURFACESELECTOR_P_H
