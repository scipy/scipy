/****************************************************************************
**
** Copyright (C) 2020 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DRENDER_QRENDERCAPABILITIES_P_H
#define QT3DRENDER_QRENDERCAPABILITIES_P_H

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

#include <QtCore/private/qobject_p.h>
#include <Qt3DRender/qrendercapabilities.h>
#include <Qt3DRender/private/qt3drender_global_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class Q_3DRENDERSHARED_PRIVATE_EXPORT QRenderCapabilitiesPrivate : public QObjectPrivate
{
public:
    QRenderCapabilitiesPrivate();

    Q_DECLARE_PUBLIC(QRenderCapabilities)
    static const QRenderCapabilitiesPrivate *get(const QRenderCapabilities *q);

    bool m_valid;
    QRenderCapabilities::API m_api = QRenderCapabilities::OpenGL;
    QRenderCapabilities::Profile m_profile = QRenderCapabilities::NoProfile;
    int m_majorVersion = 0;
    int m_minorVersion = 0;
    QStringList m_extensions;
    QString m_vendor;
    QString m_renderer;
    QString m_version;
    QString m_glslVersion;
    int m_maxSamples = 0;
    int m_maxTextureSize = 0;
    int m_maxTextureUnits = 0;
    int m_maxTextureLayers = 0;
    bool m_supportsUBO = false;
    int m_maxUBOSize = 0;
    int m_maxUBOBindings = 0;
    bool m_supportsSSBO = false;
    int m_maxSSBOSize = 0;
    int m_maxSSBOBindings = 0;
    bool m_supportsImageStore = false;
    int m_maxImageUnits = 0;
    bool m_supportCompute = false;
    int m_maxWorkGroupCount[3] = { 0, 0, 0 };
    int m_maxWorkGroupSize[3] = { 0, 0, 0 };
    int m_maxComputeInvocations = 0;
    int m_maxComputeSharedMemorySize = 0;

    QString toString() const;
};

} // namespace Qt3Drender

QT_END_NAMESPACE

#endif // QT3DRENDER_QRENDERCAPABILITIES_P_H
