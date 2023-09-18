/****************************************************************************
**
** Copyright (C) 2014 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DRENDER_QRENDERASPECT_H
#define QT3DRENDER_QRENDERASPECT_H

#include <Qt3DRender/qt3drender_global.h>
#include <Qt3DCore/qabstractaspect.h>

QT_BEGIN_NAMESPACE


class QOpenGLContext;

namespace Qt3DRender {

#if defined(QT_BUILD_INTERNAL)
class TestAspect;
#endif

namespace Render {
class Renderer;
class QRenderPlugin;
}

class QRenderAspectPrivate;

#if defined(QT_BUILD_INTERNAL)
class QRenderAspectTester;
#endif

class Q_3DRENDERSHARED_EXPORT QRenderAspect : public Qt3DCore::QAbstractAspect
{
    Q_OBJECT
public:
    enum RenderType {
        Synchronous,
        Threaded
    };

    explicit QRenderAspect(QObject *parent = nullptr);
    explicit QRenderAspect(RenderType type, QObject *parent = nullptr);
    ~QRenderAspect();

protected:
    QRenderAspect(QRenderAspectPrivate &dd, QObject *parent);
    Q_DECLARE_PRIVATE(QRenderAspect)

private:
    QVector<Qt3DCore::QAspectJobPtr> jobsToExecute(qint64 time) override;

    QVariant executeCommand(const QStringList &args) override;

    void onRegistered() override;
    void onUnregistered() override;

    void onEngineStartup() override;

    friend class Render::Renderer;
    friend class Render::QRenderPlugin;
#if defined(QT_BUILD_INTERNAL)
    friend class QRenderAspectTester;
    friend class TestAspect;
#endif
};

}

QT_END_NAMESPACE

#endif // QT3DRENDER_QRENDERASPECT_H
