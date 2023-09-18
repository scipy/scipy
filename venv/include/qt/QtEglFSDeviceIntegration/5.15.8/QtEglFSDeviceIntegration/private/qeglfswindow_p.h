/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the plugins of the Qt Toolkit.
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

#ifndef QEGLFSWINDOW_H
#define QEGLFSWINDOW_H

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

#include "qeglfsglobal_p.h"
#include "qeglfsintegration_p.h"
#include "qeglfsscreen_p.h"

#include <qpa/qplatformwindow.h>
#ifndef QT_NO_OPENGL
# include <QtPlatformCompositorSupport/private/qopenglcompositor_p.h>
#endif

QT_BEGIN_NAMESPACE

class QOpenGLCompositorBackingStore;
class QPlatformTextureList;

#ifndef QT_NO_OPENGL
class Q_EGLFS_EXPORT QEglFSWindow : public QPlatformWindow, public QOpenGLCompositorWindow
#else
class Q_EGLFS_EXPORT QEglFSWindow : public QPlatformWindow
#endif
{
public:
    QEglFSWindow(QWindow *w);
    ~QEglFSWindow();

    void create();
    void destroy();

    void setGeometry(const QRect &) override;
    QRect geometry() const override;
    void setVisible(bool visible) override;
    void requestActivateWindow() override;
    void raise() override;
    void lower() override;

    void propagateSizeHints() override { }
    void setMask(const QRegion &) override { }
    bool setKeyboardGrabEnabled(bool) override { return false; }
    bool setMouseGrabEnabled(bool) override { return false; }
    void setOpacity(qreal) override;
    WId winId() const override;

    QSurfaceFormat format() const override;

    EGLNativeWindowType eglWindow() const;
    EGLSurface surface() const;
    QEglFSScreen *screen() const override;
#if QT_CONFIG(vulkan)
    virtual void *vulkanSurfacePtr() { return nullptr; }
#endif

    bool hasNativeWindow() const { return m_flags.testFlag(HasNativeWindow); }

    void invalidateSurface() override;
    virtual void resetSurface();

#ifndef QT_NO_OPENGL
    QOpenGLCompositorBackingStore *backingStore() { return m_backingStore; }
    void setBackingStore(QOpenGLCompositorBackingStore *backingStore) { m_backingStore = backingStore; }
#endif
    bool isRaster() const;
#ifndef QT_NO_OPENGL
    QWindow *sourceWindow() const override;
    const QPlatformTextureList *textures() const override;
    void endCompositing() override;
#endif

protected:
#ifndef QT_NO_OPENGL
    QOpenGLCompositorBackingStore *m_backingStore;
    QOpenGLContext *m_rasterCompositingContext;
#endif
    WId m_winId;

    EGLSurface m_surface;
    EGLNativeWindowType m_window;

    EGLConfig m_config;
    QSurfaceFormat m_format;

    enum Flag {
        Created = 0x01,
        HasNativeWindow = 0x02
    };
    Q_DECLARE_FLAGS(Flags, Flag)
    Flags m_flags;
};

QT_END_NAMESPACE

#endif // QEGLFSWINDOW_H
