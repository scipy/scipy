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

#ifndef QEGLFSINTEGRATION_H
#define QEGLFSINTEGRATION_H

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
#include <QtCore/QPointer>
#include <QtCore/QVariant>
#include <QtGui/QWindow>
#include <qpa/qplatformintegration.h>
#include <qpa/qplatformnativeinterface.h>
#include <qpa/qplatformscreen.h>

QT_BEGIN_NAMESPACE

class QEglFSWindow;
class QEglFSContext;
class QFbVtHandler;
class QEvdevKeyboardManager;

class Q_EGLFS_EXPORT QEglFSIntegration : public QPlatformIntegration, public QPlatformNativeInterface
{
public:
    QEglFSIntegration();

    void initialize() override;
    void destroy() override;

    EGLDisplay display() const { return m_display; }

    QAbstractEventDispatcher *createEventDispatcher() const override;
    QPlatformFontDatabase *fontDatabase() const override;
    QPlatformServices *services() const override;
    QPlatformInputContext *inputContext() const override { return m_inputContext; }
    QPlatformTheme *createPlatformTheme(const QString &name) const override;

    QPlatformWindow *createPlatformWindow(QWindow *window) const override;
    QPlatformBackingStore *createPlatformBackingStore(QWindow *window) const override;
#ifndef QT_NO_OPENGL
    QPlatformOpenGLContext *createPlatformOpenGLContext(QOpenGLContext *context) const override;
    QPlatformOffscreenSurface *createPlatformOffscreenSurface(QOffscreenSurface *surface) const override;
#endif
#if QT_CONFIG(vulkan)
    QPlatformVulkanInstance *createPlatformVulkanInstance(QVulkanInstance *instance) const override;
#endif
    bool hasCapability(QPlatformIntegration::Capability cap) const override;

    QPlatformNativeInterface *nativeInterface() const override;

    // QPlatformNativeInterface
    void *nativeResourceForIntegration(const QByteArray &resource) override;
    void *nativeResourceForScreen(const QByteArray &resource, QScreen *screen) override;
    void *nativeResourceForWindow(const QByteArray &resource, QWindow *window) override;
#ifndef QT_NO_OPENGL
    void *nativeResourceForContext(const QByteArray &resource, QOpenGLContext *context) override;
#endif
    NativeResourceForContextFunction nativeResourceFunctionForContext(const QByteArray &resource) override;

    QFunctionPointer platformFunction(const QByteArray &function) const override;

    QFbVtHandler *vtHandler() { return m_vtHandler.data(); }

    QPointer<QWindow> pointerWindow() { return m_pointerWindow; }
    void setPointerWindow(QWindow *pointerWindow) { m_pointerWindow = pointerWindow; }

private:
    EGLNativeDisplayType nativeDisplay() const;
    void createInputHandlers();
    static void loadKeymapStatic(const QString &filename);
    static void switchLangStatic();

    EGLDisplay m_display;
    QPlatformInputContext *m_inputContext;
    QScopedPointer<QPlatformFontDatabase> m_fontDb;
    QScopedPointer<QPlatformServices> m_services;
    QScopedPointer<QFbVtHandler> m_vtHandler;
    QEvdevKeyboardManager *m_kbdMgr;
    QPointer<QWindow> m_pointerWindow;
    bool m_disableInputHandlers;
};

QT_END_NAMESPACE

#endif // QEGLFSINTEGRATION_H
