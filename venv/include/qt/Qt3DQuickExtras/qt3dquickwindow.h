/****************************************************************************
**
** Copyright (C) 2016 Klaralvdalens Datakonsult AB (KDAB).
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** BSD License Usage
** Alternatively, you may use this file under the terms of the BSD license
** as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of The Qt Company Ltd nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QT3DEXTRAS_QUICK_QT3DQUICKWINDOW_H
#define QT3DEXTRAS_QUICK_QT3DQUICKWINDOW_H

#include <Qt3DQuickExtras/qt3dquickextras_global.h>
#include <QtCore/QUrl>
#include <QtCore/qpointer.h>
#include <QtQuick/QQuickWindow>

QT_BEGIN_NAMESPACE

class QQmlIncubationController;

namespace Qt3DCore {
class QAbstractAspect;
namespace Quick {
class QQmlAspectEngine;
}
}

namespace Qt3DRender {
class QRenderAspect;
class QCamera;
}

namespace Qt3DInput {
class QInputAspect;
}

namespace Qt3DLogic {
class QLogicAspect;
}

namespace Qt3DExtras {

namespace Quick {

class Qt3DQuickWindowPrivate;

class Q_3DQUICKEXTRASSHARED_EXPORT Qt3DQuickWindow : public QWindow
{
    Q_OBJECT
    Q_PROPERTY(CameraAspectRatioMode cameraAspectRatioMode READ cameraAspectRatioMode WRITE setCameraAspectRatioMode NOTIFY cameraAspectRatioModeChanged)

public:
    explicit Qt3DQuickWindow(QWindow *parent = nullptr);
    ~Qt3DQuickWindow();

    void registerAspect(Qt3DCore::QAbstractAspect *aspect);
    void registerAspect(const QString &name);

    void setSource(const QUrl &source);
    Qt3DCore::Quick::QQmlAspectEngine *engine() const;

    enum CameraAspectRatioMode {
        AutomaticAspectRatio,
        UserAspectRatio
    };
    Q_ENUM(CameraAspectRatioMode); // LCOV_EXCL_LINE

    void setCameraAspectRatioMode(CameraAspectRatioMode mode);
    CameraAspectRatioMode cameraAspectRatioMode() const;

Q_SIGNALS:
    void cameraAspectRatioModeChanged(CameraAspectRatioMode mode);

protected:
    void showEvent(QShowEvent *e) override;

private:
    void onSceneCreated(QObject *rootObject);
    void setWindowSurface(QObject *rootObject);
    void setCameraAspectModeHelper();
    void updateCameraAspectRatio();

    Q_DECLARE_PRIVATE(Qt3DQuickWindow)
};

} // Quick

} // Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DEXTRAS_QUICK_QT3DQUICKWINDOW_H
