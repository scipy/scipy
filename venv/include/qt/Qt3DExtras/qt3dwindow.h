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

#ifndef QT3DWINDOW_H
#define QT3DWINDOW_H

#include <Qt3DExtras/qt3dextras_global.h>
#include <Qt3DRender/qrenderapi.h>
#include <QtGui/QWindow>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {
class QAspectEngine;
class QAbstractAspect;
class QEntity;
}

namespace Qt3DRender {
class QCamera;
class QFrameGraphNode;
class QRenderAspect;
class QRenderSettings;
}

namespace Qt3DExtras {
class QForwardRenderer;
}

namespace Qt3DInput {
class QInputAspect;
class QInputSettings;
}

namespace Qt3DLogic {
class QLogicAspect;
}

namespace Qt3DExtras {

class Qt3DWindowPrivate;

class Q_3DEXTRASSHARED_EXPORT Qt3DWindow : public QWindow
{
    Q_OBJECT
public:
    Qt3DWindow(QScreen *screen = nullptr, Qt3DRender::API = Qt3DRender::API::OpenGL);
    ~Qt3DWindow();

    void registerAspect(Qt3DCore::QAbstractAspect *aspect);
    void registerAspect(const QString &name);

    void setRootEntity(Qt3DCore::QEntity *root);

    void setActiveFrameGraph(Qt3DRender::QFrameGraphNode *activeFrameGraph);
    Qt3DRender::QFrameGraphNode *activeFrameGraph() const;
    Qt3DExtras::QForwardRenderer *defaultFrameGraph() const;

    Qt3DRender::QCamera *camera() const;
    Qt3DRender::QRenderSettings *renderSettings() const;

public Q_SLOTS:

Q_SIGNALS:

protected:
    void showEvent(QShowEvent *e) override;
    void resizeEvent(QResizeEvent *) override;
    bool event(QEvent *e) override;

private:
    Q_DECLARE_PRIVATE(Qt3DWindow)
};

Q_3DEXTRASSHARED_EXPORT
void setupWindowSurface(QWindow* window, Qt3DRender::API) noexcept;

} // Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DWINDOW_H
