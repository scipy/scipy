/****************************************************************************
**
** Copyright (C) 2020 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of Qt Quick 3D.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QQUICK3DEFFECT_H
#define QQUICK3DEFFECT_H

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

#include <QtQuick3D/qtquick3dglobal.h>
#include <QtQuick3D/private/qquick3dobject_p.h>
#include <QtQuick3D/private/qquick3dtexture_p.h>

#include <QtQuick3DRender/private/qssgrenderbasetypes_p.h>

#include <QtQuick3DRuntimeRender/private/qssgrenderdynamicobjectsystemcommands_p.h>

#include <QtCore/qvector.h>

#include <QtQuick3D/private/qquick3dshaderutils_p.h>

QT_BEGIN_NAMESPACE

class Q_QUICK3D_EXPORT QQuick3DEffect : public QQuick3DObject
{
    Q_OBJECT
    Q_PROPERTY(QQmlListProperty<QQuick3DShaderUtilsRenderPass> passes READ passes)
public:
    explicit QQuick3DEffect(QQuick3DObject *parent = nullptr);

    QQmlListProperty<QQuick3DShaderUtilsRenderPass> passes();

    // Passes
    static void qmlAppendPass(QQmlListProperty<QQuick3DShaderUtilsRenderPass> *list,
                              QQuick3DShaderUtilsRenderPass *pass);
    static QQuick3DShaderUtilsRenderPass *qmlPassAt(QQmlListProperty<QQuick3DShaderUtilsRenderPass> *list,
                                                    int index);
    static int qmlPassCount(QQmlListProperty<QQuick3DShaderUtilsRenderPass> *list);
    static void qmlPassClear(QQmlListProperty<QQuick3DShaderUtilsRenderPass> *list);

    void setDynamicTextureMap(QQuick3DTexture *textureMap, const QByteArray &name);

protected:
    QSSGRenderGraphObject *updateSpatialNode(QSSGRenderGraphObject *node) override;
    void itemChange(QQuick3DObject::ItemChange , const QQuick3DObject::ItemChangeData &) override;

private Q_SLOTS:
    void onPropertyDirty();
    void onTextureDirty(QQuick3DShaderUtilsTextureInput *texture);
private:
    enum Dirty {
        TextureDirty = 0x1,
        PropertyDirty = 0x2
    };

    void markDirty(QQuick3DEffect::Dirty type);

    quint32 m_dirtyAttributes = 0xffffffff;

    void updateSceneManager(const QSharedPointer<QQuick3DSceneManager> &sceneManager);

    friend class QQuick3DSceneRenderer;
    QVector<QQuick3DShaderUtilsRenderPass *> m_passes;
    QVector<QQuick3DTexture *> m_dynamicTextureMaps;
    ConnectionMap m_connections;
};

QT_END_NAMESPACE

#endif // QQUICK3DEFFECT_H
