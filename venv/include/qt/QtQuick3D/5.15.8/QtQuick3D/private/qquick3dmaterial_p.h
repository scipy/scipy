/****************************************************************************
**
** Copyright (C) 2019 The Qt Company Ltd.
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

#ifndef QSSGMATERIAL_H
#define QSSGMATERIAL_H

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

#include <QtQuick3D/qquick3dobject.h>
#include <QtQuick3D/private/qquick3dtexture_p.h>

#include <QtCore/QVector>

QT_BEGIN_NAMESPACE

class QQuick3DSceneManager;
class Q_QUICK3D_EXPORT QQuick3DMaterial : public QQuick3DObject
{
    Q_OBJECT
    Q_PROPERTY(QQuick3DTexture *lightmapIndirect READ lightmapIndirect WRITE setLightmapIndirect NOTIFY lightmapIndirectChanged)
    Q_PROPERTY(QQuick3DTexture *lightmapRadiosity READ lightmapRadiosity WRITE setLightmapRadiosity NOTIFY lightmapRadiosityChanged)
    Q_PROPERTY(QQuick3DTexture *lightmapShadow READ lightmapShadow WRITE setLightmapShadow NOTIFY lightmapShadowChanged)
    Q_PROPERTY(QQuick3DTexture *lightProbe READ lightProbe WRITE setLightProbe NOTIFY lightProbeChanged)

    Q_PROPERTY(QQuick3DTexture *displacementMap READ displacementMap WRITE setDisplacementMap NOTIFY displacementMapChanged)
    Q_PROPERTY(float displacementAmount READ displacementAmount WRITE setDisplacementAmount NOTIFY displacementAmountChanged)
    Q_PROPERTY(CullMode cullMode READ cullMode WRITE setCullMode NOTIFY cullModeChanged)

public:
    enum CullMode {
        BackFaceCulling = 1,
        FrontFaceCulling = 2,
        NoCulling = 3,
    };
    Q_ENUM(CullMode)

    enum TextureChannelMapping {
        R = 0,
        G,
        B,
        A,
    };
    Q_ENUM(TextureChannelMapping)

    ~QQuick3DMaterial() override;

    QQuick3DTexture *lightmapIndirect() const;
    QQuick3DTexture *lightmapRadiosity() const;
    QQuick3DTexture *lightmapShadow() const;
    QQuick3DTexture *lightProbe() const;

    QQuick3DTexture *displacementMap() const;
    float displacementAmount() const;
    CullMode cullMode() const;

public Q_SLOTS:
    void setLightmapIndirect(QQuick3DTexture *lightmapIndirect);
    void setLightmapRadiosity(QQuick3DTexture *lightmapRadiosity);
    void setLightmapShadow(QQuick3DTexture *lightmapShadow);
    void setLightProbe(QQuick3DTexture *lightProbe);

    void setDisplacementMap(QQuick3DTexture *displacementMap);
    void setDisplacementAmount(float displacementAmount);
    void setCullMode(CullMode cullMode);

Q_SIGNALS:
    void lightmapIndirectChanged(QQuick3DTexture *lightmapIndirect);
    void lightmapRadiosityChanged(QQuick3DTexture *lightmapRadiosity);
    void lightmapShadowChanged(QQuick3DTexture *lightmapShadow);
    void lightProbeChanged(QQuick3DTexture *lightProbe);

    void displacementMapChanged(QQuick3DTexture *displacementMap);
    void displacementAmountChanged(float displacementAmount);
    void cullModeChanged(CullMode cullMode);

protected:
    explicit QQuick3DMaterial(QQuick3DObjectPrivate &dd, QQuick3DObject *parent = nullptr);
    QSSGRenderGraphObject *updateSpatialNode(QSSGRenderGraphObject *node) override;
    void itemChange(ItemChange, const ItemChangeData &) override;
public:
    void setDynamicTextureMap(QQuick3DTexture *textureMap, const QByteArray &name);
private:
    void updateSceneManager(const QSharedPointer<QQuick3DSceneManager> &sceneManager);
    QQuick3DTexture *m_lightmapIndirect = nullptr;
    QQuick3DTexture *m_lightmapRadiosity = nullptr;
    QQuick3DTexture *m_lightmapShadow = nullptr;
    QQuick3DTexture *m_iblProbe = nullptr;

    QQuick3DTexture *m_displacementMap = nullptr;
    float m_displacementAmount = 0.0f;
    CullMode m_cullMode = CullMode::BackFaceCulling;

    QHash<QByteArray, QMetaObject::Connection> m_connections;
    QVector<QQuick3DTexture *> m_dynamicTextureMaps;
};

QT_END_NAMESPACE

#endif // QSSGMATERIAL_H
