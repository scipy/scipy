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

#ifndef QSSGCUSTOMMATERIAL_H
#define QSSGCUSTOMMATERIAL_H

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

#include <QtQuick3D/private/qquick3dmaterial_p.h>
#include <QtCore/qvector.h>

#include <QtQuick3DRender/private/qssgrenderbasetypes_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendergraphobject_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderdynamicobjectsystemcommands_p.h>

#include <QtQuick3D/private/qquick3dshaderutils_p.h>


QT_BEGIN_NAMESPACE

class Q_QUICK3D_EXPORT QQuick3DCustomMaterial : public QQuick3DMaterial
{
    Q_OBJECT
    Q_PROPERTY(bool hasTransparency READ hasTransparency WRITE setHasTransparency NOTIFY hasTransparencyChanged)
    Q_PROPERTY(bool hasRefraction READ hasRefraction WRITE setHasRefraction NOTIFY hasRefractionChanged)
    Q_PROPERTY(bool alwaysDirty READ alwaysDirty WRITE setAlwaysDirty NOTIFY alwaysDirtyChanged)
    Q_PROPERTY(QQuick3DShaderUtilsShaderInfo *shaderInfo READ shaderInfo WRITE setShaderInfo)
    Q_PROPERTY(QQmlListProperty<QQuick3DShaderUtilsRenderPass> passes READ passes)

public:
    explicit QQuick3DCustomMaterial(QQuick3DObject *parent = nullptr);
    ~QQuick3DCustomMaterial() override;

    bool hasTransparency() const;
    bool hasRefraction() const;
    bool alwaysDirty() const;

    QQuick3DShaderUtilsShaderInfo *shaderInfo() const;
    QQmlListProperty<QQuick3DShaderUtilsRenderPass> passes();

public Q_SLOTS:
    void setHasTransparency(bool hasTransparency);
    void setHasRefraction(bool hasRefraction);
    void setShaderInfo(QQuick3DShaderUtilsShaderInfo *shaderInfo);
    void setAlwaysDirty(bool alwaysDirty);

Q_SIGNALS:
    void hasTransparencyChanged(bool hasTransparency);
    void hasRefractionChanged(bool hasRefraction);
    void alwaysDirtyChanged(bool alwaysDirty);

protected:
    QSSGRenderGraphObject *updateSpatialNode(QSSGRenderGraphObject *node) override;
    void markAllDirty() override;

private Q_SLOTS:
    void onPropertyDirty();
    void onTextureDirty(QQuick3DShaderUtilsTextureInput *texture);

private:
    enum Dirty {
        TextureDirty = 0x1,
        PropertyDirty = 0x2
    };

    // Passes
    static void qmlAppendPass(QQmlListProperty<QQuick3DShaderUtilsRenderPass> *list, QQuick3DShaderUtilsRenderPass *pass);
    static QQuick3DShaderUtilsRenderPass *qmlPassAt(QQmlListProperty<QQuick3DShaderUtilsRenderPass> *list, int index);
    static int qmlPassCount(QQmlListProperty<QQuick3DShaderUtilsRenderPass> *list);
    static void qmlPassClear(QQmlListProperty<QQuick3DShaderUtilsRenderPass> *list);

    void markDirty(QQuick3DCustomMaterial::Dirty type)
    {
        if (!(m_dirtyAttributes & quint32(type))) {
            m_dirtyAttributes |= quint32(type);
            update();
        }
    }

    quint32 m_dirtyAttributes = 0xffffffff;
    bool m_hasTransparency = false;
    bool m_hasRefraction = false;
    QQuick3DShaderUtilsShaderInfo *m_shaderInfo = nullptr;
    QVector<QQuick3DShaderUtilsRenderPass *> m_passes;
    bool m_alwaysDirty = false;
};

QT_END_NAMESPACE

#endif // QSSGCUSTOMMATERIAL_H
