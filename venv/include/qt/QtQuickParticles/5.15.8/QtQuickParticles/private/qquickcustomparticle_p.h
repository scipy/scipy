/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtQml module of the Qt Toolkit.
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

#ifndef CUSTOM_PARTICLE_H
#define CUSTOM_PARTICLE_H

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
#include "qquickparticlepainter_p.h"
#include <private/qquickopenglshadereffectnode_p.h>
#include <private/qquickopenglshadereffect_p.h>

QT_BEGIN_NAMESPACE

class QSGNode;
struct PlainVertices;

class QQuickShaderEffectMaterialObject;

//Genealogy: Hybrid of UltraParticle and ShaderEffect
class QQuickCustomParticle : public QQuickParticlePainter
{
    Q_OBJECT
    Q_PROPERTY(QByteArray fragmentShader READ fragmentShader WRITE setFragmentShader NOTIFY fragmentShaderChanged)
    Q_PROPERTY(QByteArray vertexShader READ vertexShader WRITE setVertexShader NOTIFY vertexShaderChanged)
    QML_NAMED_ELEMENT(CustomParticle)

public:
    explicit QQuickCustomParticle(QQuickItem* parent=0);
    ~QQuickCustomParticle();

    QByteArray fragmentShader() const { return m_common.source.sourceCode[Key::FragmentShader]; }
    void setFragmentShader(const QByteArray &code);

    QByteArray vertexShader() const { return m_common.source.sourceCode[Key::VertexShader]; }
    void setVertexShader(const QByteArray &code);

Q_SIGNALS:
    void fragmentShaderChanged();
    void vertexShaderChanged();

protected:
    void initialize(int gIdx, int pIdx) override;
    void commit(int gIdx, int pIdx) override;

    QSGNode *updatePaintNode(QSGNode *, UpdatePaintNodeData *) override;
    QQuickOpenGLShaderEffectNode *prepareNextFrame(QQuickOpenGLShaderEffectNode *rootNode);
    void reset() override;
    void resize(int oldCount, int newCount);
    void componentComplete() override;
    QQuickOpenGLShaderEffectNode *buildCustomNodes();

    void sceneGraphInvalidated() override;
    void itemChange(ItemChange change, const ItemChangeData &value) override;

private Q_SLOTS:
    void sourceDestroyed(QObject *object);

private:
    void propertyChanged(int mappedId);

    typedef QQuickOpenGLShaderEffectMaterialKey Key;
    typedef QQuickOpenGLShaderEffectMaterial::UniformData UniformData;

    void buildData(QQuickOpenGLShaderEffectNode *rootNode);
    void updateVertexShader();

    QQuickOpenGLShaderEffectCommon m_common;
    const QMetaObject *m_myMetaObject;

    QHash<int, QQuickOpenGLShaderEffectNode*> m_nodes;
    qreal m_lastTime;

    uint m_dirtyUniforms : 1;
    uint m_dirtyUniformValues : 1;
    uint m_dirtyTextureProviders : 1;
    uint m_dirtyProgram : 1;
};

QT_END_NAMESPACE

#endif //HEADER_GUARD
