/****************************************************************************
**
** Copyright (C) 2008-2012 NVIDIA Corporation.
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

#ifndef QSSG_RENDER_LIGHT_CONSTANT_PROPERTIES
#define QSSG_RENDER_LIGHT_CONSTANT_PROPERTIES

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

#include <QtQuick3DRender/private/qssgrendershaderprogram_p.h>

QT_BEGIN_NAMESPACE

struct QSSGLightConstants
{
    QSSGRenderCachedShaderProperty<QVector4D> position;
    QSSGRenderCachedShaderProperty<QVector4D> direction;
    QSSGRenderCachedShaderProperty<QVector4D> up;
    QSSGRenderCachedShaderProperty<QVector4D> right;
    QSSGRenderCachedShaderProperty<QVector4D> diffuse;
    QSSGRenderCachedShaderProperty<QVector4D> ambient;
    QSSGRenderCachedShaderProperty<QVector4D> specular;
    QSSGRenderCachedShaderProperty<float> coneAngle;
    QSSGRenderCachedShaderProperty<float> innerConeAngle;
    QSSGRenderCachedShaderProperty<float> constantAttenuation;
    QSSGRenderCachedShaderProperty<float> linearAttenuation;
    QSSGRenderCachedShaderProperty<float> quadraticAttenuation;
    QSSGRenderCachedShaderProperty<float> range;
    QSSGRenderCachedShaderProperty<float> width;
    QSSGRenderCachedShaderProperty<float> height;
    QSSGRenderCachedShaderProperty<QVector4D> shadowControls;
    QSSGRenderCachedShaderProperty<QMatrix4x4> shadowView;
    QSSGRenderCachedShaderProperty<qint32> shadowIdx;
    QSSGRenderCachedShaderProperty<QVector3D> attenuation;

    QSSGLightConstants(const QByteArray &lightRef, const QSSGRef<QSSGRenderShaderProgram> &shader);

    template<typename LightProps>
    void updateLights(LightProps &props)
    {
        position.set(props.position);
        direction.set(props.direction);
        up.set(props.up);
        right.set(props.right);
        diffuse.set(props.diffuse);
        ambient.set(props.ambient);
        specular.set(props.specular);
        coneAngle.set(props.coneAngle);
        innerConeAngle.set(props.innerConeAngle);
        constantAttenuation.set(props.constantAttenuation);
        linearAttenuation.set(props.linearAttenuation);
        quadraticAttenuation.set(props.quadraticAttenuation);
        range.set(props.range);
        width.set(props.width);
        height.set(props.height);
        shadowControls.set(props.shadowControls);
        shadowView.set(QMatrix4x4(props.shadowView));
        shadowIdx.set(props.shadowIdx);
        attenuation.set(QVector3D(props.constantAttenuation, props.linearAttenuation, props.quadraticAttenuation));
    }
};

template<typename GeneratedShader>
struct QSSGLightConstantProperties
{

    QSSGLightConstantProperties(GeneratedShader *shader, bool packed) : m_lightCount("lightCount", shader->m_shader)
    {
        m_constants.resize(shader->m_lights.size());
        for (int i = 0; i < shader->m_lights.size(); ++i) {
            QByteArray lref;
            if (packed) {
                lref += "light_";
                lref += QByteArray::number(i);
                lref += "_";
            } else {
                lref += "lights[";
                lref += QByteArray::number(i);
                lref += "].";
            }
            m_constants[i] = new QSSGLightConstants(lref, shader->m_shader);
        }
        m_lightCount.set(shader->m_lights.size());
        m_lightCountInt = shader->m_lights.size();
    }

    QSSGLightConstantProperties(const QByteArray &lseed, const QByteArray &lcount, GeneratedShader *shader, bool packed, int count)
        : m_lightCount(lcount, shader->m_shader)
    {
        m_constants.resize(count);
        for (int i = 0; i < count; ++i) {
            QByteArray lref = lseed;
            if (packed) {
                lref += "_";
                lref += QByteArray::number(i);
                lref += "_";
            } else {
                lref += "[";
                lref += QByteArray::number(i);
                lref += "].";
            }
            m_constants[i] = new QSSGLightConstants(lref, shader->m_shader);
        }
        m_lightCount.set(count);
        m_lightCountInt = count;
    }

    ~QSSGLightConstantProperties() { qDeleteAll(m_constants); }

    void updateLights(const QSSGRef<GeneratedShader> &shader)
    {
        for (int i = 0; i < m_constants.size(); ++i)
            m_constants[i]->updateLights(shader->m_lights[i].lightData);
    }
    template<typename LightProps>
    void updateLights(const QVector<QSSGRef<LightProps>> &props)
    {
        for (int i = 0; i < m_constants.size(); ++i)
            m_constants[i]->updateLights(props[i]->m_lightData);
    }

    QVector<QSSGLightConstants *> m_constants;
    QSSGRenderCachedShaderProperty<qint32> m_lightCount;
    int m_lightCountInt;
};

QT_END_NAMESPACE

#endif
