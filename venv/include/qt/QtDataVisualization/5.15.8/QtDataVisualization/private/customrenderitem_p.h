/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Data Visualization module of the Qt Toolkit.
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

//
//  W A R N I N G
//  -------------
//
// This file is not part of the QtDataVisualization API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#ifndef CUSTOMRENDERITEM_P_H
#define CUSTOMRENDERITEM_P_H

#include "abstractrenderitem_p.h"
#include "objecthelper_p.h"
#include <QtGui/QRgb>
#include <QtGui/QImage>
#include <QtGui/QColor>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class QCustom3DItem;
class Abstract3DRenderer;

class CustomRenderItem : public AbstractRenderItem
{
public:
    CustomRenderItem();
    virtual ~CustomRenderItem();

    inline void setTexture(GLuint texture) { m_texture = texture; }
    inline GLuint texture() const { return m_texture; }
    void setMesh(const QString &meshFile);
    inline ObjectHelper *mesh() const { return m_object; }
    inline void setScaling(const QVector3D &scaling) { m_scaling = scaling; }
    inline const QVector3D &scaling() const { return m_scaling; }
    inline void setOrigScaling(const QVector3D &scaling) { m_origScaling = scaling; }
    inline const QVector3D &origScaling() const { return m_origScaling; }
    inline void setPosition(const QVector3D &position) { m_position = position; }
    inline const QVector3D &position() const { return m_position; }
    inline void setOrigPosition(const QVector3D &position) { m_origPosition = position; }
    inline const QVector3D &origPosition() const { return m_origPosition; }
    inline void setPositionAbsolute(bool absolute) { m_positionAbsolute = absolute; }
    inline bool isPositionAbsolute() const { return m_positionAbsolute; }
    inline void setScalingAbsolute(bool absolute) { m_scalingAbsolute = absolute; }
    inline bool isScalingAbsolute() const { return m_scalingAbsolute; }
    inline void setBlendNeeded(bool blend) { m_needBlend = blend; }
    inline bool isBlendNeeded() const { return m_needBlend; }
    inline void setVisible(bool visible) { m_visible = visible; }
    inline bool isVisible() const { return m_visible; }
    inline void setItemPointer(QCustom3DItem *item) { m_item = item; }
    inline QCustom3DItem *itemPointer() const { return m_item; }
    inline void setValid(bool valid) { m_valid = valid; }
    inline bool isValid() const { return m_valid; }
    inline void setIndex(int index) { m_index = index; }
    inline int index() const { return m_index; }
    inline void setShadowCasting(bool shadowCasting) { m_shadowCasting = shadowCasting; }
    inline bool isShadowCasting() const { return m_shadowCasting; }
    inline void setFacingCamera(bool facing) { m_isFacingCamera = facing; }
    inline bool isFacingCamera() const { return m_isFacingCamera; }
    inline void setRenderer(Abstract3DRenderer *renderer) { m_renderer = renderer; }
    inline void setLabelItem(bool isLabel) { m_labelItem = isLabel; }
    inline bool isLabel() const { return m_labelItem; }

    // Volume specific
    inline void setTextureWidth(int width) { m_textureWidth = width; setSliceIndexX(m_sliceIndexX); }
    inline int textureWidth() const { return m_textureWidth; }
    inline void setTextureHeight(int height) { m_textureHeight = height; setSliceIndexY(m_sliceIndexY); }
    inline int textureHeight() const { return m_textureHeight; }
    inline void setTextureDepth(int depth) { m_textureDepth = depth; setSliceIndexZ(m_sliceIndexZ); }
    inline int textureDepth() const { return m_textureDepth; }
    inline int textureSize() const { return m_textureWidth * m_textureHeight * m_textureDepth; }
    inline void setColorTable(const QVector<QVector4D> &colors) { m_colorTable = colors; }
    void setColorTable(const QVector<QRgb> &colors);
    inline const QVector<QVector4D> &colorTable() const { return m_colorTable; }
    inline void setVolume(bool volume) { m_isVolume = volume; }
    inline bool isVolume() const { return m_isVolume; }
    inline void setTextureFormat(QImage::Format format) { m_textureFormat = format; }
    inline QImage::Format textureFormat() const { return m_textureFormat; }
    inline void setSliceIndexX(int index)
    {
        m_sliceIndexX = index;
        m_sliceFractions.setX((float(index) + 0.5f) / float(m_textureWidth) * 2.0 - 1.0);
    }
    inline void setSliceIndexY(int index)
    {
        m_sliceIndexY = index;
        m_sliceFractions.setY((float(index) + 0.5f) / float(m_textureHeight) * 2.0 - 1.0);
    }
    inline void setSliceIndexZ(int index)
    {
        m_sliceIndexZ = index;
        m_sliceFractions.setZ((float(index) + 0.5f) / float(m_textureDepth) * 2.0 - 1.0);
    }
    inline const QVector3D &sliceFractions() const { return m_sliceFractions; }
    inline int sliceIndexX() const { return m_sliceIndexX; }
    inline int sliceIndexY() const { return m_sliceIndexY; }
    inline int sliceIndexZ() const { return m_sliceIndexZ; }
    inline void setAlphaMultiplier(float mult) { m_alphaMultiplier = mult; }
    inline float alphaMultiplier() const { return m_alphaMultiplier; }
    inline void setPreserveOpacity(bool enable) { m_preserveOpacity = enable; }
    inline bool preserveOpacity() const { return m_preserveOpacity; }
    inline void setUseHighDefShader(bool enable) { m_useHighDefShader = enable; }
    inline bool useHighDefShader() const {return m_useHighDefShader; }
    void setMinBounds(const QVector3D &bounds);
    inline const QVector3D &minBounds() const { return m_minBounds; }
    void setMaxBounds(const QVector3D &bounds);
    inline const QVector3D &maxBounds() const { return m_maxBounds; }
    inline const QVector3D &minBoundsNormal() const { return m_minBoundsNormal; }
    inline const QVector3D &maxBoundsNormal() const { return m_maxBoundsNormal; }
    inline void setDrawSlices(bool enable) { m_drawSlices = enable; }
    inline bool drawSlices() const {return m_drawSlices; }
    inline void setDrawSliceFrames(bool enable) { m_drawSliceFrames = enable; }
    inline bool drawSliceFrames() const {return m_drawSliceFrames; }
    void setSliceFrameColor(const QColor &color);
    inline const QVector4D &sliceFrameColor() const { return m_sliceFrameColor; }
    inline void setSliceFrameWidths(const QVector3D &widths) { m_sliceFrameWidths = widths * 2.0f; }
    inline const QVector3D &sliceFrameWidths() const { return m_sliceFrameWidths; }
    inline void setSliceFrameGaps(const QVector3D &gaps) { m_sliceFrameGaps = gaps * 2.0f; }
    inline const QVector3D &sliceFrameGaps() const { return m_sliceFrameGaps; }
    inline void setSliceFrameThicknesses(const QVector3D &thicknesses) { m_sliceFrameThicknesses = thicknesses; }
    inline const QVector3D &sliceFrameThicknesses() const { return m_sliceFrameThicknesses; }

private:
    Q_DISABLE_COPY(CustomRenderItem)

    GLuint m_texture;
    QVector3D m_scaling;
    QVector3D m_origScaling;
    QVector3D m_position;
    QVector3D m_origPosition;
    bool m_positionAbsolute;
    bool m_scalingAbsolute;
    ObjectHelper *m_object; // shared reference
    bool m_needBlend;
    bool m_visible;
    bool m_valid;
    int m_index;
    bool m_shadowCasting;
    bool m_isFacingCamera;
    QCustom3DItem *m_item;
    Abstract3DRenderer *m_renderer;
    bool m_labelItem;

    // Volume specific
    int m_textureWidth;
    int m_textureHeight;
    int m_textureDepth;
    QVector<QVector4D> m_colorTable;
    bool m_isVolume;
    QImage::Format m_textureFormat;
    int m_sliceIndexX;
    int m_sliceIndexY;
    int m_sliceIndexZ;
    QVector3D m_sliceFractions;
    float m_alphaMultiplier;
    bool m_preserveOpacity;
    bool m_useHighDefShader;
    QVector3D m_minBounds;
    QVector3D m_maxBounds;
    QVector3D m_minBoundsNormal;
    QVector3D m_maxBoundsNormal;
    bool m_drawSlices;
    bool m_drawSliceFrames;
    QVector4D m_sliceFrameColor;
    QVector3D m_sliceFrameWidths;
    QVector3D m_sliceFrameGaps;
    QVector3D m_sliceFrameThicknesses;
};
typedef QHash<QCustom3DItem *, CustomRenderItem *> CustomRenderItemArray;

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
