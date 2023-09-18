/****************************************************************************
**
** Copyright (C) 2016 Klaralvdalens Datakonsult AB (KDAB).
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
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

#ifndef QT3DRENDER_RENDER_TEXTURE_H
#define QT3DRENDER_RENDER_TEXTURE_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of other Qt classes.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <Qt3DRender/private/backendnode_p.h>
#include <Qt3DRender/private/handle_types_p.h>
#include <Qt3DRender/private/qabstracttexture_p.h>
#include <Qt3DRender/qtexture.h>
#include <Qt3DRender/qtexturedata.h>
#include <Qt3DRender/qtexturegenerator.h>
#include <QOpenGLContext>
#include <QMutex>

QT_BEGIN_NAMESPACE

class QOpenGLTexture;

namespace Qt3DRender {

class QAbstractTexture;

namespace Render {

class TextureManager;
class TextureImageManager;

/**
 * General, constant properties of a texture
 */
struct TextureProperties
{
    int width = 1;
    int height = 1;
    int depth = 1;
    int layers = 1;
    int mipLevels = 1;
    int samples = 1;
    QAbstractTexture::Target target = QAbstractTexture::TargetAutomatic;
    QAbstractTexture::TextureFormat format = QAbstractTexture::NoFormat;
    bool generateMipMaps = false;
    QAbstractTexture::Status status = QAbstractTexture::None;

    bool operator==(const TextureProperties &o) const {
        return (width == o.width) && (height == o.height) && (depth == o.depth)
            && (layers == o.layers) && (mipLevels == o.mipLevels) && (target == o.target)
            && (format == o.format) && (generateMipMaps == o.generateMipMaps)
            && (samples == o.samples) && (status == o.status);
    }
    inline bool operator!=(const TextureProperties &o) const { return !(*this == o); }
};


/**
 *   Texture parameters that are independent of texture data and that may
 *   change without the re-uploading the texture
 */
struct TextureParameters
{
    QAbstractTexture::Filter magnificationFilter = QAbstractTexture::Nearest;
    QAbstractTexture::Filter minificationFilter = QAbstractTexture::Nearest;
    QTextureWrapMode::WrapMode wrapModeX = QTextureWrapMode::ClampToEdge;
    QTextureWrapMode::WrapMode wrapModeY = QTextureWrapMode::ClampToEdge;
    QTextureWrapMode::WrapMode wrapModeZ = QTextureWrapMode::ClampToEdge;
    float maximumAnisotropy = 1.0f;
    QAbstractTexture::ComparisonFunction comparisonFunction = QAbstractTexture::CompareLessEqual;
    QAbstractTexture::ComparisonMode comparisonMode = QAbstractTexture::CompareNone;

    bool operator==(const TextureParameters &o) const {
        return (magnificationFilter == o.magnificationFilter) && (minificationFilter == o.minificationFilter)
            && (wrapModeX == o.wrapModeX) && (wrapModeY == o.wrapModeY) && (wrapModeZ == o.wrapModeZ)
            && (maximumAnisotropy == o.maximumAnisotropy)
            && (comparisonFunction == o.comparisonFunction) && (comparisonMode == o.comparisonMode);
    }
    inline bool operator!=(const TextureParameters &o) const { return !(*this == o); }
};

/**
 *   Backend object for QAbstractTexture. Just holds texture properties and parameters.
 *   Will query the TextureImplManager for an instance of TextureImpl that matches it's
 *   properties.
 */
class Q_3DRENDERSHARED_PRIVATE_EXPORT Texture : public BackendNode
{
public:
    Texture();
    ~Texture();

    enum DirtyFlag {
        NotDirty = 0,
        DirtyProperties = (1 << 0),
        DirtyParameters = (1 << 1),
        DirtyImageGenerators = (1 << 2),
        DirtyDataGenerator = (1 << 3),
        DirtySharedTextureId = (1 << 4),
        DirtyPendingDataUpdates = (1 << 5),
    };
    Q_DECLARE_FLAGS(DirtyFlags, DirtyFlag)

    struct TextureUpdateInfo
    {
        TextureProperties properties;
        QVariant handle;
        QAbstractTexture::HandleType handleType;
    };

    void addDirtyFlag(DirtyFlags flags);
    DirtyFlags dirtyFlags();
    void unsetDirty();

    void cleanup();

    void addTextureDataUpdate(const QTextureDataUpdate &update);
    QVector<QTextureDataUpdate> takePendingTextureDataUpdates() { return std::move(m_pendingTextureDataUpdates); }

    void syncFromFrontEnd(const Qt3DCore::QNode *frontEnd, bool firstTime) override;

    inline const TextureProperties& properties() const { return m_properties; }
    inline const TextureParameters& parameters() const { return m_parameters; }
    inline const Qt3DCore::QNodeIdVector textureImageIds() const { return m_textureImageIds; }
    inline const QTextureGeneratorPtr& dataGenerator() const { return m_dataFunctor; }
    inline int sharedTextureId() const { return m_sharedTextureId; }

    void setDataGenerator(const QTextureGeneratorPtr &generator);
    bool isValid(TextureImageManager *manager) const;
private:
    void initializeFromPeer(const Qt3DCore::QNodeCreatedChangeBasePtr &change) final;

    DirtyFlags m_dirty;
    TextureProperties m_properties;
    TextureParameters m_parameters;
    int m_sharedTextureId;

    QTextureGeneratorPtr m_dataFunctor;
    Qt3DCore::QNodeIdVector m_textureImageIds;

    QMutex m_flagsMutex;
    QVector<QTextureDataUpdate> m_pendingTextureDataUpdates;
};

class Q_AUTOTEST_EXPORT TextureFunctor : public Qt3DCore::QBackendNodeMapper
{
public:
    explicit TextureFunctor(AbstractRenderer *renderer,
                            TextureManager *textureNodeManager);
    Qt3DCore::QBackendNode *create(const Qt3DCore::QNodeCreatedChangeBasePtr &change) const final;
    Qt3DCore::QBackendNode *get(Qt3DCore::QNodeId id) const final;
    void destroy(Qt3DCore::QNodeId id) const final;

private:
    AbstractRenderer *m_renderer;
    TextureManager *m_textureNodeManager;
};

#ifndef QT_NO_DEBUG_STREAM
inline QDebug operator<<(QDebug dbg, const Texture &texture)
{
    QDebugStateSaver saver(dbg);
    dbg << "QNodeId =" << texture.peerId() << "imageCount =" << texture.textureImageIds().size() << Qt::endl;
    return dbg;
}
#endif

} // namespace Render
} // namespace Qt3DRender

QT_END_NAMESPACE

Q_DECLARE_METATYPE(Qt3DRender::Render::Texture*) // LCOV_EXCL_LINE
Q_DECLARE_METATYPE(Qt3DRender::Render::TextureProperties) // LCOV_EXCL_LINE

#endif // QT3DRENDER_RENDER_TEXTURE_H
