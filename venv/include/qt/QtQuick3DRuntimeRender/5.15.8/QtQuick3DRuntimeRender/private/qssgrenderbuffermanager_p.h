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

#ifndef QSSG_RENDER_BUFFER_MANAGER_H
#define QSSG_RENDER_BUFFER_MANAGER_H

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

#include <QtQuick3DRuntimeRender/private/qtquick3druntimerenderglobal_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderimagetexturedata_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendermesh_p.h>
#include <QtQuick3DAssetImport/private/qssgmeshutilities_p.h>
#include <QtQuick3DUtils/private/qssgperftimer_p.h>
#include <QtQuick3DUtils/private/qssgbounds3_p.h>

#include <QtCore/qmutex.h>

QT_BEGIN_NAMESPACE

struct QSSGRenderMesh;
struct QSSGLoadedTexture;
class QSSGRenderContext;
class QSSGInputStreamFactory;
struct QSSGMeshBVH;
namespace QSSGMeshUtilities {
    struct MultiLoadResult;
}

class Q_QUICK3DRUNTIMERENDER_EXPORT QSSGBufferManager
{
public:
    QAtomicInt ref;
private:
    typedef QSet<QString> StringSet;
    typedef QHash<QString, QSSGRenderImageTextureData> ImageMap;
    typedef QHash<QSGTexture *, QSSGRenderImageTextureData> QSGImageMap;
    typedef QHash<QSSGRenderMeshPath, QSSGRenderMesh *> MeshMap;
    typedef QHash<QString, QString> AliasImageMap;

    QSSGRef<QSSGRenderContext> context;
    QSSGRef<QSSGInputStreamFactory> inputStreamFactory;
    QSSGPerfTimer *perfTimer;
    ImageMap imageMap;
    QSGImageMap qsgImageMap;
    QMutex loadedImageSetMutex;
    StringSet loadedImageSet;
    AliasImageMap aliasImageMap;
    MeshMap meshMap;
    QVector<QSSGRenderVertexBufferEntry> entryBuffer;
    bool gpuSupportsDXT;

    void clear();

    QSSGMeshUtilities::MultiLoadResult loadPrimitive(const QString &inRelativePath) const;
    QVector<QVector3D> createPackedPositionDataArray(
            const QSSGMeshUtilities::MultiLoadResult &inResult) const;
    static void releaseMesh(QSSGRenderMesh &inMesh);
    static void releaseTexture(QSSGRenderImageTextureData &inEntry);

public:
    QSSGBufferManager(const QSSGRef<QSSGRenderContext> &inRenderContext,
                        const QSSGRef<QSSGInputStreamFactory> &inInputStreamFactory,
                        QSSGPerfTimer *inTimer);
    ~QSSGBufferManager();

    void setImageHasTransparency(const QString &inSourcePath, bool inHasTransparency);
    bool getImageHasTransparency(const QString &inSourcePath) const;
    void setImageTransparencyToFalseIfNotSet(const QString &inSourcePath);
    void setInvertImageUVCoords(const QString &inSourcePath, bool inShouldInvertCoords);

    // Returns true if this image has been loaded into memory
    // This call is threadsafe.  Nothing else on this object is guaranteed to be.
    bool isImageLoaded(const QString &inSourcePath);

    // Alias one image path with another image path.  Optionally this object will ignore the
    // call if
    // the source path is already loaded.  Aliasing is currently used to allow a default image
    // to be shown
    // in place of an image that is loading offline.
    // Returns true if the image was aliased, false otherwise.
    bool aliasImagePath(const QString &inSourcePath, const QString &inAliasPath, bool inIgnoreIfLoaded);
    void unaliasImagePath(const QString &inSourcePath);

    // Returns the given source path unless the source path is aliased; in which case returns
    // the aliased path.
    QString getImagePath(const QString &inSourcePath) const;
    // Returns a texture and a boolean indicating if this texture has transparency in it or not.
    // Can't name this LoadImage because that gets mangled by windows to LoadImageA (uggh)
    // In some cases we need to only scan particular images for transparency.
    QSSGRenderImageTextureData loadRenderImage(const QString &inImagePath,
                                                 const QSSGRef<QSSGLoadedTexture> &inTexture,
                                                 bool inForceScanForTransparency = false,
                                                 bool inBsdfMipmaps = false);
    QSSGRenderImageTextureData loadRenderImage(const QString &inSourcePath,
                                                 const QSSGRenderTextureFormat &inFormat = QSSGRenderTextureFormat::Unknown,
                                                 bool inForceScanForTransparency = false,
                                                 bool inBsdfMipmaps = false);
    QSSGRenderImageTextureData loadRenderImage(QSGTexture *qsgTexture);
    QSSGRenderMesh *getMesh(const QSSGRenderMeshPath &inSourcePath) const;
    QSSGRenderMesh *loadMesh(const QSSGRenderMeshPath &inSourcePath);
    QSSGRenderMesh *loadCustomMesh(const QSSGRenderMeshPath &inSourcePath,
                                   QSSGMeshUtilities::Mesh *mesh,
                                   bool update = false);
    QSSGMeshBVH *loadMeshBVH(const QSSGRenderMeshPath &inSourcePath);
    QSSGMeshUtilities::MultiLoadResult loadMeshData(const QSSGRenderMeshPath &inSourcePath) const;

    QSSGRenderMesh *createMesh(const QString &inSourcePath,
                                         quint8 *inVertData,
                                         quint32 inNumVerts,
                                         quint32 inVertStride,
                                         quint32 *inIndexData,
                                         quint32 inIndexCount,
                                         QSSGBounds3 inBounds);
    QSSGRenderMesh *createRenderMesh(const QSSGMeshUtilities::MultiLoadResult &result,
                                     const QSSGRenderMeshPath &inSourcePath);

    // Remove *all* buffers from the buffer manager;

    void invalidateBuffer(const QString &inSourcePath);

};
QT_END_NAMESPACE

#endif
