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

#ifndef QSSG_RENDER_THREADED_IMAGE_LOADER_H
#define QSSG_RENDER_THREADED_IMAGE_LOADER_H

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

#include <QtQuick3DUtils/private/qssgdataref_p.h>
#include <QtQuick3DRender/private/qssgrenderbasetypes_p.h>

QT_BEGIN_NAMESPACE
enum class ImageLoadResult
{
    Succeeded,
    Failed,
};

class IImageLoadListener
{
public:
    QAtomicInt ref;
    virtual ~IImageLoadListener() {}
    virtual void OnImageLoadComplete(QString inPath, ImageLoadResult inResult) = 0;
    virtual void OnImageBatchComplete(quint64 inBatch) = 0;
};

typedef quint32 TImageBatchId;

class QSSGInputStreamFactory;
class QSSGBufferManager;
class QSSGAbstractThreadPool;
class QSSGPerfTimer;
class IImageBatchLoader
{
public:
    QAtomicInt ref;
    virtual ~IImageBatchLoader() {}
    // Returns an ID to the load request.  Request a block of images to be loaded.
    // Also takes an image that the buffer system will return when requested for the given
    // source paths
    // until said path is loaded.
    // An optional listener can be passed in to get callbacks about the batch.
    virtual TImageBatchId loadImageBatch(QSSGDataView<QString> inSourcePaths,
                                         QString inImageTillLoaded,
                                         IImageLoadListener *inListener,
                                         QSSGRenderContextType type) = 0;
    // Blocks if any of the images in the batch are in flight
    virtual void cancelImageBatchLoading(TImageBatchId inBatchId) = 0;
    // Blocks if the image is currently in-flight
    virtual void cancelImageLoading(QString inSourcePath) = 0;
    // Block until every image in the batch is loaded.
    virtual void blockUntilLoaded(TImageBatchId inId) = 0;

    // These are called by the render context, users don't need to call this.
    virtual void beginFrame() = 0;
    virtual void endFrame() = 0;

    static QSSGRef<IImageBatchLoader> createBatchLoader(const QSSGRef<QSSGInputStreamFactory> &inFactory,
                                                          const QSSGRef<QSSGBufferManager> &inBufferManager,
                                                          const QSSGRef<QSSGAbstractThreadPool> &inThreadPool,
                                                          QSSGPerfTimer *inTimer);
};
QT_END_NAMESPACE

#endif
