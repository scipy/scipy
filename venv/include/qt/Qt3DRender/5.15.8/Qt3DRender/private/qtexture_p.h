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

#ifndef QT3DRENDER_QTEXTURE_P_H
#define QT3DRENDER_QTEXTURE_P_H

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

#include <Qt3DCore/QNodeId>
#include <Qt3DCore/private/qdownloadhelperservice_p.h>
#include <Qt3DRender/private/qabstracttexture_p.h>
#include <Qt3DRender/qtexturegenerator.h>
#include <Qt3DRender/qtexture.h>
#include <Qt3DRender/private/qt3drender_global_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class Q_3DRENDERSHARED_PRIVATE_EXPORT QTextureLoaderPrivate : public QAbstractTexturePrivate
{
public:
    QTextureLoaderPrivate();

    Q_DECLARE_PUBLIC(QTextureLoader)

    void setScene(Qt3DCore::QScene *scene) override;
    void updateGenerator();

    QUrl m_source;
    bool m_mirrored;
};

class QTextureFromSourceGenerator;
typedef QSharedPointer<QTextureFromSourceGenerator> QTextureFromSourceGeneratorPtr;

class Q_AUTOTEST_EXPORT TextureDownloadRequest : public Qt3DCore::QDownloadRequest
{
public:
    TextureDownloadRequest(const QTextureFromSourceGeneratorPtr &functor,
                           const QUrl &url,
                           Qt3DCore::QAspectEngine *engine,
                           Qt3DCore::QNodeId texNodeId);

    void onCompleted() override;

private:
    QTextureFromSourceGeneratorPtr m_functor;
    Qt3DCore::QAspectEngine *m_engine;
    Qt3DCore::QNodeId m_texNodeId;
};

class Q_AUTOTEST_EXPORT QTextureFromSourceGenerator : public QTextureGenerator,
                                                      public QEnableSharedFromThis<QTextureFromSourceGenerator>
{
public:
    explicit QTextureFromSourceGenerator(QTextureLoader *textureLoader,
                                         Qt3DCore::QAspectEngine *engine,
                                         Qt3DCore::QNodeId textureId);

    QTextureFromSourceGenerator(const QTextureFromSourceGenerator &other);

    QTextureDataPtr operator ()() override;
    bool operator ==(const QTextureGenerator &other) const override;
    inline QAbstractTexture::Status status() const { return m_status; }

    QT3D_FUNCTOR(QTextureFromSourceGenerator)

    QUrl url() const;
    bool isMirrored() const;

private:
    friend class TextureDownloadRequest;

    QUrl m_url;
    QAbstractTexture::Status m_status;
    bool m_mirrored;

    QByteArray m_sourceData;
    Qt3DCore::QNodeId m_texture;
    Qt3DCore::QAspectEngine *m_engine;

    // Options that can be overridden on TextureLoader when loading
    QAbstractTexture::TextureFormat m_format;
};

class Q_AUTOTEST_EXPORT TextureLoadingHelper
{
public:
    static QTextureImageDataPtr loadTextureData(const QUrl &source, bool allow3D, bool mirrored);
    static QTextureImageDataPtr loadTextureData(QIODevice *data, const QString& suffix,
                                                bool allow3D, bool mirrored);
};

} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_QTEXTURE_P_H
