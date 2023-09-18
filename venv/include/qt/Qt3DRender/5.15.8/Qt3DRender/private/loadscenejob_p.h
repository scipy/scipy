/****************************************************************************
**
** Copyright (C) 2014 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DRENDER_RENDER_LOADSCENEJOB_H
#define QT3DRENDER_RENDER_LOADSCENEJOB_H

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

#include <Qt3DCore/qaspectjob.h>
#include <Qt3DCore/private/qaspectjob_p.h>
#include <Qt3DCore/qnodeid.h>
#include <Qt3DRender/qsceneloader.h>
#include <QSharedPointer>
#include <QUrl>
#include <functional>
#include <Qt3DRender/private/qt3drender_global_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QSceneImporter;

namespace Render {

class Scene;
class NodeManagers;

class LoadSceneJob;

class Q_AUTOTEST_EXPORT LoadSceneJobPrivate : public Qt3DCore::QAspectJobPrivate
{
public:
    explicit LoadSceneJobPrivate(LoadSceneJob *q): q_ptr(q) {}
    ~LoadSceneJobPrivate() override {}

    void postFrame(Qt3DCore::QAspectManager *manager) override;

    Qt3DCore::QEntity *m_sceneSubtree = nullptr;
    QSceneLoader::Status m_status = QSceneLoader::None;

    Q_DECLARE_PUBLIC(LoadSceneJob)
private:
    LoadSceneJob *q_ptr;
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT LoadSceneJob : public Qt3DCore::QAspectJob
{
public:
    explicit LoadSceneJob(const QUrl &source, Qt3DCore::QNodeId sceneComponent);
    void setData(const QByteArray &data);
    void setNodeManagers(NodeManagers *managers) { m_managers = managers; }
    void setSceneImporters(const QList<QSceneImporter *> sceneImporters) { m_sceneImporters = sceneImporters; }

    NodeManagers *nodeManagers() const;
    QList<QSceneImporter *> sceneImporters() const;
    QUrl source() const;
    Qt3DCore::QNodeId sceneComponentId() const;

    void run() override;

private:
    QUrl m_source;
    QByteArray m_data;
    Qt3DCore::QNodeId m_sceneComponent;
    NodeManagers *m_managers;
    QList<QSceneImporter *> m_sceneImporters;

    Qt3DCore::QEntity *tryLoadScene(QSceneLoader::Status &finalStatus,
                                    const QStringList &extensions,
                                    const std::function<void (QSceneImporter *)> &importerSetupFunc);
    Q_DECLARE_PRIVATE(LoadSceneJob)
};

typedef QSharedPointer<LoadSceneJob> LoadSceneJobPtr;

} // namespace Render

} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // LOADSCENEJOB_H
