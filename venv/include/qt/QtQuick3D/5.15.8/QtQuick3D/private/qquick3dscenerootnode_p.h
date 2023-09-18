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

#ifndef QQUICK3DSCENEROOTNODE_P_H
#define QQUICK3DSCENEROOTNODE_P_H

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

#include <QtQuick3D/private/qquick3dnode_p.h>

QT_BEGIN_NAMESPACE
class QQuick3DViewport;

class QQuick3DSceneRootNode : public QQuick3DNode
{
public:
    explicit QQuick3DSceneRootNode(QQuick3DViewport *view3D = nullptr, QQuick3DNode *parent = nullptr);
    ~QQuick3DSceneRootNode() override;

    QQuick3DViewport *view3D();

private:
    QQuick3DViewport *m_view3D = nullptr;
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QQuick3DSceneRootNode)

#endif // QQUICK3DSCENEROOTNODE_P_H
