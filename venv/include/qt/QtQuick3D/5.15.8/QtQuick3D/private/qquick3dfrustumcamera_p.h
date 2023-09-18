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

#ifndef QSSGFRUSTUMCAMERA_H
#define QSSGFRUSTUMCAMERA_H

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

#include <QtQuick3D/private/qquick3dperspectivecamera_p.h>

QT_BEGIN_NAMESPACE

struct QSSGRenderCamera;
class Q_QUICK3D_EXPORT QQuick3DFrustumCamera : public QQuick3DPerspectiveCamera
{
    Q_OBJECT
    Q_PROPERTY(float top READ top WRITE setTop NOTIFY topChanged)
    Q_PROPERTY(float bottom READ bottom WRITE setBottom NOTIFY bottomChanged)
    Q_PROPERTY(float right READ right WRITE setRight NOTIFY rightChanged)
    Q_PROPERTY(float left READ right WRITE setLeft NOTIFY leftChanged)

public:
    QQuick3DFrustumCamera();

    float top() const;
    float bottom() const;
    float right() const;
    float left() const;

public Q_SLOTS:
    void setTop(float top);
    void setBottom(float bottom);
    void setRight(float right);
    void setLeft(float left);

Q_SIGNALS:
    void topChanged();
    void bottomChanged();
    void rightChanged();
    void leftChanged();

protected:
    bool checkSpatialNode(QSSGRenderCamera *camera) override;

private:
    float m_top = 0.0f;
    float m_bottom = 0.0f;
    float m_right = 0.0f;
    float m_left = 0.0f;
};

QT_END_NAMESPACE

#endif // QSSGFRUSTUMCAMERA_H
