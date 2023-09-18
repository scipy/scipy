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

#ifndef QT3DRENDER_QPICKINGSETTINGS_H
#define QT3DRENDER_QPICKINGSETTINGS_H

#include <Qt3DCore/qnode.h>
#include <Qt3DRender/qt3drender_global.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QPickingSettingsPrivate;

// TO DO: Qt 6 -> Make this a QObject

class Q_3DRENDERSHARED_EXPORT QPickingSettings : public Qt3DCore::QNode
{
    Q_OBJECT
    Q_PROPERTY(PickMethod pickMethod READ pickMethod WRITE setPickMethod NOTIFY pickMethodChanged)
    Q_PROPERTY(PickResultMode pickResultMode READ pickResultMode WRITE setPickResultMode NOTIFY pickResultModeChanged)
    Q_PROPERTY(FaceOrientationPickingMode faceOrientationPickingMode READ faceOrientationPickingMode WRITE setFaceOrientationPickingMode NOTIFY faceOrientationPickingModeChanged)
    Q_PROPERTY(float worldSpaceTolerance READ worldSpaceTolerance WRITE setWorldSpaceTolerance NOTIFY worldSpaceToleranceChanged REVISION 10)
public:
    explicit QPickingSettings(Qt3DCore::QNode *parent = nullptr);
    ~QPickingSettings();

    enum PickMethod {
        BoundingVolumePicking = 0x00,
        TrianglePicking = 0x01,
        LinePicking = 0x02,
        PointPicking = 0x04,
        PrimitivePicking = TrianglePicking | LinePicking | PointPicking
    };
    Q_ENUM(PickMethod) // LCOV_EXCL_LINE

    enum PickResultMode {
        NearestPick,
        AllPicks,
        NearestPriorityPick
    };
    Q_ENUM(PickResultMode) // LCOV_EXCL_LINE

    enum FaceOrientationPickingMode {
        FrontFace = 0x01,
        BackFace = 0x02,
        FrontAndBackFace = 0x03
    };
    Q_ENUM(FaceOrientationPickingMode) // LCOV_EXCL_LINE

    PickMethod pickMethod() const;
    PickResultMode pickResultMode() const;
    FaceOrientationPickingMode faceOrientationPickingMode() const;
    float worldSpaceTolerance() const;

public Q_SLOTS:
    void setPickMethod(PickMethod pickMethod);
    void setPickResultMode(PickResultMode pickResultMode);
    void setFaceOrientationPickingMode(FaceOrientationPickingMode faceOrientationPickingMode);
    void setWorldSpaceTolerance(float worldSpaceTolerance);

Q_SIGNALS:
    void pickMethodChanged(QPickingSettings::PickMethod pickMethod);
    void pickResultModeChanged(QPickingSettings::PickResultMode pickResult);
    void faceOrientationPickingModeChanged(QPickingSettings::FaceOrientationPickingMode faceOrientationPickingMode);
    void worldSpaceToleranceChanged(float worldSpaceTolerance);

protected:
    Q_DECLARE_PRIVATE(QPickingSettings)
    explicit QPickingSettings(QPickingSettingsPrivate &dd, Qt3DCore::QNode *parent = nullptr);
};

} // namespace Qt3Drender

QT_END_NAMESPACE

#endif // QT3DRENDER_QPICKINGSETTINGS_H
