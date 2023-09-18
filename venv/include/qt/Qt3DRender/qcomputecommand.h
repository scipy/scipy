/****************************************************************************
**
** Copyright (C) 2015 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DRENDER_QCOMPUTECOMMAND_H
#define QT3DRENDER_QCOMPUTECOMMAND_H

#include <Qt3DCore/qcomponent.h>
#include <Qt3DRender/qt3drender_global.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QComputeCommandPrivate;

class Q_3DRENDERSHARED_EXPORT QComputeCommand : public Qt3DCore::QComponent
{
    Q_OBJECT
    Q_PROPERTY(int workGroupX READ workGroupX WRITE setWorkGroupX NOTIFY workGroupXChanged)
    Q_PROPERTY(int workGroupY READ workGroupY WRITE setWorkGroupY NOTIFY workGroupYChanged)
    Q_PROPERTY(int workGroupZ READ workGroupZ WRITE setWorkGroupZ NOTIFY workGroupZChanged)
    Q_PROPERTY(RunType runType READ runType WRITE setRunType NOTIFY runTypeChanged REVISION 13)

public:
    enum RunType {
        Continuous = 0,
        Manual
    };
    Q_ENUM(RunType)

    explicit QComputeCommand(Qt3DCore::QNode *parent = nullptr);
    ~QComputeCommand();

    int workGroupX() const;
    int workGroupY() const;
    int workGroupZ() const;
    RunType runType() const;

public Q_SLOTS:
    void setWorkGroupX(int workGroupX);
    void setWorkGroupY(int workGroupY);
    void setWorkGroupZ(int workGroupZ);
    Q_REVISION(13) void setRunType(RunType runType);

    Q_REVISION(13) void trigger(int frameCount = 1);
    Q_REVISION(13) void trigger(int workGroupX, int workGroupY, int workGroupZ, int frameCount = 1);

Q_SIGNALS:
    void workGroupXChanged();
    void workGroupYChanged();
    void workGroupZChanged();
    void runTypeChanged();

private:
    Q_DECLARE_PRIVATE(QComputeCommand)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

} // Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_QCOMPUTECOMMAND_H
