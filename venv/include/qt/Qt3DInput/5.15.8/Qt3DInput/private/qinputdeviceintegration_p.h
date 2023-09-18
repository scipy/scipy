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

#ifndef QT3DINPUT_QINPUTDEVICEINTEGRATION_P_H
#define QT3DINPUT_QINPUTDEVICEINTEGRATION_P_H

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

#include <Qt3DInput/qt3dinput_global.h>
#include <Qt3DCore/qaspectjob.h>
#include <Qt3DCore/qnodeid.h>
#include <QtCore/QObject>

#include <Qt3DInput/private/qabstractphysicaldevicebackendnode_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {
class QBackendNodeMapper;
typedef QSharedPointer<QBackendNodeMapper> QBackendNodeMapperPtr;
}


namespace Qt3DInput {

class QInputAspect;
class QAbstractPhysicalDevice;
class QInputDeviceIntegrationPrivate;

class Q_3DINPUTSHARED_PRIVATE_EXPORT QInputDeviceIntegration : public QObject
{
    Q_OBJECT
protected:
    explicit QInputDeviceIntegration(QObject *parent = nullptr);
    explicit QInputDeviceIntegration(QInputDeviceIntegrationPrivate &dd, QObject *parent = nullptr);

    template<class Frontend>
    void registerBackendType(const Qt3DCore::QBackendNodeMapperPtr &functor)
    {
        registerBackendType(Frontend::staticMetaObject, functor);
    }

    void registerBackendType(const QMetaObject &metaObject, const Qt3DCore::QBackendNodeMapperPtr &functor);

public:
    void initialize(Qt3DInput::QInputAspect *aspect);

    virtual QVector<Qt3DCore::QAspectJobPtr> jobsToExecute(qint64 time) = 0;
    virtual QAbstractPhysicalDevice *createPhysicalDevice(const QString &name) = 0;
    virtual QVector<Qt3DCore::QNodeId> physicalDevices() const = 0;
    virtual QAbstractPhysicalDeviceBackendNode *physicalDevice(Qt3DCore::QNodeId id) const = 0;
    virtual QStringList deviceNames() const = 0;

protected:
    QInputAspect *inputAspect() const;

private:
    virtual void onInitialize() = 0;

    Q_DECLARE_PRIVATE(QInputDeviceIntegration)
};

} // namespace Qt3DInput

QT_END_NAMESPACE

#endif // QT3DINPUT_QINPUTDEVICEINTEGRATION_P_H
