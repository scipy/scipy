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

#ifndef QT3DCORE_QSERVICELOCATOR_H
#define QT3DCORE_QSERVICELOCATOR_H

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

#include <Qt3DCore/qt3dcore_global.h>
#include <QtCore/qobject.h>
#include <QtCore/qscopedpointer.h>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

class QAbstractServiceProviderPrivate;

class Q_3DCORESHARED_EXPORT QAbstractServiceProvider : public QObject
{
    Q_OBJECT
public:
    virtual ~QAbstractServiceProvider();

    int type() const;
    QString description() const;

protected:
    explicit QAbstractServiceProvider(int type, const QString &description, QObject* parent = nullptr);
    explicit QAbstractServiceProvider(QAbstractServiceProviderPrivate &dd, QObject* parent = nullptr);

private:
    Q_DISABLE_COPY(QAbstractServiceProvider)
    Q_DECLARE_PRIVATE(QAbstractServiceProvider)
};

class QAbstractFrameAdvanceService;
class QOpenGLInformationService;
class QSystemInformationService;
class QServiceLocatorPrivate;
class QEventFilterService;
class QDownloadHelperService;
class QAspectEngine;

class Q_3DCORESHARED_EXPORT QServiceLocator
{
public:
    QServiceLocator(QAspectEngine *aspectEngine = nullptr);
    ~QServiceLocator();

    enum ServiceType {
        SystemInformation,
        OpenGLInformation,
        CollisionService,
        FrameAdvanceService,
        EventFilterService,
        DownloadHelperService,
#if !defined(Q_QDOC)
        DefaultServiceCount, // Add additional default services before here
#endif
        UserService = 256
    };

    void registerServiceProvider(int serviceType, QAbstractServiceProvider *provider);
    void unregisterServiceProvider(int serviceType);

    int serviceCount() const;

    template<class T>
    T *service(int serviceType)
    {
        QAbstractServiceProvider *service_ = _q_getServiceHelper(serviceType);
        return static_cast<T *>(service_);
    }

    // Convenience accessors for Qt3D provided services
    QSystemInformationService *systemInformation();
    QOpenGLInformationService *openGLInformation();
    QAbstractFrameAdvanceService *frameAdvanceService();
    QEventFilterService *eventFilterService();
    QDownloadHelperService *downloadHelperService();

private:
    Q_DISABLE_COPY(QServiceLocator)
    Q_DECLARE_PRIVATE(QServiceLocator)
    QScopedPointer<QServiceLocatorPrivate> d_ptr;

    QAbstractServiceProvider *_q_getServiceHelper(int type);
};

}

QT_END_NAMESPACE

#endif // QT3DCORE_QSERVICELOCATOR_H
