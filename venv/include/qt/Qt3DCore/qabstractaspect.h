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

#ifndef QT3DCORE_QABSTRACTASPECT_H
#define QT3DCORE_QABSTRACTASPECT_H

#include <Qt3DCore/qt3dcore_global.h>
#include <Qt3DCore/qnodeid.h>
#include <Qt3DCore/qscenechange.h>
#include <QtCore/QObject>
#include <QtCore/QSharedPointer>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

class QAspectEngine;
class QAspectJob;
class QAspectManager;
class QNode;
class QEntity;
class QAbstractAspectPrivate;
class QBackendNodeMapper;

typedef QSharedPointer<QAspectJob> QAspectJobPtr;
typedef QSharedPointer<QBackendNodeMapper> QBackendNodeMapperPtr;

class Q_3DCORESHARED_EXPORT QAbstractAspect : public QObject
{
    Q_OBJECT

public:
    explicit QAbstractAspect(QObject *parent = nullptr);
    ~QAbstractAspect();

    void scheduleSingleShotJob(const Qt3DCore::QAspectJobPtr &job);

protected:
    explicit QAbstractAspect(QAbstractAspectPrivate &dd, QObject *parent = nullptr);

    QNodeId rootEntityId() const Q_DECL_NOEXCEPT;

    template<class Frontend>
    void registerBackendType(const QBackendNodeMapperPtr &functor);
    template<class Frontend, bool supportsSyncing>
    void registerBackendType(const QBackendNodeMapperPtr &functor);
    void registerBackendType(const QMetaObject &obj, const QBackendNodeMapperPtr &functor);
    void registerBackendType(const QMetaObject &obj, const QBackendNodeMapperPtr &functor, bool supportsSyncing);
    template<class Frontend>
    void unregisterBackendType();
    void unregisterBackendType(const QMetaObject &);

private:
    virtual QVariant executeCommand(const QStringList &args);

    virtual QVector<QAspectJobPtr> jobsToExecute(qint64 time);

    virtual void onRegistered();
    virtual void onUnregistered();

    virtual void onEngineStartup();
    virtual void onEngineShutdown();

    Q_DECLARE_PRIVATE(QAbstractAspect)
    friend class QAspectEngine;
    friend class QAspectManager;
};

template<class Frontend>
void QAbstractAspect::registerBackendType(const QBackendNodeMapperPtr &functor)
{
    registerBackendType(Frontend::staticMetaObject, functor);
}

template<class Frontend, bool supportsSyncing>
void QAbstractAspect::registerBackendType(const QBackendNodeMapperPtr &functor)
{
    registerBackendType(Frontend::staticMetaObject, functor, supportsSyncing);
}

template<class Frontend>
void QAbstractAspect::unregisterBackendType()
{
    unregisterBackendType(Frontend::staticMetaObject);
}

} // namespace Qt3DCore

QT_END_NAMESPACE

#define QT3D_REGISTER_NAMESPACED_ASPECT(name, AspectNamespace, AspectType) \
    QT_BEGIN_NAMESPACE \
    namespace Qt3DCore { \
        typedef QAbstractAspect *(*AspectCreateFunction)(QObject *); \
        QT_DEPRECATED Q_3DCORESHARED_EXPORT void qt3d_QAspectFactory_addDefaultFactory(const QString &, const QMetaObject *, AspectCreateFunction); \
        Q_3DCORESHARED_EXPORT void qt3d_QAspectFactory_addDefaultFactory(const QLatin1String &, const QMetaObject *, AspectCreateFunction); \
    } \
    QT_END_NAMESPACE \
    namespace { \
    Qt3DCore::QAbstractAspect *qt3d_ ## AspectType ## _createFunction(QObject *parent) \
    { \
        using namespace AspectNamespace; \
        return new AspectType(parent); \
    } \
    \
    void qt3d_ ## AspectType ## _registerFunction() \
    { \
        using namespace AspectNamespace; \
        qt3d_QAspectFactory_addDefaultFactory(QLatin1String(name), &AspectType::staticMetaObject, qt3d_ ## AspectType ## _createFunction); \
    } \
    \
    Q_CONSTRUCTOR_FUNCTION(qt3d_ ## AspectType ## _registerFunction) \
    }

#define QT3D_REGISTER_ASPECT(name, AspectType) \
    QT3D_REGISTER_NAMESPACED_ASPECT(name, QT_PREPEND_NAMESPACE(Qt3DCore), AspectType)

#endif // QT3DCORE_ABSTRACTASPECT_H
