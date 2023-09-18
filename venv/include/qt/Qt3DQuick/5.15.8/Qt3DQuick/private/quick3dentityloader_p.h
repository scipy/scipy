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

#ifndef QT3D_QUICK_QUICK3DENTITYLOADER_P_H
#define QT3D_QUICK_QUICK3DENTITYLOADER_P_H

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

#include <Qt3DCore/QEntity>
#include <QtCore/QObject>
#include <QtCore/QUrl>

#include <Qt3DQuick/private/qt3dquick_global_p.h>


QT_BEGIN_NAMESPACE

class QQmlComponent;

namespace Qt3DCore {

class QEntity;

namespace Quick {

class Quick3DEntityLoaderPrivate;

class Q_3DQUICKSHARED_PRIVATE_EXPORT Quick3DEntityLoader : public QEntity
{
    Q_OBJECT
    Q_PROPERTY(QObject *entity READ entity NOTIFY entityChanged)
    Q_PROPERTY(QUrl source READ source WRITE setSource NOTIFY sourceChanged)
    Q_PROPERTY(Status status READ status NOTIFY statusChanged)
    Q_PROPERTY(QQmlComponent *sourceComponent READ sourceComponent WRITE setSourceComponent NOTIFY sourceComponentChanged REVISION 12)
public:
    enum Status {
        Null = 0,
        Loading,
        Ready,
        Error
    };
    Q_ENUM(Status)

    explicit Quick3DEntityLoader(QNode *parent = 0);
    ~Quick3DEntityLoader();

    QObject *entity() const;

    QUrl source() const;
    void setSource(const QUrl &url);

    QQmlComponent *sourceComponent() const;
    void setSourceComponent(QQmlComponent *components);

    Status status() const;

Q_SIGNALS:
    void entityChanged();
    void sourceChanged();
    void sourceComponentChanged();
    void statusChanged(Status status);

private:
    Q_DECLARE_PRIVATE(Quick3DEntityLoader)
    Q_PRIVATE_SLOT(d_func(), void _q_componentStatusChanged(QQmlComponent::Status))
};

} // namespace Quick
} // namespace Qt3DCore

QT_END_NAMESPACE

#endif // QT3D_QUICK_QUICK3DENTITYLOADER_P_H
