/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL3$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see http://www.qt.io/terms-conditions. For further
** information use the contact form at http://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPLv3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or later as published by the Free
** Software Foundation and appearing in the file LICENSE.GPL included in
** the packaging of this file. Please review the following information to
** ensure the GNU General Public License version 2.0 requirements will be
** met: http://www.gnu.org/licenses/gpl-2.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QT3DSCENE2D_QUICK_QUICK3DSCENE2D_P_H
#define QT3DSCENE2D_QUICK_QUICK3DSCENE2D_P_H

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

#include <QtQml/qqmllist.h>

#include <Qt3DQuickScene2D/qscene2d.h>

#include <private/qt3dquickscene2d_global_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {
namespace Render {
namespace Quick {

class Q_3DQUICKSCENE2DSHARED_EXPORT QQuick3DScene2D : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QQmlListProperty<Qt3DCore::QEntity> entities READ entities)

public:
    QQuick3DScene2D(QObject *parent = nullptr);

    inline Qt3DRender::Quick::QScene2D *parentScene2D() const
    {
        return qobject_cast<Qt3DRender::Quick::QScene2D *>(parent());
    }

    QQmlListProperty<Qt3DCore::QEntity> entities();

private:

    static void appendEntity(QQmlListProperty<Qt3DCore::QEntity> *list, Qt3DCore::QEntity *entity);
    static Qt3DCore::QEntity *entityAt(QQmlListProperty<Qt3DCore::QEntity> *list, int index);
    static int entityCount(QQmlListProperty<Qt3DCore::QEntity> *list);
    static void clearEntities(QQmlListProperty<Qt3DCore::QEntity> *list);
};

} // namespace Quick
} // namespace Render
} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DANIMATION_QUICK_QUICK3DMORPHINGANIMATION_P_H

