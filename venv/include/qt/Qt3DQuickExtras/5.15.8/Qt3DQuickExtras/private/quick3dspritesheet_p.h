/****************************************************************************
**
** Copyright (C) 2017 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DEXTRAS_EXTRAS_QUICK_QUICK3DSPRITESHEET_P_H
#define QT3DEXTRAS_EXTRAS_QUICK_QUICK3DSPRITESHEET_P_H

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

#include <Qt3DQuickExtras/qt3dquickextras_global.h>
#include <Qt3DExtras/qspritesheet.h>
#include <QtQml/QQmlListProperty>

QT_BEGIN_NAMESPACE

namespace Qt3DExtras {
namespace Extras {
namespace Quick {

class Q_3DQUICKEXTRASSHARED_EXPORT Quick3DSpriteSheet : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QQmlListProperty<Qt3DExtras::QSpriteSheetItem> sprites READ sprites CONSTANT)
    Q_CLASSINFO("DefaultProperty", "sprites")
public:
    explicit Quick3DSpriteSheet(QObject *parent = 0);
    ~Quick3DSpriteSheet();

    QQmlListProperty<Qt3DExtras::QSpriteSheetItem> sprites();
    inline QSpriteSheet *parentSpriteSheet() const { return qobject_cast<QSpriteSheet *>(parent()); }

private:
    static void appendSprite(QQmlListProperty<Qt3DExtras::QSpriteSheetItem> *list, Qt3DExtras::QSpriteSheetItem *state);
    static Qt3DExtras::QSpriteSheetItem *spriteAt(QQmlListProperty<Qt3DExtras::QSpriteSheetItem> *list, int index);
    static int spriteCount(QQmlListProperty<Qt3DExtras::QSpriteSheetItem> *list);
    static void clearSprites(QQmlListProperty<Qt3DExtras::QSpriteSheetItem> *list);
};

} // namespace Quick
} // namespace Extras
} // namespace Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DEXTRAS_EXTRAS_QUICK_QUICK3DSPRITESHEET_P_H
