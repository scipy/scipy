/****************************************************************************
**
** Copyright (C) 2016 Klaralvdalens Datakonsult AB (KDAB).
** Contact: http://www.qt-project.org/legal
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

#ifndef QT3DCORE_QNODECREATEDCHANGE_H
#define QT3DCORE_QNODECREATEDCHANGE_H

#include <Qt3DCore/qscenechange.h>
#include <Qt3DCore/qt3dcore_global.h>
#include <QtCore/qsharedpointer.h>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

class QNode;
class QNodeCreatedChangeBasePrivate;

class Q_3DCORESHARED_EXPORT QNodeCreatedChangeBase : public QSceneChange
{
public:
    Q3D_DECL_DEPRECATED explicit QNodeCreatedChangeBase(const QNode *node);
    ~QNodeCreatedChangeBase();

    QNodeId parentId() const Q_DECL_NOTHROW;
    const QMetaObject *metaObject() const Q_DECL_NOTHROW;
    bool isNodeEnabled() const Q_DECL_NOTHROW;

protected:
    QNodeCreatedChangeBase(QNodeCreatedChangeBasePrivate &dd, const QNode *node);

private:
    Q_DECLARE_PRIVATE(QNodeCreatedChangeBase)
};

typedef QSharedPointer<QNodeCreatedChangeBase> QNodeCreatedChangeBasePtr;

template<typename T>
class QNodeCreatedChange : public QNodeCreatedChangeBase
{
public:
    QT_WARNING_PUSH
    QT_WARNING_DISABLE_DEPRECATED
    explicit QNodeCreatedChange(const QNode *_node)
        : QNodeCreatedChangeBase(_node)
        , data()
    {
    }
    QT_WARNING_POP

    T data;
};

template<typename T>
using QNodeCreatedChangePtr = QSharedPointer<QNodeCreatedChange<T>>;

} // namespace Qt3DCore

QT_END_NAMESPACE

#endif // QT3DCORE_QNODECREATEDCHANGE_H
