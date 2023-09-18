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

#ifndef QT3DCORE_QNODECREATEDCHANGEGENERATOR_H
#define QT3DCORE_QNODECREATEDCHANGEGENERATOR_H

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

#include <Qt3DCore/qnode.h>
#include <Qt3DCore/qnodecreatedchange.h>
#include <Qt3DCore/qt3dcore_global.h>
#include <QtCore/qvector.h>

#include <Qt3DCore/private/qnode_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

class Q_3DCORE_PRIVATE_EXPORT QNodeCreatedChangeGenerator
{
public:
    QNodeCreatedChangeGenerator(QNode *rootNode);

    QVector<QNodeCreatedChangeBasePtr> creationChanges() const { return m_creationChanges; }

private:
    void createCreationChange(QNode *node)
    {
        QT_WARNING_PUSH
        QT_WARNING_DISABLE_DEPRECATED
        const auto creationChange = node->createNodeCreationChange();
        m_creationChanges.push_back(creationChange);

        // Store the metaobject of the node in the QNode so that we have it available
        // to us during destruction in the QNode destructor. This allows us to send
        // the QNodeId and the metaobject as typeinfo to the backend aspects so they
        // in turn can find the correct QBackendNodeMapper object to handle the destruction
        // of the corresponding backend nodes.
        QNodePrivate *d = QNodePrivate::get(node);
        d->m_typeInfo = const_cast<QMetaObject*>(creationChange->metaObject());

        // Mark this node as having been handled for creation so that it is picked up
        d->m_hasBackendNode = true;
        QT_WARNING_POP
    }

    QVector<QNodeCreatedChangeBasePtr> m_creationChanges;
};

} // namespace Qt3DCore

QT_END_NAMESPACE

#endif // QT3DCORE_QNODECREATEDCHANGEGENERATOR_H
