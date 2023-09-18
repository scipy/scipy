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

#ifndef QT3DCORE_QENTITY_H
#define QT3DCORE_QENTITY_H

#include <Qt3DCore/qnode.h>
#include <Qt3DCore/qt3dcore_global.h>
#include <QtCore/QMetaType>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

class QComponent;
class QEntityPrivate;

typedef QVector<QComponent*> QComponentVector;

class Q_3DCORESHARED_EXPORT QEntity : public QNode
{
    Q_OBJECT
public:
    explicit QEntity(QNode *parent = nullptr);
    virtual ~QEntity();

    QComponentVector components() const;

    template<class T>
    QVector<T *> componentsOfType() const
    {
        QVector<T*> matchComponents;
        const QComponentVector comps = this->components();
        for (QComponent *component : comps) {
            T *typedComponent = qobject_cast<T*>(component);
            if (typedComponent != nullptr)
                matchComponents.append(typedComponent);
        }
        return matchComponents;
    }

    void addComponent(QComponent *comp);
    void removeComponent(QComponent *comp);

    QEntity *parentEntity() const;

protected:
    explicit QEntity(QEntityPrivate &dd, QNode *parent = nullptr);

private Q_SLOTS:
    void onParentChanged(QObject *);

private:
    Q_DECLARE_PRIVATE(QEntity)

    QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

typedef QSharedPointer<QEntity> QEntityPtr;

} // namespace Qt3DCore

QT_END_NAMESPACE

#endif // QT3DCORE_QENTITY_H
