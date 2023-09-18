/****************************************************************************
**
** Copyright (C) 2017 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DCORE_QNODECOMMAND_H
#define QT3DCORE_QNODECOMMAND_H

#include <Qt3DCore/qt3dcore_global.h>
#include <Qt3DCore/qscenechange.h>

#include <QtCore/qsharedpointer.h>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

class QNodeCommandPrivate;

class Q_3DCORESHARED_EXPORT QNodeCommand : public QSceneChange
{
public:
#if defined(Q_ATOMIC_INT64_IS_SUPPORTED)
    typedef quint64 CommandId;
#else
    typedef quint32 CommandId;
#endif

    Q3D_DECL_DEPRECATED explicit QNodeCommand(QNodeId id);
    ~QNodeCommand();

    CommandId commandId() const;

    QString name() const;
    void setName(const QString &name);
    QVariant data() const;
    void setData(const QVariant &data);
    CommandId inReplyTo() const;
    void setReplyToCommandId(CommandId id);

protected:
    explicit QNodeCommand(QNodeCommandPrivate &dd, QNodeId id);

private:
    Q_DECLARE_PRIVATE(QNodeCommand)
};

typedef QSharedPointer<QNodeCommand> QNodeCommandPtr;

} // namespace Qt3DCore

QT_END_NAMESPACE

#endif // QT3DCORE_QNODECOMMAND_H
