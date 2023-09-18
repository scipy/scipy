/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtScxml module of the Qt Toolkit.
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

#ifndef QSCXMLINVOKABLESERVICE_P_H
#define QSCXMLINVOKABLESERVICE_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include "qscxmlinvokableservice.h"
#include <QtCore/private/qobject_p.h>

QT_BEGIN_NAMESPACE

class QScxmlInvokableServicePrivate : public QObjectPrivate
{
public:
    QScxmlInvokableServicePrivate(QScxmlStateMachine *parentStateMachine);

    QString calculateId(QScxmlStateMachine *parent,
                        const QScxmlExecutableContent::InvokeInfo &invokeInfo, bool *ok) const;
    QVariantMap calculateData(QScxmlStateMachine *parent,
                              const QVector<QScxmlExecutableContent::ParameterInfo> &parameters,
                              const QVector<QScxmlExecutableContent::StringId> &names,
                              bool *ok) const;

    QScxmlStateMachine *parentStateMachine;
};

class QScxmlInvokableServiceFactoryPrivate : public QObjectPrivate
{
public:
    QScxmlInvokableServiceFactoryPrivate(
            const QScxmlExecutableContent::InvokeInfo &invokeInfo,
            const QVector<QScxmlExecutableContent::StringId> &names,
            const QVector<QScxmlExecutableContent::ParameterInfo> &parameters);

    QScxmlExecutableContent::InvokeInfo invokeInfo;
    QVector<QScxmlExecutableContent::StringId> names;
    QVector<QScxmlExecutableContent::ParameterInfo> parameters;
};

class Q_SCXML_EXPORT QScxmlScxmlService: public QScxmlInvokableService
{
    Q_OBJECT
    Q_DECLARE_PRIVATE(QScxmlInvokableService)
    Q_PROPERTY(QScxmlStateMachine *stateMachine READ stateMachine CONSTANT)
public:
    QScxmlScxmlService(QScxmlStateMachine *stateMachine,
                       QScxmlStateMachine *parentStateMachine,
                       QScxmlInvokableServiceFactory *parent);
    ~QScxmlScxmlService();

    bool start() override;
    QString id() const override;
    QString name() const override;
    void postEvent(QScxmlEvent *event) override;
    QScxmlStateMachine *stateMachine() const;

private:
    QScxmlStateMachine *m_stateMachine;
};

class QScxmlStaticScxmlServiceFactoryPrivate : public QScxmlInvokableServiceFactoryPrivate
{
public:
    QScxmlStaticScxmlServiceFactoryPrivate(
            const QMetaObject *metaObject,
            const QScxmlExecutableContent::InvokeInfo &invokeInfo,
            const QVector<QScxmlExecutableContent::StringId> &names,
            const QVector<QScxmlExecutableContent::ParameterInfo> &parameters);

    const QMetaObject *metaObject;
};

QScxmlScxmlService *invokeDynamicScxmlService(const QString &sourceUrl,
                                              QScxmlStateMachine *parentStateMachine,
                                              QScxmlInvokableServiceFactory *factory);
QScxmlScxmlService *invokeStaticScxmlService(QScxmlStateMachine *childStateMachine,
                                             QScxmlStateMachine *parentStateMachine,
                                             QScxmlInvokableServiceFactory *factory);
QString calculateSrcexpr(QScxmlStateMachine *parent, QScxmlExecutableContent::EvaluatorId srcexpr,
                         bool *ok);

QT_END_NAMESPACE

#endif // QSCXMLINVOKABLESERVICE_P_H
