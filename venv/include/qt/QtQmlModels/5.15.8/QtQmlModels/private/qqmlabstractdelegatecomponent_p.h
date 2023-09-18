/****************************************************************************
**
** Copyright (C) 2019 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtQml module of the Qt Toolkit.
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

#ifndef QQMLABSTRACTDELEGATECOMPONENT_P_H
#define QQMLABSTRACTDELEGATECOMPONENT_P_H

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

#include <private/qtqmlmodelsglobal_p.h>
#include <private/qqmlcomponentattached_p.h>
#include <qqmlcomponent.h>

QT_REQUIRE_CONFIG(qml_delegate_model);

QT_BEGIN_NAMESPACE

// TODO: consider making QQmlAbstractDelegateComponent public API
class QQmlAdaptorModel;
class Q_QMLMODELS_PRIVATE_EXPORT QQmlAbstractDelegateComponent : public QQmlComponent
{
    Q_OBJECT
    QML_NAMED_ELEMENT(AbstractDelegateComponent)
    QML_UNCREATABLE("Cannot create instance of abstract class AbstractDelegateComponent.")

public:
    QQmlAbstractDelegateComponent(QObject *parent = nullptr);
    ~QQmlAbstractDelegateComponent() override;

    virtual QQmlComponent *delegate(QQmlAdaptorModel *adaptorModel, int row, int column = 0) const = 0;

signals:
    void delegateChanged();

protected:
    QVariant value(QQmlAdaptorModel *adaptorModel,int row, int column, const QString &role) const;
};

QT_END_NAMESPACE

#endif // QQMLABSTRACTDELEGATECOMPONENT_P_H
