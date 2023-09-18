/****************************************************************************
**
** Copyright (C) 2016 Research In Motion.
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

#ifndef QQMLINSTANTIATOR_P_H
#define QQMLINSTANTIATOR_P_H

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

#include <QtQml/qqmlcomponent.h>
#include <QtQml/qqmlparserstatus.h>
#include <QtQmlModels/private/qtqmlmodelsglobal_p.h>

QT_REQUIRE_CONFIG(qml_object_model);

QT_BEGIN_NAMESPACE

class QQmlInstantiatorPrivate;
class Q_QMLMODELS_PRIVATE_EXPORT QQmlInstantiator : public QObject, public QQmlParserStatus
{
    Q_OBJECT
    Q_INTERFACES(QQmlParserStatus)

    Q_PROPERTY(bool active READ isActive WRITE setActive NOTIFY activeChanged)
    Q_PROPERTY(bool asynchronous READ isAsync WRITE setAsync NOTIFY asynchronousChanged)
    Q_PROPERTY(QVariant model READ model WRITE setModel NOTIFY modelChanged)
    Q_PROPERTY(int count READ count NOTIFY countChanged)
    Q_PROPERTY(QQmlComponent *delegate READ delegate WRITE setDelegate NOTIFY delegateChanged)
    Q_PROPERTY(QObject *object READ object NOTIFY objectChanged)
    Q_CLASSINFO("DefaultProperty", "delegate")
    QML_NAMED_ELEMENT(Instantiator)
    QML_ADDED_IN_MINOR_VERSION(14)

public:
    QQmlInstantiator(QObject *parent = nullptr);
    ~QQmlInstantiator();

    bool isActive() const;
    void setActive(bool newVal);

    bool isAsync() const;
    void setAsync(bool newVal);

    int count() const;

    QQmlComponent* delegate();
    void setDelegate(QQmlComponent* c);

    QVariant model() const;
    void setModel(const QVariant &v);

    QObject *object() const;

    Q_INVOKABLE QObject *objectAt(int index) const;

    void classBegin() override;
    void componentComplete() override;

Q_SIGNALS:
    void modelChanged();
    void delegateChanged();
    void countChanged();
    void objectChanged();
    void activeChanged();
    void asynchronousChanged();

    void objectAdded(int index, QObject* object);
    void objectRemoved(int index, QObject* object);

private:
    Q_DISABLE_COPY(QQmlInstantiator)
    Q_DECLARE_PRIVATE(QQmlInstantiator)
    Q_PRIVATE_SLOT(d_func(), void _q_createdItem(int, QObject *))
    Q_PRIVATE_SLOT(d_func(), void _q_modelUpdated(const QQmlChangeSet &, bool))
};

QT_END_NAMESPACE

#endif // QQMLCREATOR_P_H
