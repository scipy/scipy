/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
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

#ifndef QQUICKLISTMODELWORKERAGENT_P_H
#define QQUICKLISTMODELWORKERAGENT_P_H

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

#include <qtqmlmodelsglobal_p.h>

#include <QEvent>
#include <QMutex>
#include <QWaitCondition>
#include <QtQml/qqml.h>

#include <private/qv4engine_p.h>

QT_REQUIRE_CONFIG(qml_list_model);

QT_BEGIN_NAMESPACE


class QQmlListModel;

class QQmlListModelWorkerAgent : public QObject
{
    Q_OBJECT
    Q_PROPERTY(int count READ count)
    Q_PROPERTY(QV4::ExecutionEngine *engine READ engine WRITE setEngine NOTIFY engineChanged)
    QML_ANONYMOUS

public:
    QQmlListModelWorkerAgent(QQmlListModel *);
    ~QQmlListModelWorkerAgent();

    QV4::ExecutionEngine *engine() const;
    void setEngine(QV4::ExecutionEngine *eng);

    Q_INVOKABLE void addref();
    Q_INVOKABLE void release();

    int count() const;

    Q_INVOKABLE void clear();
    Q_INVOKABLE void remove(QQmlV4Function *args);
    Q_INVOKABLE void append(QQmlV4Function *args);
    Q_INVOKABLE void insert(QQmlV4Function *args);
    Q_INVOKABLE QJSValue get(int index) const;
    Q_INVOKABLE void set(int index, const QJSValue &value);
    Q_INVOKABLE void setProperty(int index, const QString& property, const QVariant& value);
    Q_INVOKABLE void move(int from, int to, int count);
    Q_INVOKABLE void sync();

    void modelDestroyed();

signals:
    void engineChanged(QV4::ExecutionEngine *engine);

protected:
    bool event(QEvent *) override;

private:
    friend class QQuickWorkerScriptEnginePrivate;
    friend class QQmlListModel;

    struct Sync : public QEvent {
        Sync(QQmlListModel *l)
            : QEvent(QEvent::User)
            , list(l)
        {}
        ~Sync();
        QQmlListModel *list;
    };

    QAtomicInt m_ref;
    QQmlListModel *m_orig;
    QQmlListModel *m_copy;
    QMutex mutex;
    QWaitCondition syncDone;
};

QT_END_NAMESPACE

#endif // QQUICKLISTMODELWORKERAGENT_P_H

