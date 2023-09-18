/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the QtLocation module of the Qt Toolkit.
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

#ifndef QPARAMETERIZABLEOBJECT_H
#define QPARAMETERIZABLEOBJECT_H

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

#include <QtLocation/private/qlocationglobal_p.h>
#include <QObject>
#include <QQmlListProperty>
#include <QQmlParserStatus>

QT_BEGIN_NAMESPACE
class QGeoMapParameter;
class Q_LOCATION_PRIVATE_EXPORT QParameterizableObject : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QObject *parent READ parent NOTIFY parentChanged DESIGNABLE false FINAL)
    Q_PROPERTY(QQmlListProperty<QObject> quickChildren READ declarativeChildren DESIGNABLE false)
    Q_CLASSINFO("DefaultProperty", "quickChildren")

public:
    explicit QParameterizableObject(QObject *parent = nullptr);

    template <typename T = QObject>
    QList<T*> quickChildren() const
    {
        QList<T*> res;
        for (auto kid : qAsConst(m_children)) {
            auto val = qobject_cast<T*>(kid);
            if (val)
                res.push_back(val);
        }
        return res;
    }

Q_SIGNALS:
    void parentChanged();

protected:
    virtual void appendChild(QObject *v);
    virtual void clearChildren();

    static void append(QQmlListProperty<QObject> *p, QObject *v);
    static int count(QQmlListProperty<QObject> *p);
    static QObject *at(QQmlListProperty<QObject> *p, int idx);
    static void clear(QQmlListProperty<QObject> *p);

    QQmlListProperty<QObject> declarativeChildren();
    QList<QObject*> m_children;
};
QT_END_NAMESPACE

#endif // QPARAMETERIZABLEOBJECT_H
