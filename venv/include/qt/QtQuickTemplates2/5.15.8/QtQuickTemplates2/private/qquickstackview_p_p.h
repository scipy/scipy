/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the Qt Quick Templates 2 module of the Qt Toolkit.
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

#ifndef QQUICKSTACKVIEW_P_P_H
#define QQUICKSTACKVIEW_P_P_H

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

#include <QtQuickTemplates2/private/qquickstackview_p.h>
#include <QtQuickTemplates2/private/qquickcontrol_p_p.h>
#include <QtQuick/private/qquickitemviewtransition_p.h>
#include <QtQuick/private/qquickitemchangelistener_p.h>
#include <QtQml/private/qv4value_p.h>
#include <QtCore/qset.h>
#include <QtCore/qstack.h>

QT_BEGIN_NAMESPACE

class QQmlContextData;
class QQuickStackElement;
struct QQuickStackTransition;

class QQuickStackViewPrivate : public QQuickControlPrivate, public QQuickItemViewTransitionChangeListener
{
    Q_DECLARE_PUBLIC(QQuickStackView)

public:
    static QQuickStackViewPrivate *get(QQuickStackView *view)
    {
        return view->d_func();
    }

    void warn(const QString &error);
    void warnOfInterruption(const QString &attemptedOperation);

    void setCurrentItem(QQuickStackElement *element);

    QList<QQuickStackElement *> parseElements(int from, QQmlV4Function *args, QStringList *errors);
    QQuickStackElement *findElement(QQuickItem *item) const;
    QQuickStackElement *findElement(const QV4::Value &value) const;
    QQuickStackElement *createElement(const QV4::Value &value, QQmlContextData *context, QString *error);
    bool pushElements(const QList<QQuickStackElement *> &elements);
    bool pushElement(QQuickStackElement *element);
    bool popElements(QQuickStackElement *element);
    bool replaceElements(QQuickStackElement *element, const QList<QQuickStackElement *> &elements);

    void ensureTransitioner();
    void startTransition(const QQuickStackTransition &first, const QQuickStackTransition &second, bool immediate);
    void completeTransition(QQuickStackElement *element, QQuickTransition *transition, QQuickStackView::Status status);

    void viewItemTransitionFinished(QQuickItemViewTransitionableItem *item) override;
    void setBusy(bool busy);
    void depthChange(int newDepth, int oldDepth);

    bool busy = false;
    bool modifyingElements = false;
    QString operation;
    QJSValue initialItem;
    QQuickItem *currentItem = nullptr;
    QSet<QQuickStackElement*> removing;
    QList<QQuickStackElement*> removed;
    QStack<QQuickStackElement *> elements;
    QQuickItemViewTransitioner *transitioner = nullptr;
};

class QQuickStackViewAttachedPrivate : public QObjectPrivate, public QQuickItemChangeListener
{
    Q_DECLARE_PUBLIC(QQuickStackViewAttached)

public:
    static QQuickStackViewAttachedPrivate *get(QQuickStackViewAttached *attached)
    {
        return attached->d_func();
    }

    void itemParentChanged(QQuickItem *item, QQuickItem *parent) override;

    bool explicitVisible = false;
    QQuickStackElement *element = nullptr;
};

QT_END_NAMESPACE

#endif // QQUICKSTACKVIEW_P_P_H
