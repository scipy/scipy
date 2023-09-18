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

#ifndef QQUICKOVERLAY_P_H
#define QQUICKOVERLAY_P_H

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

#include <QtQuick/qquickitem.h>
#include <QtQuickTemplates2/private/qquickabstractbutton_p.h>

QT_BEGIN_NAMESPACE

class QQmlComponent;
class QQuickOverlayPrivate;
class QQuickOverlayAttached;
class QQuickOverlayAttachedPrivate;

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickOverlay : public QQuickItem
{
    Q_OBJECT
    Q_PROPERTY(QQmlComponent *modal READ modal WRITE setModal NOTIFY modalChanged FINAL)
    Q_PROPERTY(QQmlComponent *modeless READ modeless WRITE setModeless NOTIFY modelessChanged FINAL)

public:
    explicit QQuickOverlay(QQuickItem *parent = nullptr);
    ~QQuickOverlay();

    QQmlComponent *modal() const;
    void setModal(QQmlComponent *modal);

    QQmlComponent *modeless() const;
    void setModeless(QQmlComponent *modeless);

    static QQuickOverlay *overlay(QQuickWindow *window);

    static QQuickOverlayAttached *qmlAttachedProperties(QObject *object);

Q_SIGNALS:
    void modalChanged();
    void modelessChanged();
    void pressed();
    void released();

protected:
    void itemChange(ItemChange change, const ItemChangeData &data) override;
    void geometryChanged(const QRectF &newGeometry, const QRectF &oldGeometry) override;

    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
#if QT_CONFIG(quicktemplates2_multitouch)
    void touchEvent(QTouchEvent *event) override;
#endif
#if QT_CONFIG(wheelevent)
    void wheelEvent(QWheelEvent *event) override;
#endif
    bool childMouseEventFilter(QQuickItem *item, QEvent *event) override;
    bool eventFilter(QObject *object, QEvent *event) override;

private:
    Q_DISABLE_COPY(QQuickOverlay)
    Q_DECLARE_PRIVATE(QQuickOverlay)
};

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickOverlayAttached : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QQuickOverlay *overlay READ overlay NOTIFY overlayChanged FINAL)
    Q_PROPERTY(QQmlComponent *modal READ modal WRITE setModal NOTIFY modalChanged FINAL)
    Q_PROPERTY(QQmlComponent *modeless READ modeless WRITE setModeless NOTIFY modelessChanged FINAL)

public:
    explicit QQuickOverlayAttached(QObject *parent = nullptr);

    QQuickOverlay *overlay() const;

    QQmlComponent *modal() const;
    void setModal(QQmlComponent *modal);

    QQmlComponent *modeless() const;
    void setModeless(QQmlComponent *modeless);

Q_SIGNALS:
    void overlayChanged();
    void modalChanged();
    void modelessChanged();
    void pressed();
    void released();

private:
    Q_DISABLE_COPY(QQuickOverlayAttached)
    Q_DECLARE_PRIVATE(QQuickOverlayAttached)
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QQuickOverlay)
QML_DECLARE_TYPEINFO(QQuickOverlay, QML_HAS_ATTACHED_PROPERTIES)

#endif // QQUICKOVERLAY_P_H
