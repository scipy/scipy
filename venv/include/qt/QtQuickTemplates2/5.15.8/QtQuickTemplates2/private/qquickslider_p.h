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

#ifndef QQUICKSLIDER_P_H
#define QQUICKSLIDER_P_H

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

#include <QtQuickTemplates2/private/qquickcontrol_p.h>

QT_BEGIN_NAMESPACE

class QQuickSliderPrivate;

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickSlider : public QQuickControl
{
    Q_OBJECT
    Q_PROPERTY(qreal from READ from WRITE setFrom NOTIFY fromChanged FINAL)
    Q_PROPERTY(qreal to READ to WRITE setTo NOTIFY toChanged FINAL)
    Q_PROPERTY(qreal value READ value WRITE setValue NOTIFY valueChanged FINAL)
    Q_PROPERTY(qreal position READ position NOTIFY positionChanged FINAL)
    Q_PROPERTY(qreal visualPosition READ visualPosition NOTIFY visualPositionChanged FINAL)
    Q_PROPERTY(qreal stepSize READ stepSize WRITE setStepSize NOTIFY stepSizeChanged FINAL)
    Q_PROPERTY(SnapMode snapMode READ snapMode WRITE setSnapMode NOTIFY snapModeChanged FINAL)
    Q_PROPERTY(bool pressed READ isPressed WRITE setPressed NOTIFY pressedChanged FINAL)
    Q_PROPERTY(Qt::Orientation orientation READ orientation WRITE setOrientation NOTIFY orientationChanged FINAL)
    Q_PROPERTY(QQuickItem *handle READ handle WRITE setHandle NOTIFY handleChanged FINAL)
    Q_PROPERTY(bool live READ live WRITE setLive NOTIFY liveChanged FINAL REVISION 2)
    // 2.3 (Qt 5.10)
    Q_PROPERTY(bool horizontal READ isHorizontal NOTIFY orientationChanged FINAL REVISION 3)
    Q_PROPERTY(bool vertical READ isVertical NOTIFY orientationChanged FINAL REVISION 3)
    // 2.5 (Qt 5.12)
    Q_PROPERTY(qreal touchDragThreshold READ touchDragThreshold WRITE setTouchDragThreshold RESET resetTouchDragThreshold NOTIFY touchDragThresholdChanged FINAL REVISION 5)
    Q_PROPERTY(qreal implicitHandleWidth READ implicitHandleWidth NOTIFY implicitHandleWidthChanged FINAL REVISION 5)
    Q_PROPERTY(qreal implicitHandleHeight READ implicitHandleHeight NOTIFY implicitHandleHeightChanged FINAL REVISION 5)
    Q_CLASSINFO("DeferredPropertyNames", "background,handle")

public:
    explicit QQuickSlider(QQuickItem *parent = nullptr);
    ~QQuickSlider();

    qreal from() const;
    void setFrom(qreal from);

    qreal to() const;
    void setTo(qreal to);

    qreal value() const;
    void setValue(qreal value);

    qreal position() const;
    qreal visualPosition() const;

    qreal stepSize() const;
    void setStepSize(qreal step);

    enum SnapMode {
        NoSnap,
        SnapAlways,
        SnapOnRelease
    };
    Q_ENUM(SnapMode)

    SnapMode snapMode() const;
    void setSnapMode(SnapMode mode);

    bool isPressed() const;
    void setPressed(bool pressed);

    Qt::Orientation orientation() const;
    void setOrientation(Qt::Orientation orientation);

    QQuickItem *handle() const;
    void setHandle(QQuickItem *handle);

    // 2.1 (Qt 5.8)
    Q_REVISION(1) Q_INVOKABLE qreal valueAt(qreal position) const;

    // 2.2 (Qt 5.9)
    bool live() const;
    void setLive(bool live);

    // 2.3 (Qt 5.10)
    bool isHorizontal() const;
    bool isVertical() const;

    // 2.5 (Qt 5.12)
    qreal touchDragThreshold() const;
    void setTouchDragThreshold(qreal touchDragThreshold);
    void resetTouchDragThreshold();

    qreal implicitHandleWidth() const;
    qreal implicitHandleHeight() const;

public Q_SLOTS:
    void increase();
    void decrease();

Q_SIGNALS:
    void fromChanged();
    void toChanged();
    void valueChanged();
    void positionChanged();
    void visualPositionChanged();
    void stepSizeChanged();
    void snapModeChanged();
    void pressedChanged();
    void orientationChanged();
    void handleChanged();
    // 2.2 (Qt 5.9)
    Q_REVISION(2) void moved();
    Q_REVISION(2) void liveChanged();
    // 2.5 (Qt 5.12)
    Q_REVISION(5) void touchDragThresholdChanged();
    Q_REVISION(5) void implicitHandleWidthChanged();
    Q_REVISION(5) void implicitHandleHeightChanged();

protected:
    void keyPressEvent(QKeyEvent *event) override;
    void keyReleaseEvent(QKeyEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
#if QT_CONFIG(quicktemplates2_multitouch)
    void touchEvent(QTouchEvent *event) override;
#endif
#if QT_CONFIG(wheelevent)
    void wheelEvent(QWheelEvent *event) override;
#endif

    void mirrorChange() override;
    void componentComplete() override;

#if QT_CONFIG(accessibility)
    void accessibilityActiveChanged(bool active) override;
    QAccessible::Role accessibleRole() const override;
#endif

private:
    Q_DISABLE_COPY(QQuickSlider)
    Q_DECLARE_PRIVATE(QQuickSlider)
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QQuickSlider)

#endif // QQUICKSLIDER_P_H
