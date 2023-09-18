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

#ifndef QQUICKSCROLLBAR_P_H
#define QQUICKSCROLLBAR_P_H

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

class QQuickScrollBarAttached;
class QQuickScrollBarPrivate;

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickScrollBar : public QQuickControl
{
    Q_OBJECT
    Q_PROPERTY(qreal size READ size WRITE setSize NOTIFY sizeChanged FINAL)
    Q_PROPERTY(qreal position READ position WRITE setPosition NOTIFY positionChanged FINAL)
    Q_PROPERTY(qreal stepSize READ stepSize WRITE setStepSize NOTIFY stepSizeChanged FINAL)
    Q_PROPERTY(bool active READ isActive WRITE setActive NOTIFY activeChanged FINAL)
    Q_PROPERTY(bool pressed READ isPressed WRITE setPressed NOTIFY pressedChanged FINAL)
    Q_PROPERTY(Qt::Orientation orientation READ orientation WRITE setOrientation NOTIFY orientationChanged FINAL)
    // 2.2 (Qt 5.9)
    Q_PROPERTY(SnapMode snapMode READ snapMode WRITE setSnapMode NOTIFY snapModeChanged FINAL REVISION 2)
    Q_PROPERTY(bool interactive READ isInteractive WRITE setInteractive RESET resetInteractive NOTIFY interactiveChanged FINAL REVISION 2)
    Q_PROPERTY(Policy policy READ policy WRITE setPolicy NOTIFY policyChanged FINAL REVISION 2)
    // 2.3 (Qt 5.10)
    Q_PROPERTY(bool horizontal READ isHorizontal NOTIFY orientationChanged FINAL REVISION 3)
    Q_PROPERTY(bool vertical READ isVertical NOTIFY orientationChanged FINAL REVISION 3)
    // 2.4 (Qt 5.11)
    Q_PROPERTY(qreal minimumSize READ minimumSize WRITE setMinimumSize NOTIFY minimumSizeChanged FINAL REVISION 4)
    Q_PROPERTY(qreal visualSize READ visualSize NOTIFY visualSizeChanged FINAL REVISION 4)
    Q_PROPERTY(qreal visualPosition READ visualPosition NOTIFY visualPositionChanged FINAL REVISION 4)

public:
    explicit QQuickScrollBar(QQuickItem *parent = nullptr);

    static QQuickScrollBarAttached *qmlAttachedProperties(QObject *object);

    qreal size() const;
    qreal position() const;

    qreal stepSize() const;
    void setStepSize(qreal step);

    bool isActive() const;
    void setActive(bool active);

    bool isPressed() const;
    void setPressed(bool pressed);

    Qt::Orientation orientation() const;
    void setOrientation(Qt::Orientation orientation);

    // 2.2 (Qt 5.9)
    enum SnapMode {
        NoSnap,
        SnapAlways,
        SnapOnRelease
    };
    Q_ENUM(SnapMode)

    SnapMode snapMode() const;
    void setSnapMode(SnapMode mode);

    bool isInteractive() const;
    void setInteractive(bool interactive);
    void resetInteractive();

    enum Policy {
        AsNeeded = Qt::ScrollBarAsNeeded,
        AlwaysOff = Qt::ScrollBarAlwaysOff,
        AlwaysOn = Qt::ScrollBarAlwaysOn
    };
    Q_ENUM(Policy)

    Policy policy() const;
    void setPolicy(Policy policy);

    // 2.3 (Qt 5.10)
    bool isHorizontal() const;
    bool isVertical() const;

    // 2.4 (Qt 5.11)
    qreal minimumSize() const;
    void setMinimumSize(qreal minimumSize);

    qreal visualSize() const;
    qreal visualPosition() const;

public Q_SLOTS:
    void increase();
    void decrease();
    void setSize(qreal size);
    void setPosition(qreal position);

Q_SIGNALS:
    void sizeChanged();
    void positionChanged();
    void stepSizeChanged();
    void activeChanged();
    void pressedChanged();
    void orientationChanged();
    // 2.2 (Qt 5.9)
    Q_REVISION(2) void snapModeChanged();
    Q_REVISION(2) void interactiveChanged();
    Q_REVISION(2) void policyChanged();
    // 2.4 (Qt 5.11)
    Q_REVISION(4) void minimumSizeChanged();
    Q_REVISION(4) void visualSizeChanged();
    Q_REVISION(4) void visualPositionChanged();

protected:
    void mousePressEvent(QMouseEvent *event) override;

#if QT_CONFIG(quicktemplates2_hover)
    void hoverChange() override;
#endif

#if QT_CONFIG(accessibility)
    void accessibilityActiveChanged(bool active) override;
    QAccessible::Role accessibleRole() const override;
#endif

private:
    Q_DISABLE_COPY(QQuickScrollBar)
    Q_DECLARE_PRIVATE(QQuickScrollBar)
};

class QQuickScrollBarAttachedPrivate;

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickScrollBarAttached : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QQuickScrollBar *horizontal READ horizontal WRITE setHorizontal NOTIFY horizontalChanged FINAL)
    Q_PROPERTY(QQuickScrollBar *vertical READ vertical WRITE setVertical NOTIFY verticalChanged FINAL)

public:
    explicit QQuickScrollBarAttached(QObject *parent = nullptr);
    ~QQuickScrollBarAttached();

    QQuickScrollBar *horizontal() const;
    void setHorizontal(QQuickScrollBar *horizontal);

    QQuickScrollBar *vertical() const;
    void setVertical(QQuickScrollBar *vertical);

Q_SIGNALS:
    void horizontalChanged();
    void verticalChanged();

private:
    Q_DISABLE_COPY(QQuickScrollBarAttached)
    Q_DECLARE_PRIVATE(QQuickScrollBarAttached)
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QQuickScrollBar)
QML_DECLARE_TYPEINFO(QQuickScrollBar, QML_HAS_ATTACHED_PROPERTIES)

#endif // QQUICKSCROLLBAR_P_H
