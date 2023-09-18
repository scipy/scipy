/****************************************************************************
**
** Copyright (C) 2015 Klaralvdalens Datakonsult AB (KDAB).
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
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

#ifndef QT3DINPUT_QMOUSEEVENT_H
#define QT3DINPUT_QMOUSEEVENT_H

#include <Qt3DInput/qt3dinput_global.h>
#include <QtCore/QtGlobal>
#include <QtCore/QObject>
#include <QtGui/QMouseEvent>

QT_BEGIN_NAMESPACE

namespace Qt3DInput {

class Q_3DINPUTSHARED_EXPORT QMouseEvent : public QObject
{
    Q_OBJECT
    Q_PROPERTY(int x READ x CONSTANT)
    Q_PROPERTY(int y READ y CONSTANT)
    Q_PROPERTY(bool wasHeld READ wasHeld CONSTANT)
    Q_PROPERTY(Qt3DInput::QMouseEvent::Buttons button READ button CONSTANT)
    Q_PROPERTY(int buttons READ buttons CONSTANT)
    Q_PROPERTY(Qt3DInput::QMouseEvent::Modifiers modifiers READ modifiers CONSTANT)
    Q_PROPERTY(bool accepted READ isAccepted WRITE setAccepted)

public:
    enum Buttons {
        LeftButton = Qt::LeftButton,
        RightButton = Qt::RightButton,
        MiddleButton = Qt::MiddleButton,
        BackButton = Qt::BackButton,
        NoButton = Qt::NoButton
    };
    Q_ENUM(Buttons) // LCOV_EXCL_LINE

    enum Modifiers {
        NoModifier = Qt::NoModifier,
        ShiftModifier = Qt::ShiftModifier,
        ControlModifier = Qt::ControlModifier,
        AltModifier = Qt::AltModifier,
        MetaModifier = Qt::MetaModifier,
        KeypadModifier = Qt::KeypadModifier
    };
    Q_ENUM(Modifiers) // LCOV_EXCL_LINE
    // TO DO Qt6 Modifiers -> Modifier and add Q_FLAG(Modifiers)

    explicit QMouseEvent(const QT_PREPEND_NAMESPACE(QMouseEvent) &e);
    ~QMouseEvent();

    inline int x() const { return m_event.x(); }
    inline int y() const { return m_event.y(); }
    inline bool wasHeld() const {
#if QT_CONFIG(gestures)
      return static_cast<Qt::GestureType>(m_event.type()) == Qt::TapAndHoldGesture;
#else
      return false;
#endif
    }
    Buttons button() const;
    int buttons() const;
    Modifiers modifiers() const;

    inline bool isAccepted() const { return m_event.isAccepted(); }
    inline void setAccepted(bool accepted) { m_event.setAccepted(accepted); }
    inline QEvent::Type type() const { return m_event.type(); }

private:
    QT_PREPEND_NAMESPACE(QMouseEvent) m_event;
};

typedef QSharedPointer<QMouseEvent> QMouseEventPtr;

#if QT_CONFIG(wheelevent)
class Q_3DINPUTSHARED_EXPORT QWheelEvent : public QObject
{
    Q_OBJECT
    Q_PROPERTY(int x READ x CONSTANT)
    Q_PROPERTY(int y READ y CONSTANT)
    Q_PROPERTY(QPoint angleDelta READ angleDelta CONSTANT)
    Q_PROPERTY(int buttons READ buttons CONSTANT)
    Q_PROPERTY(Qt3DInput::QWheelEvent::Modifiers modifiers READ modifiers CONSTANT)
    Q_PROPERTY(bool accepted READ isAccepted WRITE setAccepted)

public:
    enum Buttons {
        LeftButton = Qt::LeftButton,
        RightButton = Qt::RightButton,
        MiddleButton = Qt::MiddleButton,
        BackButton = Qt::BackButton,
        NoButton = Qt::NoButton
    };
    Q_ENUM(Buttons) // LCOV_EXCL_LINE

    enum Modifiers {
        NoModifier = Qt::NoModifier,
        ShiftModifier = Qt::ShiftModifier,
        ControlModifier = Qt::ControlModifier,
        AltModifier = Qt::AltModifier,
        MetaModifier = Qt::MetaModifier,
        KeypadModifier = Qt::KeypadModifier
    };
    Q_ENUM(Modifiers) // LCOV_EXCL_LINE
    // TO DO Qt6 Modifiers -> Modifier and add Q_FLAG(Modifiers)

    explicit QWheelEvent(const QT_PREPEND_NAMESPACE(QWheelEvent) &e);
    ~QWheelEvent();

    inline int x() const { return int(m_event.position().x()); }
    inline int y() const { return int(m_event.position().y()); }
    inline QPoint angleDelta() const { return m_event.angleDelta(); }
    int buttons() const;
    Modifiers modifiers() const;

    inline bool isAccepted() const { return m_event.isAccepted(); }
    inline void setAccepted(bool accepted) { m_event.setAccepted(accepted); }
    inline QEvent::Type type() const { return m_event.type(); }

private:
    QT_PREPEND_NAMESPACE(QWheelEvent) m_event;
};

typedef QSharedPointer<QWheelEvent> QWheelEventPtr;
#endif

} // namespace Qt3DInput

QT_END_NAMESPACE

Q_DECLARE_METATYPE(Qt3DInput::QMouseEvent*) // LCOV_EXCL_LINE

#if QT_CONFIG(wheelevent)
Q_DECLARE_METATYPE(Qt3DInput::QWheelEvent*) // LCOV_EXCL_LINE
#endif

#endif // QT3DINPUT_QMOUSEEVENT_H
