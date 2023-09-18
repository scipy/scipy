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

#ifndef QT3DRENDER_QPICKEVENT_H
#define QT3DRENDER_QPICKEVENT_H

#include <QtCore/QObject>
#include <QtGui/QVector3D>
#include <QtCore/QPointF>
#include <Qt3DRender/qt3drender_global.h>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {
class QEntity;
}

namespace Qt3DRender {

class QViewport;
class QPickEventPrivate;

class QPickEvent;
typedef QSharedPointer<QPickEvent> QPickEventPtr;

class Q_3DRENDERSHARED_EXPORT QPickEvent : public QObject
{
    Q_OBJECT
    Q_PROPERTY(bool accepted READ isAccepted WRITE setAccepted NOTIFY acceptedChanged)
    Q_PROPERTY(QPointF position READ position CONSTANT)
    Q_PROPERTY(float distance READ distance CONSTANT)
    Q_PROPERTY(QVector3D localIntersection READ localIntersection CONSTANT)
    Q_PROPERTY(QVector3D worldIntersection READ worldIntersection CONSTANT)
    Q_PROPERTY(Qt3DRender::QPickEvent::Buttons button READ button CONSTANT)
    Q_PROPERTY(int buttons READ buttons CONSTANT)
    Q_PROPERTY(int modifiers READ modifiers CONSTANT)
    Q_PROPERTY(Qt3DRender::QViewport *viewport READ viewport CONSTANT REVISION 14)
    Q_PROPERTY(Qt3DCore::QEntity *entity READ entity CONSTANT REVISION 14)
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

    QPickEvent();
    QPickEvent(const QPointF &position, const QVector3D& worldIntersection, const QVector3D& localIntersection, float distance);
    QPickEvent(const QPointF &position, const QVector3D& worldIntersection, const QVector3D& localIntersection, float distance, Buttons button,
               int buttons, int modifiers);
    ~QPickEvent();

    bool isAccepted() const;

public Q_SLOTS:
    void setAccepted(bool accepted);

public:
    QPointF position() const;
    float distance() const;
    QVector3D worldIntersection() const;
    QVector3D localIntersection() const;
    Buttons button() const;
    int buttons() const;
    int modifiers() const;
    QViewport *viewport() const;
    Qt3DCore::QEntity *entity() const;

Q_SIGNALS:
    void acceptedChanged(bool accepted);

protected:
    explicit QPickEvent(QObjectPrivate &dd, QObject *parent = nullptr);

private:
    Q_DECLARE_PRIVATE(QPickEvent)

    friend class QObjectPickerPrivate;
};

} // Qt3DRender

QT_END_NAMESPACE

Q_DECLARE_METATYPE(Qt3DRender::QPickEvent*) // LCOV_EXCL_LINE

#endif // QT3DRENDER_QPICKEVENT_H
