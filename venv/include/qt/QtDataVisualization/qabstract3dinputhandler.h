/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Data Visualization module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QABSTRACT3DINPUTHANDLER_H
#define QABSTRACT3DINPUTHANDLER_H

#include <QtDataVisualization/qdatavisualizationglobal.h>
#include <QtDataVisualization/q3dscene.h>
#include <QtCore/QObject>
#include <QtCore/QPoint>
#include <QtGui/QWheelEvent>
#include <QtGui/QMouseEvent>
#include <QtGui/QTouchEvent>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class QAbstract3DInputHandlerPrivate;

class QT_DATAVISUALIZATION_EXPORT QAbstract3DInputHandler : public QObject
{
    Q_OBJECT
    Q_ENUMS(InputView)
    Q_PROPERTY(InputView inputView READ inputView WRITE setInputView NOTIFY inputViewChanged)
    Q_PROPERTY(QPoint inputPosition READ inputPosition WRITE setInputPosition NOTIFY positionChanged)
    Q_PROPERTY(Q3DScene *scene READ scene WRITE setScene NOTIFY sceneChanged)

public:
    enum InputView {
        InputViewNone = 0,
        InputViewOnPrimary,
        InputViewOnSecondary
    };

protected:
    explicit QAbstract3DInputHandler(QObject *parent = nullptr);
public:
    virtual ~QAbstract3DInputHandler();

    // Input event listeners
    virtual void mouseDoubleClickEvent(QMouseEvent *event);
    virtual void touchEvent(QTouchEvent *event);
    virtual void mousePressEvent(QMouseEvent *event, const QPoint &mousePos);
    virtual void mouseReleaseEvent(QMouseEvent *event, const QPoint &mousePos);
    virtual void mouseMoveEvent(QMouseEvent *event, const QPoint &mousePos);
#if QT_CONFIG(wheelevent)
    virtual void wheelEvent(QWheelEvent *event);
#endif

    InputView inputView() const;
    void setInputView(InputView inputView);

    QPoint inputPosition() const;
    void setInputPosition(const QPoint &position);

    Q3DScene *scene() const;
    void setScene(Q3DScene *scene);

Q_SIGNALS:
    void positionChanged(const QPoint &position);
    void inputViewChanged(QAbstract3DInputHandler::InputView view);
    void sceneChanged(Q3DScene *scene);

protected:
    void setPrevDistance(int distance);
    int prevDistance() const;
    void setPreviousInputPos(const QPoint &position);
    QPoint previousInputPos() const;

private:
    Q_DISABLE_COPY(QAbstract3DInputHandler)

    QScopedPointer<QAbstract3DInputHandlerPrivate> d_ptr;

    friend class Abstract3DController;
    friend class QTouch3DInputHandlerPrivate;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
