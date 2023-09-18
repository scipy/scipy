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

#ifndef QQUICKPOPUP_P_P_H
#define QQUICKPOPUP_P_P_H

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

#include <QtQuickTemplates2/private/qquickpopup_p.h>
#include <QtQuickTemplates2/private/qquickcontrol_p.h>

#include <QtCore/private/qobject_p.h>
#include <QtQuick/qquickitem.h>
#include <QtQuick/private/qquickitemchangelistener_p.h>
#include <QtQuick/private/qquicktransitionmanager_p_p.h>

QT_BEGIN_NAMESPACE

class QQuickTransition;
class QQuickTransitionManager;
class QQuickPopup;
class QQuickPopupAnchors;
class QQuickPopupItem;
class QQuickPopupPrivate;
class QQuickPopupPositioner;

class QQuickPopupTransitionManager : public QQuickTransitionManager
{
public:
    QQuickPopupTransitionManager(QQuickPopupPrivate *popup);

    void transitionEnter();
    void transitionExit();

protected:
    void finished() override;

private:
    QQuickPopupPrivate *popup = nullptr;
};

class Q_AUTOTEST_EXPORT QQuickPopupPrivate : public QObjectPrivate, public QQuickItemChangeListener
{
    Q_DECLARE_PUBLIC(QQuickPopup)

public:
    QQuickPopupPrivate();

    static QQuickPopupPrivate *get(QQuickPopup *popup)
    {
        return popup->d_func();
    }

    QQmlListProperty<QObject> contentData();
    QQmlListProperty<QQuickItem> contentChildren();

    void init();
    void closeOrReject();
    bool tryClose(const QPointF &pos, QQuickPopup::ClosePolicy flags);

    bool contains(const QPointF &scenePos) const;

#if QT_CONFIG(quicktemplates2_multitouch)
    virtual bool acceptTouch(const QTouchEvent::TouchPoint &point);
#endif
    virtual bool blockInput(QQuickItem *item, const QPointF &point) const;

    virtual bool handlePress(QQuickItem* item, const QPointF &point, ulong timestamp);
    virtual bool handleMove(QQuickItem* item, const QPointF &point, ulong timestamp);
    virtual bool handleRelease(QQuickItem* item, const QPointF &point, ulong timestamp);
    virtual void handleUngrab();

    bool handleMouseEvent(QQuickItem *item, QMouseEvent *event);
#if QT_CONFIG(quicktemplates2_multitouch)
    bool handleTouchEvent(QQuickItem *item, QTouchEvent *event);
#endif

    void reposition();

    void createOverlay();
    void destroyOverlay();
    void toggleOverlay();
    virtual void showOverlay();
    virtual void hideOverlay();
    virtual void resizeOverlay();

    virtual bool prepareEnterTransition();
    virtual bool prepareExitTransition();
    virtual void finalizeEnterTransition();
    virtual void finalizeExitTransition();

    QMarginsF getMargins() const;

    void setTopMargin(qreal value, bool reset = false);
    void setLeftMargin(qreal value, bool reset = false);
    void setRightMargin(qreal value, bool reset = false);
    void setBottomMargin(qreal value, bool reset = false);

    QQuickPopupAnchors *getAnchors();
    virtual QQuickPopupPositioner *getPositioner();

    void setWindow(QQuickWindow *window);
    void itemDestroyed(QQuickItem *item) override;

    enum TransitionState {
        NoTransition, EnterTransition, ExitTransition
    };

    static const QQuickPopup::ClosePolicy DefaultClosePolicy;

    bool focus = false;
    bool modal = false;
    bool dim = false;
    bool hasDim = false;
    bool visible = false;
    bool complete = true;
    bool positioning = false;
    bool hasWidth = false;
    bool hasHeight = false;
    bool hasTopMargin = false;
    bool hasLeftMargin = false;
    bool hasRightMargin = false;
    bool hasBottomMargin = false;
    bool allowVerticalFlip = false;
    bool allowHorizontalFlip = false;
    bool allowVerticalMove = true;
    bool allowHorizontalMove = true;
    bool allowVerticalResize = true;
    bool allowHorizontalResize = true;
    bool hadActiveFocusBeforeExitTransition = false;
    bool interactive = true;
    bool hasClosePolicy = false;
    int touchId = -1;
    qreal x = 0;
    qreal y = 0;
    qreal effectiveX = 0;
    qreal effectiveY = 0;
    qreal margins = -1;
    qreal topMargin = 0;
    qreal leftMargin = 0;
    qreal rightMargin = 0;
    qreal bottomMargin = 0;
    QPointF pressPoint;
    TransitionState transitionState = NoTransition;
    QQuickPopup::ClosePolicy closePolicy = DefaultClosePolicy;
    QQuickItem *parentItem = nullptr;
    QQuickItem *dimmer = nullptr;
    QPointer<QQuickWindow> window;
    QQuickTransition *enter = nullptr;
    QQuickTransition *exit = nullptr;
    QQuickPopupItem *popupItem = nullptr;
    QQuickPopupPositioner *positioner = nullptr;
    QList<QQuickStateAction> enterActions;
    QList<QQuickStateAction> exitActions;
    QQuickPopupTransitionManager transitionManager;
    QQuickPopupAnchors *anchors = nullptr;
    qreal prevOpacity = 0;
    qreal prevScale = 0;

    friend class QQuickPopupTransitionManager;
};

QT_END_NAMESPACE

#endif // QQUICKPOPUP_P_P_H
