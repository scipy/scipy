/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Designer of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:GPL-EXCEPT$
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
** General Public License version 3 as published by the Free Software
** Foundation with exceptions as appearing in the file LICENSE.GPL3-EXCEPT
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of Qt Designer.  This header
// file may change from version to version without notice, or even be removed.
//
// We mean it.
//

#ifndef QDESIGNER_TOOLBAR_H
#define QDESIGNER_TOOLBAR_H

#include "shared_global_p.h"

#include <QtWidgets/qaction.h>
#include <QtWidgets/qtoolbutton.h>

#include <QtCore/qlist.h>
#include <QtCore/qpoint.h>

QT_BEGIN_NAMESPACE

class QDesignerFormWindowInterface;
class QToolBar;
class QRect;
class QAction;

namespace qdesigner_internal {

class PromotionTaskMenu;

// Special event filter for  tool bars in designer.
// Handles drag and drop to and from. Ensures that each
// child widget is  WA_TransparentForMouseEvents to enable  drag and drop.

class QDESIGNER_SHARED_EXPORT ToolBarEventFilter : public QObject {
    Q_OBJECT

public:
    static void install(QToolBar *tb);

    // Find action by position. Note that QToolBar::actionAt() will
    // not work as designer sets WA_TransparentForMouseEvents on its tool bar buttons
    // to be able to drag them. This function will return the dummy
    // sentinel action when applied to tool bars created by designer if the position matches.
    static QAction *actionAt(const QToolBar *tb, const QPoint &pos);

    static bool withinHandleArea(const QToolBar *tb, const QPoint &pos);
    static QRect handleArea(const QToolBar *tb);
    static QRect freeArea(const QToolBar *tb);

    // Utility to create an action
    static QAction *createAction(QDesignerFormWindowInterface *fw, const QString &objectName, bool separator);

    bool eventFilter (QObject *watched, QEvent *event) override;

    // Helper for task menu extension
    QList<QAction *> contextMenuActions(const QPoint &globalPos = QPoint(-1, -1));

    static ToolBarEventFilter *eventFilterOf(const QToolBar *tb);

private slots:
    void slotRemoveSelectedAction();
    void slotRemoveToolBar();
    void slotInsertSeparator();

private:
    explicit ToolBarEventFilter(QToolBar *tb);

    bool handleContextMenuEvent(QContextMenuEvent * event);
    bool handleDragEnterMoveEvent(QDragMoveEvent *event);
    bool handleDragLeaveEvent(QDragLeaveEvent *);
    bool handleDropEvent(QDropEvent *event);
    bool handleMousePressEvent(QMouseEvent *event);
    bool handleMouseReleaseEvent(QMouseEvent *event);
    bool handleMouseMoveEvent(QMouseEvent *event);

    QDesignerFormWindowInterface *formWindow() const;
    void adjustDragIndicator(const QPoint &pos);
    void hideDragIndicator();
    void startDrag(const QPoint &pos, Qt::KeyboardModifiers modifiers);
    bool withinHandleArea(const QPoint &pos) const;

    QToolBar *m_toolBar;
    PromotionTaskMenu *m_promotionTaskMenu;
    QPoint m_startPosition;
};
} // namespace qdesigner_internal

QT_END_NAMESPACE

#endif // QDESIGNER_TOOLBAR_H
