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

#ifndef QDESIGNER_MENUBAR_H
#define QDESIGNER_MENUBAR_H

#include "shared_global_p.h"

#include <QtWidgets/qaction.h>
#include <QtWidgets/qmenubar.h>

#include <QtCore/qpointer.h>
#include <QtCore/qmimedata.h>

QT_BEGIN_NAMESPACE

class QDesignerFormWindowInterface;
class QDesignerActionProviderExtension;

class QLineEdit;
class QMimeData;

namespace qdesigner_internal {
class PromotionTaskMenu;

class SpecialMenuAction: public QAction
{
    Q_OBJECT
public:
    SpecialMenuAction(QObject *parent = nullptr);
    ~SpecialMenuAction() override;
};

} // namespace qdesigner_internal

class QDESIGNER_SHARED_EXPORT QDesignerMenuBar: public QMenuBar
{
    Q_OBJECT
public:
    QDesignerMenuBar(QWidget *parent = nullptr);
    ~QDesignerMenuBar() override;

    bool eventFilter(QObject *object, QEvent *event) override;

    QDesignerFormWindowInterface *formWindow() const;
    QDesignerActionProviderExtension *actionProvider();

    void adjustSpecialActions();
    bool dragging() const;

    void moveLeft(bool ctrl = false);
    void moveRight(bool ctrl = false);
    void moveUp();
    void moveDown();

    // Helpers for MenuTaskMenu/MenuBarTaskMenu extensions
    QList<QAction *> contextMenuActions();
    void deleteMenuAction(QAction *action);

private slots:
    void deleteMenu();
    void slotRemoveMenuBar();

protected:
    void actionEvent(QActionEvent *event) override;
    void dragEnterEvent(QDragEnterEvent *event) override;
    void dragMoveEvent(QDragMoveEvent *event) override;
    void dragLeaveEvent(QDragLeaveEvent *event) override;
    void dropEvent(QDropEvent *event) override;
    void paintEvent(QPaintEvent *event) override;
    void focusOutEvent(QFocusEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;
    void keyReleaseEvent(QKeyEvent *event) override;

    bool handleEvent(QWidget *widget, QEvent *event);
    bool handleMouseDoubleClickEvent(QWidget *widget, QMouseEvent *event);
    bool handleMousePressEvent(QWidget *widget, QMouseEvent *event);
    bool handleMouseReleaseEvent(QWidget *widget, QMouseEvent *event);
    bool handleMouseMoveEvent(QWidget *widget, QMouseEvent *event);
    bool handleContextMenuEvent(QWidget *widget, QContextMenuEvent *event);
    bool handleKeyPressEvent(QWidget *widget, QKeyEvent *event);

    void startDrag(const QPoint &pos);

    enum ActionDragCheck { NoActionDrag, ActionDragOnSubMenu, AcceptActionDrag };
    ActionDragCheck checkAction(QAction *action) const;

    void adjustIndicator(const QPoint &pos);
    int findAction(const QPoint &pos) const;

    QAction *currentAction() const;
    int realActionCount() const;

    enum LeaveEditMode {
        Default = 0,
        ForceAccept
    };

    void enterEditMode();
    void leaveEditMode(LeaveEditMode mode);
    void showLineEdit();

    void showMenu(int index = -1);
    void hideMenu(int index = -1);

    QAction *safeActionAt(int index) const;

    bool swapActions(int a, int b);

private:
    void updateCurrentAction(bool selectAction);
    void movePrevious(bool ctrl);
    void moveNext(bool ctrl);

    QAction *m_addMenu;
    QPointer<QMenu> m_activeMenu;
    QPoint m_startPosition;
    int m_currentIndex = 0;
    QLineEdit *m_editor;
    bool m_dragging = false;
    int m_lastMenuActionIndex = -1;
    QPointer<QWidget> m_lastFocusWidget;
    qdesigner_internal::PromotionTaskMenu* m_promotionTaskMenu;
};

QT_END_NAMESPACE

#endif // QDESIGNER_MENUBAR_H
