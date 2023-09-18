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

#ifndef QDESIGNER_MENU_H
#define QDESIGNER_MENU_H

#include "shared_global_p.h"

#include <QtWidgets/qaction.h>
#include <QtWidgets/qmenu.h>
#include <QtGui/qpixmap.h>
#include <QtCore/qhash.h>

QT_BEGIN_NAMESPACE

class QTimer;
class QLineEdit;

class QDesignerFormWindowInterface;
class QDesignerActionProviderExtension;
class QDesignerMenu;
class QDesignerMenuBar;
class QPainter;
class QMimeData;

namespace qdesigner_internal {
    class CreateSubmenuCommand;
    class ActionInsertionCommand;
}

class QDESIGNER_SHARED_EXPORT QDesignerMenu: public QMenu
{
    Q_OBJECT
public:
    QDesignerMenu(QWidget *parent = nullptr);
    ~QDesignerMenu() override;

    bool eventFilter(QObject *object, QEvent *event) override;

    QDesignerFormWindowInterface *formWindow() const;
    QDesignerActionProviderExtension *actionProvider();

    QDesignerMenu *parentMenu() const;
    QDesignerMenuBar *parentMenuBar() const;

    void setVisible(bool visible) override;

    void adjustSpecialActions();

    void createRealMenuAction(QAction *action);
    void removeRealMenu(QAction *action);

    static void drawSelection(QPainter *p, const QRect &r);

    bool dragging() const;

    void closeMenuChain();

    void moveLeft();
    void moveRight();
    void moveUp(bool ctrl);
    void moveDown(bool ctrl);

    // Helper for MenuTaskMenu extension
    void deleteAction(QAction *a);

private slots:
    void slotAddSeparator();
    void slotRemoveSelectedAction();
    void slotShowSubMenuNow();
    void slotDeactivateNow();
    void slotAdjustSizeNow();

protected:
    void actionEvent(QActionEvent *event) override;
    void dragEnterEvent(QDragEnterEvent *event) override;
    void dragMoveEvent(QDragMoveEvent *event) override;
    void dragLeaveEvent(QDragLeaveEvent *event) override;
    void dropEvent(QDropEvent *event) override;
    void paintEvent(QPaintEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;
    void keyReleaseEvent(QKeyEvent *event) override;
    void showEvent(QShowEvent *event) override;

    bool handleEvent(QWidget *widget, QEvent *event);
    bool handleMouseDoubleClickEvent(QWidget *widget, QMouseEvent *event);
    bool handleMousePressEvent(QWidget *widget, QMouseEvent *event);
    bool handleMouseReleaseEvent(QWidget *widget, QMouseEvent *event);
    bool handleMouseMoveEvent(QWidget *widget, QMouseEvent *event);
    bool handleContextMenuEvent(QWidget *widget, QContextMenuEvent *event);
    bool handleKeyPressEvent(QWidget *widget, QKeyEvent *event);

    void startDrag(const QPoint &pos, Qt::KeyboardModifiers modifiers);

    void adjustIndicator(const QPoint &pos);
    int findAction(const QPoint &pos) const;

    QAction *currentAction() const;
    int realActionCount() const;
    enum ActionDragCheck { NoActionDrag, ActionDragOnSubMenu, AcceptActionDrag };
    ActionDragCheck checkAction(QAction *action) const;

    void showSubMenu(QAction *action);

    enum LeaveEditMode {
        Default = 0,
        ForceAccept
    };

    void enterEditMode();
    void leaveEditMode(LeaveEditMode mode);
    void showLineEdit();

    QAction *createAction(const QString &text, bool separator = false);
    QDesignerMenu *findOrCreateSubMenu(QAction *action);

    QAction *safeActionAt(int index) const;
    QAction *safeMenuAction(QDesignerMenu *menu) const;
    bool swap(int a, int b);

    void hideSubMenu();
    void deleteAction();
    void deactivateMenu();

    bool canCreateSubMenu(QAction *action) const;
    QDesignerMenu *findRootMenu() const;
    QDesignerMenu *findActivatedMenu() const;

    QRect subMenuPixmapRect(QAction *action) const;
    bool hasSubMenuPixmap(QAction *action) const;

    void selectCurrentAction();

private:
    bool hideSubMenuOnCursorKey();
    bool showSubMenuOnCursorKey();
    const QPixmap m_subMenuPixmap;

    QPoint m_startPosition;
    int m_currentIndex = 0;
    QAction *m_addItem;
    QAction *m_addSeparator;
    QHash<QAction*, QDesignerMenu*> m_subMenus;
    QTimer *m_showSubMenuTimer;
    QTimer *m_deactivateWindowTimer;
    QTimer *m_adjustSizeTimer;
    QLineEdit *m_editor;
    bool m_dragging = false;
    int m_lastSubMenuIndex = -1;

    friend class qdesigner_internal::CreateSubmenuCommand;
    friend class qdesigner_internal::ActionInsertionCommand;
};

QT_END_NAMESPACE

#endif // QDESIGNER_MENU_H
