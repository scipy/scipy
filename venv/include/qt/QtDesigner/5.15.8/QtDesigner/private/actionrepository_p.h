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

#ifndef ACTIONREPOSITORY_H
#define ACTIONREPOSITORY_H

#include "shared_global_p.h"
#include <QtCore/qmimedata.h>
#include <QtGui/qstandarditemmodel.h>
#include <QtWidgets/qtreeview.h>
#include <QtWidgets/qlistview.h>
#include <QtWidgets/qstackedwidget.h>
#include <QtGui/qicon.h>

QT_BEGIN_NAMESPACE

class QPixmap;

class QDesignerFormEditorInterface;
class QDesignerPropertySheetExtension;

namespace qdesigner_internal {

class PropertySheetKeySequenceValue;

// Shared model of actions, to be used for several views (detailed/icon view).
class QDESIGNER_SHARED_EXPORT ActionModel: public QStandardItemModel
{
    Q_OBJECT
public:
    enum Columns { NameColumn, UsedColumn, TextColumn, ShortCutColumn, CheckedColumn, ToolTipColumn, NumColumns };
    enum   { ActionRole = Qt::UserRole + 1000 };

    explicit ActionModel(QWidget *parent = nullptr);
    void initialize(QDesignerFormEditorInterface *core) { m_core = core; }

    void clearActions();
    QModelIndex addAction(QAction *a);
    // remove row
    void remove(int row);
    // update the row from the underlying action
    void update(int row);

    // return row of action or -1.
    int findAction(QAction *) const;

    QString actionName(int row) const;
    QAction *actionAt(const QModelIndex &index) const;

    QMimeData *mimeData(const QModelIndexList &indexes) const override;
    QStringList mimeTypes() const override;
    bool dropMimeData(const QMimeData *data, Qt::DropAction action, int row, int column, const QModelIndex &parent) override;

    // Find the associated menus and toolbars, ignore toolbuttons
    static QWidgetList associatedWidgets(const QAction *action);

    // Retrieve shortcut via property sheet as it is a fake property
    static PropertySheetKeySequenceValue actionShortCut(QDesignerFormEditorInterface *core, QAction *action);
    static PropertySheetKeySequenceValue actionShortCut(const QDesignerPropertySheetExtension *ps);

signals:
    void resourceImageDropped(const QString &path, QAction *action);

private:
    using QStandardItemList = QList<QStandardItem *>;

    void initializeHeaders();
    static void setItems(QDesignerFormEditorInterface *core, QAction *a,
                         const QIcon &defaultIcon,
                         QStandardItemList &sl);

    const QIcon m_emptyIcon;

    QDesignerFormEditorInterface *m_core = nullptr;
};

// Internal class that provides the detailed view of actions.
class  ActionTreeView: public QTreeView
{
    Q_OBJECT
public:
    explicit ActionTreeView(ActionModel *model, QWidget *parent = nullptr);
    QAction *currentAction() const;

public slots:
    void filter(const QString &text);

signals:
    void actionContextMenuRequested(QContextMenuEvent *event, QAction *);
    void currentActionChanged(QAction *action);
    void actionActivated(QAction *action, int column);

protected slots:
    void currentChanged(const QModelIndex &current, const QModelIndex &previous) override;

protected:
    void dragEnterEvent(QDragEnterEvent *event) override;
    void dragMoveEvent(QDragMoveEvent *event) override;
    void dropEvent(QDropEvent *event) override;
    void focusInEvent(QFocusEvent *event) override;
    void contextMenuEvent(QContextMenuEvent *event) override;
    void startDrag(Qt::DropActions supportedActions) override;

private slots:
    void slotActivated(const QModelIndex &);

private:
    ActionModel *m_model;
};

// Internal class that provides the icon view of actions.
class ActionListView: public QListView
{
    Q_OBJECT
public:
    explicit ActionListView(ActionModel *model, QWidget *parent = nullptr);
    QAction *currentAction() const;

public slots:
    void filter(const QString &text);

signals:
    void actionContextMenuRequested(QContextMenuEvent *event, QAction *);
    void currentActionChanged(QAction *action);
    void actionActivated(QAction *action);

protected slots:
    void currentChanged(const QModelIndex &current, const QModelIndex &previous) override;

protected:
    void dragEnterEvent(QDragEnterEvent *event) override;
    void dragMoveEvent(QDragMoveEvent *event) override;
    void dropEvent(QDropEvent *event) override;
    void focusInEvent(QFocusEvent *event) override;
    void contextMenuEvent(QContextMenuEvent *event) override;
    void startDrag(Qt::DropActions supportedActions) override;

private slots:
    void slotActivated(const QModelIndex &);

private:
    ActionModel *m_model;
};

// Action View that can be switched between detailed and icon view
// using a  QStackedWidget of  ActionListView / ActionTreeView
// that share the item model and the selection model.

class ActionView : public  QStackedWidget {
    Q_OBJECT
public:
    // Separate initialize() function takes core argument to make this
    // thing usable as promoted widget.
    explicit ActionView(QWidget *parent = nullptr);
    void initialize(QDesignerFormEditorInterface *core) { m_model->initialize(core); }

    // View mode
    enum { DetailedView, IconView };
    int viewMode() const;
    void setViewMode(int lm);

    void setSelectionMode(QAbstractItemView::SelectionMode sm);
    QAbstractItemView::SelectionMode selectionMode() const;

    ActionModel *model() const { return m_model; }

    QAction *currentAction() const;
    void setCurrentIndex(const QModelIndex &index);

    using ActionList = QList<QAction *>;
    ActionList selectedActions() const;
    QItemSelection selection() const;

public slots:
    void filter(const QString &text);
    void selectAll();
    void clearSelection();

signals:
    void contextMenuRequested(QContextMenuEvent *event, QAction *);
    void currentChanged(QAction *action);
    void activated(QAction *action, int column);
    void selectionChanged(const QItemSelection& selected, const QItemSelection& deselected);
    void resourceImageDropped(const QString &data, QAction *action);

private slots:
    void slotCurrentChanged(QAction *action);

private:
    ActionModel *m_model;
    ActionTreeView *m_actionTreeView;
    ActionListView *m_actionListView;
};

class QDESIGNER_SHARED_EXPORT ActionRepositoryMimeData: public QMimeData
{
    Q_OBJECT
public:
    using ActionList = QList<QAction *>;

    ActionRepositoryMimeData(const ActionList &, Qt::DropAction dropAction);
    ActionRepositoryMimeData(QAction *, Qt::DropAction dropAction);

    const ActionList &actionList() const { return m_actionList; }
    QStringList formats() const override;

    static QPixmap actionDragPixmap(const QAction *action);

    // Utility to accept with right action
    void accept(QDragMoveEvent *event) const;
private:
    const Qt::DropAction m_dropAction;
    ActionList m_actionList;
};

} // namespace qdesigner_internal

QT_END_NAMESPACE

#endif // ACTIONREPOSITORY_H
