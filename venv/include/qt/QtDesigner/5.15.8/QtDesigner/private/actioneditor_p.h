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

#ifndef ACTIONEDITOR_H
#define ACTIONEDITOR_H

#include "shared_global_p.h"
#include "shared_enums_p.h"
#include <QtDesigner/abstractactioneditor.h>

#include <QtCore/qpointer.h>

QT_BEGIN_NAMESPACE

class QDesignerPropertyEditorInterface;
class QDesignerSettingsInterface;
class QMenu;
class QActionGroup;
class QItemSelection;
class QListWidget;
class QPushButton;
class QLineEdit;
class QToolButton;

namespace qdesigner_internal {

class ActionView;
class ResourceMimeData;

class QDESIGNER_SHARED_EXPORT ActionEditor: public QDesignerActionEditorInterface
{
    Q_OBJECT
public:
    explicit ActionEditor(QDesignerFormEditorInterface *core, QWidget *parent = nullptr,
                          Qt::WindowFlags flags = {});
    ~ActionEditor() override;

    QDesignerFormWindowInterface *formWindow() const;
    void setFormWindow(QDesignerFormWindowInterface *formWindow) override;

    QDesignerFormEditorInterface *core() const override;

    QAction *actionNew() const;
    QAction *actionDelete() const;

    QString filter() const;

    void manageAction(QAction *action) override;
    void unmanageAction(QAction *action) override;

    static ObjectNamingMode objectNamingMode() { return m_objectNamingMode; }
    static void setObjectNamingMode(ObjectNamingMode n) { m_objectNamingMode = n; }

    static QString actionTextToName(const QString &text,
                                    const QString &prefix = QLatin1String("action"));

    // Utility to create a configure button with menu for usage on toolbars
    static QToolButton *createConfigureMenuButton(const QString &t, QMenu **ptrToMenu);

public slots:
    void setFilter(const QString &filter);
    void mainContainerChanged();

private slots:
    void slotCurrentItemChanged(QAction *item);
    void slotSelectionChanged(const QItemSelection& selected, const QItemSelection& deselected);
    void editAction(QAction *item, int column = -1);
    void editCurrentAction();
    void navigateToSlotCurrentAction();
    void slotActionChanged();
    void slotNewAction();
    void slotDelete();
    void resourceImageDropped(const QString &path, QAction *action);
    void slotContextMenuRequested(QContextMenuEvent *, QAction *);
    void slotViewMode(QAction *a);
    void slotSelectAssociatedWidget(QWidget *w);
#if QT_CONFIG(clipboard)
    void slotCopy();
    void slotCut();
    void slotPaste();
#endif

signals:
    void itemActivated(QAction *item, int column);
    // Context menu for item or global menu if item == 0.
    void contextMenuRequested(QMenu *menu, QAction *item);

private:
    using ActionList = QList<QAction *>;
    void deleteActions(QDesignerFormWindowInterface *formWindow, const ActionList &);
#if QT_CONFIG(clipboard)
    void copyActions(QDesignerFormWindowInterface *formWindow, const ActionList &);
#endif

    void restoreSettings();
    void saveSettings();

    void updateViewModeActions();

    static ObjectNamingMode m_objectNamingMode;

    QDesignerFormEditorInterface *m_core;
    QPointer<QDesignerFormWindowInterface> m_formWindow;
    QListWidget *m_actionGroups;

    ActionView *m_actionView;

    QAction *m_actionNew;
    QAction *m_actionEdit;
    QAction *m_actionNavigateToSlot;
#if QT_CONFIG(clipboard)
    QAction *m_actionCopy;
    QAction *m_actionCut;
    QAction *m_actionPaste;
#endif
    QAction *m_actionSelectAll;
    QAction *m_actionDelete;

    QActionGroup *m_viewModeGroup;
    QAction *m_iconViewAction;
    QAction *m_listViewAction;

    QString m_filter;
    QWidget *m_filterWidget;
};

} // namespace qdesigner_internal

QT_END_NAMESPACE

#endif // ACTIONEDITOR_H
