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

#ifndef QTRESOURCEEDITOR_H
#define QTRESOURCEEDITOR_H

#include <QtCore/qscopedpointer.h>
#include <QtWidgets/qdialog.h>

QT_BEGIN_NAMESPACE

class QtResourceModel;
class QDesignerDialogGuiInterface;
class QDesignerFormEditorInterface;

class QtResourceEditorDialog : public QDialog
{
    Q_OBJECT
public:
    QtResourceModel *model() const;
    void setResourceModel(QtResourceModel *model);

    QString selectedResource() const;

    static QString editResources(QDesignerFormEditorInterface *core, QtResourceModel *model,
                                 QDesignerDialogGuiInterface *dlgGui, QWidget *parent = nullptr);

    // Helper to display a message box with rcc logs in case of errors.
    static void displayResourceFailures(const QString &logOutput, QDesignerDialogGuiInterface *dlgGui, QWidget *parent = nullptr);

public slots:
    void accept() override;

private:
    QtResourceEditorDialog(QDesignerFormEditorInterface *core, QDesignerDialogGuiInterface *dlgGui, QWidget *parent = nullptr);
    ~QtResourceEditorDialog() override;

    QScopedPointer<class QtResourceEditorDialogPrivate> d_ptr;
    Q_DECLARE_PRIVATE(QtResourceEditorDialog)
    Q_DISABLE_COPY_MOVE(QtResourceEditorDialog)

    Q_PRIVATE_SLOT(d_func(), void slotQrcFileInserted(QtQrcFile *))
    Q_PRIVATE_SLOT(d_func(), void slotQrcFileMoved(QtQrcFile *))
    Q_PRIVATE_SLOT(d_func(), void slotQrcFileRemoved(QtQrcFile *))
    Q_PRIVATE_SLOT(d_func(), void slotResourcePrefixInserted(QtResourcePrefix *))
    Q_PRIVATE_SLOT(d_func(), void slotResourcePrefixMoved(QtResourcePrefix *))
    Q_PRIVATE_SLOT(d_func(), void slotResourcePrefixChanged(QtResourcePrefix *))
    Q_PRIVATE_SLOT(d_func(), void slotResourceLanguageChanged(QtResourcePrefix *))
    Q_PRIVATE_SLOT(d_func(), void slotResourcePrefixRemoved(QtResourcePrefix *))
    Q_PRIVATE_SLOT(d_func(), void slotResourceFileInserted(QtResourceFile *))
    Q_PRIVATE_SLOT(d_func(), void slotResourceFileMoved(QtResourceFile *))
    Q_PRIVATE_SLOT(d_func(), void slotResourceAliasChanged(QtResourceFile *))
    Q_PRIVATE_SLOT(d_func(), void slotResourceFileRemoved(QtResourceFile *))

    Q_PRIVATE_SLOT(d_func(), void slotCurrentQrcFileChanged(QListWidgetItem *))
    Q_PRIVATE_SLOT(d_func(), void slotCurrentTreeViewItemChanged(const QModelIndex &))
    Q_PRIVATE_SLOT(d_func(), void slotListWidgetContextMenuRequested(const QPoint &))
    Q_PRIVATE_SLOT(d_func(), void slotTreeViewContextMenuRequested(const QPoint &))
    Q_PRIVATE_SLOT(d_func(), void slotTreeViewItemChanged(QStandardItem *))

    Q_PRIVATE_SLOT(d_func(), void slotNewQrcFile())
    Q_PRIVATE_SLOT(d_func(), void slotImportQrcFile())
    Q_PRIVATE_SLOT(d_func(), void slotRemoveQrcFile())
    Q_PRIVATE_SLOT(d_func(), void slotMoveUpQrcFile())
    Q_PRIVATE_SLOT(d_func(), void slotMoveDownQrcFile())

    Q_PRIVATE_SLOT(d_func(), void slotNewPrefix())
    Q_PRIVATE_SLOT(d_func(), void slotAddFiles())
    Q_PRIVATE_SLOT(d_func(), void slotChangePrefix())
    Q_PRIVATE_SLOT(d_func(), void slotChangeLanguage())
    Q_PRIVATE_SLOT(d_func(), void slotChangeAlias())
    Q_PRIVATE_SLOT(d_func(), void slotClonePrefix())
    Q_PRIVATE_SLOT(d_func(), void slotRemove())
    Q_PRIVATE_SLOT(d_func(), void slotMoveUp())
    Q_PRIVATE_SLOT(d_func(), void slotMoveDown())
};

QT_END_NAMESPACE

#endif

