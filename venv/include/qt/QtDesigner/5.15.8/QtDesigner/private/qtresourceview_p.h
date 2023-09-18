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

#ifndef QTRESOURCEVIEW_H
#define QTRESOURCEVIEW_H

#include "shared_global_p.h"
#include <QtWidgets/qwidget.h>
#include <QtWidgets/qdialog.h>

QT_BEGIN_NAMESPACE

class QtResourceModel;
class QtResourceSet;
class QDesignerFormEditorInterface;
class QMimeData;

class QDESIGNER_SHARED_EXPORT QtResourceView : public QWidget
{
    Q_OBJECT
public:
    explicit QtResourceView(QDesignerFormEditorInterface *core, QWidget *parent = nullptr);
    ~QtResourceView();

    void setDragEnabled(bool dragEnabled);
    bool dragEnabled() const;

    QtResourceModel *model() const;
    void setResourceModel(QtResourceModel *model);

    QString selectedResource() const;
    void selectResource(const QString &resource);

    QString settingsKey() const;
    void setSettingsKey(const QString &key);

    bool isResourceEditingEnabled() const;
    void setResourceEditingEnabled(bool enable);

    // Helpers for handling the drag originating in QtResourceView (Xml/text)
    enum ResourceType { ResourceImage, ResourceStyleSheet, ResourceOther };
    static QString encodeMimeData(ResourceType resourceType, const QString &path);

    static bool decodeMimeData(const QMimeData *md, ResourceType *t = nullptr, QString *file = nullptr);
    static bool decodeMimeData(const QString &text, ResourceType *t = nullptr, QString *file = nullptr);

signals:
    void resourceSelected(const QString &resource);
    void resourceActivated(const QString &resource);

protected:
    bool event(QEvent *event) override;

private:

    QScopedPointer<class QtResourceViewPrivate> d_ptr;
    Q_DECLARE_PRIVATE(QtResourceView)
    Q_DISABLE_COPY_MOVE(QtResourceView)
    Q_PRIVATE_SLOT(d_func(), void slotResourceSetActivated(QtResourceSet *))
    Q_PRIVATE_SLOT(d_func(), void slotCurrentPathChanged(QTreeWidgetItem *))
    Q_PRIVATE_SLOT(d_func(), void slotCurrentResourceChanged(QListWidgetItem *))
    Q_PRIVATE_SLOT(d_func(), void slotResourceActivated(QListWidgetItem *))
    Q_PRIVATE_SLOT(d_func(), void slotEditResources())
    Q_PRIVATE_SLOT(d_func(), void slotReloadResources())
#if QT_CONFIG(clipboard)
    Q_PRIVATE_SLOT(d_func(), void slotCopyResourcePath())
#endif
    Q_PRIVATE_SLOT(d_func(), void slotListWidgetContextMenuRequested(const QPoint &pos))
    Q_PRIVATE_SLOT(d_func(), void slotFilterChanged(const QString &pattern))
};

class QDESIGNER_SHARED_EXPORT  QtResourceViewDialog : public QDialog
{
    Q_OBJECT
public:
    explicit QtResourceViewDialog(QDesignerFormEditorInterface *core, QWidget *parent = nullptr);
    ~QtResourceViewDialog() override;

    QString selectedResource() const;
    void selectResource(const QString &path);

    bool isResourceEditingEnabled() const;
    void setResourceEditingEnabled(bool enable);

private:
    QScopedPointer<class QtResourceViewDialogPrivate> d_ptr;
    Q_DECLARE_PRIVATE(QtResourceViewDialog)
    Q_DISABLE_COPY_MOVE(QtResourceViewDialog)
    Q_PRIVATE_SLOT(d_func(), void slotResourceSelected(const QString &))
};

QT_END_NAMESPACE

#endif
