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

#ifndef PLUGINDIALOG_H
#define PLUGINDIALOG_H

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

#include "ui_plugindialog.h"

QT_BEGIN_NAMESPACE

class QDesignerFormEditorInterface;

namespace qdesigner_internal {

class PluginDialog : public QDialog
{
    Q_OBJECT
public:
    explicit PluginDialog(QDesignerFormEditorInterface *core, QWidget *parent = nullptr);

private slots:
    void updateCustomWidgetPlugins();

private:
    void populateTreeWidget();
    QTreeWidgetItem* setTopLevelItem(const QString &itemName);
    QTreeWidgetItem* setPluginItem(QTreeWidgetItem *topLevelItem,
                                   const QString &itemName, const QFont &font);
    void setItem(QTreeWidgetItem *pluginItem, const QString &name,
                 const QString &toolTip, const QString &whatsThis, const QIcon &icon);

    QDesignerFormEditorInterface *m_core;
    Ui::PluginDialog ui;
    QIcon interfaceIcon;
    QIcon featureIcon;
};

}

QT_END_NAMESPACE

#endif
