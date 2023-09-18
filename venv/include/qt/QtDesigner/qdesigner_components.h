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

#ifndef QDESIGNER_COMPONENTS_H
#define QDESIGNER_COMPONENTS_H

#include <QtDesigner/qdesigner_components_global.h>
#include <QtCore/qglobal.h>

QT_BEGIN_NAMESPACE

class QObject;
class QWidget;

class QDesignerFormEditorInterface;
class QDesignerWidgetBoxInterface;
class QDesignerPropertyEditorInterface;
class QDesignerObjectInspectorInterface;
class QDesignerActionEditorInterface;

class QDESIGNER_COMPONENTS_EXPORT QDesignerComponents
{
public:
    static void initializeResources();
    static void initializePlugins(QDesignerFormEditorInterface *core);

    static QDesignerFormEditorInterface *createFormEditor(QObject *parent);
    static QDesignerWidgetBoxInterface *createWidgetBox(QDesignerFormEditorInterface *core, QWidget *parent);
    static QDesignerPropertyEditorInterface *createPropertyEditor(QDesignerFormEditorInterface *core, QWidget *parent);
    static QDesignerObjectInspectorInterface *createObjectInspector(QDesignerFormEditorInterface *core, QWidget *parent);
    static QDesignerActionEditorInterface *createActionEditor(QDesignerFormEditorInterface *core, QWidget *parent);

    static QObject *createTaskMenu(QDesignerFormEditorInterface *core, QObject *parent);
    static QWidget *createResourceEditor(QDesignerFormEditorInterface *core, QWidget *parent);
    static QWidget *createSignalSlotEditor(QDesignerFormEditorInterface *core, QWidget *parent);
};

QT_END_NAMESPACE

#endif // QDESIGNER_COMPONENTS_H
