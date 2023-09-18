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

#ifndef ABSTRACTNEWFORMWIDGET_H
#define ABSTRACTNEWFORMWIDGET_H

#include <QtDesigner/sdk_global.h>

#include <QtWidgets/qwidget.h>

QT_BEGIN_NAMESPACE

class QDesignerFormEditorInterface;

class QDESIGNER_SDK_EXPORT QDesignerNewFormWidgetInterface : public QWidget
{
    Q_DISABLE_COPY(QDesignerNewFormWidgetInterface)
    Q_OBJECT
public:
    explicit QDesignerNewFormWidgetInterface(QWidget *parent = nullptr);
    virtual ~QDesignerNewFormWidgetInterface();

    virtual bool hasCurrentTemplate() const = 0;
    virtual QString currentTemplate(QString *errorMessage = nullptr) = 0;

    static QDesignerNewFormWidgetInterface *createNewFormWidget(QDesignerFormEditorInterface *core, QWidget *parent = nullptr);

Q_SIGNALS:
    void templateActivated();
    void currentTemplateChanged(bool templateSelected);
};

QT_END_NAMESPACE

#endif // ABSTRACTNEWFORMWIDGET_H
