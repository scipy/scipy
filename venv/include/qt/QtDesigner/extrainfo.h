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

#ifndef EXTRAINFO_H
#define EXTRAINFO_H

#include <QtDesigner/sdk_global.h>
#include <QtDesigner/extension.h>

QT_BEGIN_NAMESPACE

class DomWidget;
class DomUI;
class QWidget;

class QDesignerFormEditorInterface;

class QDESIGNER_SDK_EXPORT QDesignerExtraInfoExtension
{
public:
    virtual ~QDesignerExtraInfoExtension() {}

    virtual QDesignerFormEditorInterface *core() const = 0;
    virtual QWidget *widget() const = 0;

    virtual bool saveUiExtraInfo(DomUI *ui) = 0;
    virtual bool loadUiExtraInfo(DomUI *ui) = 0;

    virtual bool saveWidgetExtraInfo(DomWidget *ui_widget) = 0;
    virtual bool loadWidgetExtraInfo(DomWidget *ui_widget) = 0;

    QString workingDirectory() const;
    void setWorkingDirectory(const QString &workingDirectory);

private:
    QString m_workingDirectory;
};
Q_DECLARE_EXTENSION_INTERFACE(QDesignerExtraInfoExtension, "org.qt-project.Qt.Designer.ExtraInfo.2")

QT_END_NAMESPACE

#endif // EXTRAINFO_H
