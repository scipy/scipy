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

#ifndef ABSTRACTOBJECTINSPECTOR_H
#define ABSTRACTOBJECTINSPECTOR_H

#include <QtDesigner/sdk_global.h>

#include <QtWidgets/qwidget.h>

QT_BEGIN_NAMESPACE

class QDesignerFormEditorInterface;
class QDesignerFormWindowInterface;

class QDESIGNER_SDK_EXPORT QDesignerObjectInspectorInterface: public QWidget
{
    Q_OBJECT
public:
    explicit QDesignerObjectInspectorInterface(QWidget *parent, Qt::WindowFlags flags = Qt::WindowFlags());
    virtual ~QDesignerObjectInspectorInterface();

    virtual QDesignerFormEditorInterface *core() const;

public Q_SLOTS:
    virtual void setFormWindow(QDesignerFormWindowInterface *formWindow) = 0;
};

QT_END_NAMESPACE

#endif // ABSTRACTOBJECTINSPECTOR_H
