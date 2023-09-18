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

#ifndef ABSTRACTPROPERTYEDITOR_H
#define ABSTRACTPROPERTYEDITOR_H

#include <QtDesigner/sdk_global.h>

#include <QtWidgets/qwidget.h>

QT_BEGIN_NAMESPACE

class QDesignerFormEditorInterface;
class QString;
class QVariant;

class QDESIGNER_SDK_EXPORT QDesignerPropertyEditorInterface: public QWidget
{
    Q_OBJECT
public:
    explicit QDesignerPropertyEditorInterface(QWidget *parent, Qt::WindowFlags flags = Qt::WindowFlags());
    virtual ~QDesignerPropertyEditorInterface();

    virtual QDesignerFormEditorInterface *core() const;

    virtual bool isReadOnly() const = 0;
    virtual QObject *object() const = 0;

    virtual QString currentPropertyName() const = 0;

Q_SIGNALS:
    void propertyChanged(const QString &name, const QVariant &value);

public Q_SLOTS:
    virtual void setObject(QObject *object) = 0;
    virtual void setPropertyValue(const QString &name, const QVariant &value, bool changed = true) = 0;
    virtual void setReadOnly(bool readOnly) = 0;
};

QT_END_NAMESPACE

#endif // ABSTRACTPROPERTYEDITOR_H
