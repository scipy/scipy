/****************************************************************************
**
** Copyright (C) 2018 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Designer of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** BSD License Usage
** Alternatively, you may use this file under the terms of the BSD license
** as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of The Qt Company Ltd nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef FORMBUILDER_H
#define FORMBUILDER_H

#include "uilib_global.h"
#include "abstractformbuilder.h"

QT_BEGIN_NAMESPACE
#if 0
// pragma for syncqt, don't remove.

#pragma qt_class(QFormBuilder)
#endif

class QDesignerCustomWidgetInterface;

#ifdef QFORMINTERNAL_NAMESPACE
namespace QFormInternal
{
#endif

class QDESIGNER_UILIB_EXPORT QFormBuilder: public QAbstractFormBuilder
{
public:
    QFormBuilder();
    ~QFormBuilder() override;

    QStringList pluginPaths() const;

    void clearPluginPaths();
    void addPluginPath(const QString &pluginPath);
    void setPluginPath(const QStringList &pluginPaths);

    QList<QDesignerCustomWidgetInterface*> customWidgets() const;

protected:
    QWidget *create(DomUI *ui, QWidget *parentWidget) override;
    QWidget *create(DomWidget *ui_widget, QWidget *parentWidget) override;
    QLayout *create(DomLayout *ui_layout, QLayout *layout, QWidget *parentWidget) override;
    QLayoutItem *create(DomLayoutItem *ui_layoutItem, QLayout *layout, QWidget *parentWidget) override;
    QAction *create(DomAction *ui_action, QObject *parent) override;
    QActionGroup *create(DomActionGroup *ui_action_group, QObject *parent) override;

    QWidget *createWidget(const QString &widgetName, QWidget *parentWidget, const QString &name) override;
    QLayout *createLayout(const QString &layoutName, QObject *parent, const QString &name) override;

    void createConnections(DomConnections *connections, QWidget *widget) override;

    bool addItem(DomLayoutItem *ui_item, QLayoutItem *item, QLayout *layout) override;
    bool addItem(DomWidget *ui_widget, QWidget *widget, QWidget *parentWidget) override;

    virtual void updateCustomWidgets();
    void applyProperties(QObject *o, const QList<DomProperty*> &properties) override;

    static QWidget *widgetByName(QWidget *topLevel, const QString &name);

private:
};

#ifdef QFORMINTERNAL_NAMESPACE
}
#endif

QT_END_NAMESPACE

#endif // FORMBUILDER_H
