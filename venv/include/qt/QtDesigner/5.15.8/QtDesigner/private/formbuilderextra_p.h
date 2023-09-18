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

#ifndef ABSTRACTFORMBUILDERPRIVATE_H
#define ABSTRACTFORMBUILDERPRIVATE_H

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

#include "uilib_global.h"

#include <QtCore/qhash.h>
#include <QtCore/qpointer.h>
#include <QtCore/qstringlist.h>
#include <QtCore/qmap.h>
#include <QtCore/qdir.h>
#include <QtGui/qpalette.h>

QT_BEGIN_NAMESPACE

class QDesignerCustomWidgetInterface;
class QObject;
class QVariant;
class QWidget;
class QObject;
class QLabel;
class QButtonGroup;
class QBoxLayout;
class QGridLayout;
class QAction;
class QActionGroup;

#ifdef QFORMINTERNAL_NAMESPACE
namespace QFormInternal
{
#endif

class DomBrush;
class DomButtonGroups;
class DomButtonGroup;
class DomColorGroup;
class DomCustomWidget;
class DomPalette;
class DomProperty;
class DomUI;

class QAbstractFormBuilder;
class QResourceBuilder;
class QTextBuilder;

class QDESIGNER_UILIB_EXPORT QFormBuilderExtra
{
public:
    QFormBuilderExtra();
    ~QFormBuilderExtra();

    struct CustomWidgetData {
        CustomWidgetData();
        explicit CustomWidgetData(const DomCustomWidget *dc);

        QString addPageMethod;
        QString script;
        QString baseClass;
        bool isContainer = false;
    };

    void clear();

    DomUI *readUi(QIODevice *dev);
    static QString msgInvalidUiFile();

    bool applyPropertyInternally(QObject *o, const QString &propertyName, const QVariant &value);

    enum BuddyMode { BuddyApplyAll, BuddyApplyVisibleOnly };

    void applyInternalProperties() const;
    static bool applyBuddy(const QString &buddyName, BuddyMode applyMode, QLabel *label);

    const QPointer<QWidget> &parentWidget() const;
    bool parentWidgetIsSet() const;
    void setParentWidget(const QPointer<QWidget> &w);

    void setProcessingLayoutWidget(bool processing);
    bool processingLayoutWidget() const;

    void setResourceBuilder(QResourceBuilder *builder);
    QResourceBuilder *resourceBuilder() const;

    void setTextBuilder(QTextBuilder *builder);
    QTextBuilder *textBuilder() const;

    void storeCustomWidgetData(const QString &className, const DomCustomWidget *d);
    QString customWidgetAddPageMethod(const QString &className) const;
    QString customWidgetBaseClass(const QString &className) const;
    bool isCustomWidgetContainer(const QString &className) const;

    // --- Hash used in creating button groups on demand. Store a map of name and pair of dom group and real group
    void registerButtonGroups(const DomButtonGroups *groups);

    using ButtonGroupEntry = QPair<DomButtonGroup *, QButtonGroup*>;
    using ButtonGroupHash = QHash<QString, ButtonGroupEntry>;
    const ButtonGroupHash &buttonGroups() const { return m_buttonGroups; }
    ButtonGroupHash &buttonGroups()  { return m_buttonGroups; }

    // return stretch as a comma-separated list
    static QString boxLayoutStretch(const QBoxLayout*);
    // apply stretch
    static bool setBoxLayoutStretch(const QString &, QBoxLayout*);
    static void clearBoxLayoutStretch(QBoxLayout*);

    static QString gridLayoutRowStretch(const QGridLayout *);
    static bool setGridLayoutRowStretch(const QString &, QGridLayout *);
    static void clearGridLayoutRowStretch(QGridLayout *);

    static QString gridLayoutColumnStretch(const QGridLayout *);
    static bool setGridLayoutColumnStretch(const QString &, QGridLayout *);
    static void clearGridLayoutColumnStretch(QGridLayout *);

    // return the row/column sizes as  comma-separated lists
    static QString gridLayoutRowMinimumHeight(const QGridLayout *);
    static bool setGridLayoutRowMinimumHeight(const QString &, QGridLayout *);
    static void clearGridLayoutRowMinimumHeight(QGridLayout *);

    static QString gridLayoutColumnMinimumWidth(const QGridLayout *);
    static bool setGridLayoutColumnMinimumWidth(const QString &, QGridLayout *);
    static void clearGridLayoutColumnMinimumWidth(QGridLayout *);

    static void setPixmapProperty(DomProperty *p, const QPair<QString, QString> &ip);
    static QPalette loadPalette(const DomPalette *dom);
    static void setupColorGroup(QPalette *palette, QPalette::ColorGroup colorGroup,
                                const DomColorGroup *group);
    static DomColorGroup *saveColorGroup(const QPalette &palette,
                                         QPalette::ColorGroup colorGroup);
    static DomPalette *savePalette(const QPalette &palette);
    static QBrush setupBrush(const DomBrush *brush);
    static DomBrush *saveBrush(const QBrush &br);

    QStringList m_pluginPaths;
    QMap<QString, QDesignerCustomWidgetInterface*> m_customWidgets;

    QHash<QObject*, bool> m_laidout;
    QHash<QString, QAction*> m_actions;
    QHash<QString, QActionGroup*> m_actionGroups;
    int m_defaultMargin;
    int m_defaultSpacing;
    QDir m_workingDirectory;
    QString m_errorString;
    QString m_language;

private:
    void clearResourceBuilder();
    void clearTextBuilder();

    using BuddyHash = QHash<QLabel*, QString>;
    BuddyHash m_buddies;

    QHash<QString, CustomWidgetData> m_customWidgetDataHash;

    ButtonGroupHash m_buttonGroups;

    bool m_layoutWidget = false;
    QResourceBuilder *m_resourceBuilder = nullptr;
    QTextBuilder *m_textBuilder = nullptr;

    QPointer<QWidget> m_parentWidget;
    bool m_parentWidgetIsSet = false;
};

void uiLibWarning(const QString &message);

// Struct with static accessor that provides most strings used in the form builder.
struct QDESIGNER_UILIB_EXPORT QFormBuilderStrings {
    QFormBuilderStrings();

    static const QFormBuilderStrings &instance();

    const QString buddyProperty;
    const QString cursorProperty;
    const QString objectNameProperty;
    const QString trueValue;
    const QString falseValue;
    const QString horizontalPostFix;
    const QString separator;
    const QString defaultTitle;
    const QString titleAttribute;
    const QString labelAttribute;
    const QString toolTipAttribute;
    const QString whatsThisAttribute;
    const QString flagsAttribute;
    const QString iconAttribute;
    const QString pixmapAttribute;
    const QString textAttribute;
    const QString currentIndexProperty;
    const QString toolBarAreaAttribute;
    const QString toolBarBreakAttribute;
    const QString dockWidgetAreaAttribute;
    const QString marginProperty;
    const QString spacingProperty;
    const QString leftMarginProperty;
    const QString topMarginProperty;
    const QString rightMarginProperty;
    const QString bottomMarginProperty;
    const QString horizontalSpacingProperty;
    const QString verticalSpacingProperty;
    const QString sizeHintProperty;
    const QString sizeTypeProperty;
    const QString orientationProperty;
    const QString styleSheetProperty;
    const QString qtHorizontal;
    const QString qtVertical;
    const QString currentRowProperty;
    const QString tabSpacingProperty;
    const QString qWidgetClass;
    const QString lineClass;
    const QString geometryProperty;
    const QString scriptWidgetVariable;
    const QString scriptChildWidgetsVariable;

    using RoleNName = QPair<Qt::ItemDataRole, QString>;
    QList<RoleNName> itemRoles;
    QHash<QString, Qt::ItemDataRole> treeItemRoleHash;

    // first.first is primary role, first.second is shadow role.
    // Shadow is used for either the translation source or the designer
    // representation of the string value.
    using TextRoleNName = QPair<QPair<Qt::ItemDataRole, Qt::ItemDataRole>, QString>;
    QList<TextRoleNName> itemTextRoles;
    QHash<QString, QPair<Qt::ItemDataRole, Qt::ItemDataRole> > treeItemTextRoleHash;
};
#ifdef QFORMINTERNAL_NAMESPACE
}
#endif

QT_END_NAMESPACE

#endif // ABSTRACTFORMBUILDERPRIVATE_H
