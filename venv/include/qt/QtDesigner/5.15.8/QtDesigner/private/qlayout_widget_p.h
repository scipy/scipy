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

#ifndef QLAYOUT_WIDGET_H
#define QLAYOUT_WIDGET_H

#include "shared_global_p.h"

#include <QtDesigner/layoutdecoration.h>

#include <QtCore/qpointer.h>
#include <QtCore/qvariant.h>
#include <QtWidgets/qwidget.h>
#include <QtWidgets/qlayout.h>

QT_BEGIN_NAMESPACE

class QDesignerFormWindowInterface;
class QDesignerFormEditorInterface;
class QGridLayout;
class QFormLayout;

namespace qdesigner_internal {
// ---- LayoutProperties: Helper struct that stores all layout-relevant properties
//      with functions to retrieve and apply to property sheets. Can be used to store the state
//      for undo commands and while rebuilding layouts.
struct QDESIGNER_SHARED_EXPORT LayoutProperties
{
    LayoutProperties();
    void clear();

    enum Margins { LeftMargin, TopMargin, RightMargin, BottomMargin, MarginCount };
    enum Spacings { Spacing, HorizSpacing, VertSpacing, SpacingsCount };

    enum PropertyMask {
        ObjectNameProperty  = 0x1,
        LeftMarginProperty = 0x2, TopMarginProperty = 0x4, RightMarginProperty = 0x8, BottomMarginProperty = 0x10,
        SpacingProperty = 0x20, HorizSpacingProperty = 0x40, VertSpacingProperty = 0x80,
        SizeConstraintProperty = 0x100,
        FieldGrowthPolicyProperty = 0x200, RowWrapPolicyProperty = 0x400, LabelAlignmentProperty = 0x0800, FormAlignmentProperty = 0x1000,
        BoxStretchProperty = 0x2000, GridRowStretchProperty = 0x4000, GridColumnStretchProperty = 0x8000,
        GridRowMinimumHeightProperty = 0x10000, GridColumnMinimumWidthProperty = 0x20000,
        AllProperties = 0xFFFF};

    // return a PropertyMask of visible properties
    static int visibleProperties(const QLayout *layout);

    // Retrieve from /apply to sheet: A property mask is returned indicating the properties found in the sheet
    int fromPropertySheet(const QDesignerFormEditorInterface *core, QLayout *l, int mask = AllProperties);
    int toPropertySheet(const QDesignerFormEditorInterface *core, QLayout *l, int mask = AllProperties, bool applyChanged = true) const;

    int m_margins[MarginCount];
    bool m_marginsChanged[MarginCount];

    int m_spacings[SpacingsCount];
    bool m_spacingsChanged[SpacingsCount];

    QVariant m_objectName; // receives a PropertySheetStringValue
    bool m_objectNameChanged;
    QVariant m_sizeConstraint;
    bool m_sizeConstraintChanged;

    bool m_fieldGrowthPolicyChanged;
    QVariant m_fieldGrowthPolicy;
    bool m_rowWrapPolicyChanged;
    QVariant m_rowWrapPolicy;
    bool m_labelAlignmentChanged;
    QVariant m_labelAlignment;
    bool m_formAlignmentChanged;
    QVariant m_formAlignment;

    bool m_boxStretchChanged;
    QVariant m_boxStretch;

    bool m_gridRowStretchChanged;
    QVariant m_gridRowStretch;

    bool m_gridColumnStretchChanged;
    QVariant m_gridColumnStretch;

    bool m_gridRowMinimumHeightChanged;
    QVariant m_gridRowMinimumHeight;

    bool m_gridColumnMinimumWidthChanged;
    QVariant  m_gridColumnMinimumWidth;
};

// -- LayoutHelper: For use with the 'insert widget'/'delete widget' command,
//    able to store and restore states.
//    This could become part of 'QDesignerLayoutDecorationExtensionV2',
//    but to keep any existing old extensions working, it is provided as
//    separate class with a factory function.
class LayoutHelper {
protected:
    LayoutHelper();

public:
    Q_DISABLE_COPY(LayoutHelper)

    virtual ~LayoutHelper();

    static LayoutHelper *createLayoutHelper(int type);

    static int indexOf(const QLayout *lt, const QWidget *widget);

    // Return area of an item (x == columns)
    QRect itemInfo(QLayout *lt, const QWidget *widget) const;

    virtual QRect itemInfo(QLayout *lt, int index) const = 0;
    virtual void insertWidget(QLayout *lt, const QRect &info, QWidget *w) = 0;
    virtual void removeWidget(QLayout *lt, QWidget *widget) = 0;
    // Since 4.5: The 'morphing' feature requires an API for replacing widgets on layouts.
    virtual void replaceWidget(QLayout *lt, QWidget *before, QWidget *after) = 0;

    // Simplify a grid, remove empty columns, rows within the rectangle
    // The DeleteWidget command restricts the area to be simplified.
    virtual bool canSimplify(const QDesignerFormEditorInterface *core, const QWidget *widgetWithManagedLayout, const QRect &restrictionArea) const = 0;
    virtual void simplify(const QDesignerFormEditorInterface *core, QWidget *widgetWithManagedLayout, const QRect &restrictionArea) = 0;

    // Push and pop a state. Can be used for implementing undo for
    // simplify/row, column insertion commands, provided that
    // the widgets remain the same.
    virtual void pushState(const QDesignerFormEditorInterface *core, const QWidget *widgetWithManagedLayout)  = 0;
    virtual void popState(const QDesignerFormEditorInterface *core, QWidget *widgetWithManagedLayout) = 0;
};

// Base class for layout decoration extensions.
class QDESIGNER_SHARED_EXPORT QLayoutSupport: public QObject, public QDesignerLayoutDecorationExtension
{
    Q_OBJECT
    Q_INTERFACES(QDesignerLayoutDecorationExtension)

protected:
    QLayoutSupport(QDesignerFormWindowInterface *formWindow, QWidget *widget, LayoutHelper *helper, QObject *parent = nullptr);

public:
    ~QLayoutSupport() override;

    inline QDesignerFormWindowInterface *formWindow() const   { return m_formWindow; }

    // DecorationExtension V2
    LayoutHelper* helper() const                              { return m_helper; }

    // DecorationExtension
    int currentIndex() const override                  { return m_currentIndex; }

    InsertMode currentInsertMode() const override      { return m_currentInsertMode; }

    QPair<int, int> currentCell() const  override      { return m_currentCell; }

    int findItemAt(const QPoint &pos) const override;
    int indexOf(QWidget *widget) const override;
    int indexOf(QLayoutItem *item) const override;

    void adjustIndicator(const QPoint &pos, int index) override;

    QWidgetList widgets(QLayout *layout) const override;

    // Pad empty cells with dummy spacers. Called by layouting commands.
    static void createEmptyCells(QGridLayout *gridLayout);
    // remove dummy spacers in the area. Returns false if there are non-empty items in the way
    static bool removeEmptyCells(QGridLayout *gridLayout, const QRect &area);
    static void createEmptyCells(QFormLayout *formLayout); // ditto.
    static bool removeEmptyCells(QFormLayout *formLayout, const QRect &area);

    // grid helpers: find item index
    static int findItemAt(QGridLayout *, int row, int column);
    // grid helpers: Quick check whether simplify should be enabled for grids. May return false positives.
    static bool canSimplifyQuickCheck(const QGridLayout *);
    static bool canSimplifyQuickCheck(const QFormLayout *fl);
    // Factory function, create layout support according to layout type of widget
    static QLayoutSupport *createLayoutSupport(QDesignerFormWindowInterface *formWindow, QWidget *widget, QObject *parent = nullptr);

protected:
    // figure out insertion position and mode from indicator on empty cell if supported
    virtual void setCurrentCellFromIndicatorOnEmptyCell(int index) = 0;
    // figure out insertion position and mode from indicator
    virtual void setCurrentCellFromIndicator(Qt::Orientation indicatorOrientation, int index, int increment) = 0;

    // Overwrite to return the extended geometry of an item, that is,
    // if it is a border item, include the widget border for the indicator to work correctly
    virtual QRect extendedGeometry(int index) const = 0;
    virtual bool supportsIndicatorOrientation(Qt::Orientation indicatorOrientation) const = 0;

    QRect itemInfo(int index) const override;
    QLayout *layout() const;
    QGridLayout *gridLayout() const;
    QWidget *widget() const              { return m_widget; }

    void setInsertMode(InsertMode im);
    void setCurrentCell(const QPair<int, int> &cell);

private:
    enum Indicator { LeftIndicator, TopIndicator, RightIndicator, BottomIndicator, NumIndicators };

    void hideIndicator(Indicator i);
    void showIndicator(Indicator i, const QRect &geometry, const QPalette &);

    QDesignerFormWindowInterface *m_formWindow;
    LayoutHelper* m_helper;

    QPointer<QWidget> m_widget;
    QPointer<QWidget> m_indicators[NumIndicators];
    int m_currentIndex;
    InsertMode m_currentInsertMode;
    QPair<int, int> m_currentCell;
};
} // namespace qdesigner_internal

// Red layout widget.
class QDESIGNER_SHARED_EXPORT QLayoutWidget: public QWidget
{
    Q_OBJECT
public:
    explicit QLayoutWidget(QDesignerFormWindowInterface *formWindow, QWidget *parent = nullptr);

    int layoutLeftMargin() const;
    void setLayoutLeftMargin(int layoutMargin);

    int layoutTopMargin() const;
    void setLayoutTopMargin(int layoutMargin);

    int layoutRightMargin() const;
    void setLayoutRightMargin(int layoutMargin);

    int layoutBottomMargin() const;
    void setLayoutBottomMargin(int layoutMargin);

    inline QDesignerFormWindowInterface *formWindow() const    { return m_formWindow; }

protected:
    bool event(QEvent *e) override;
    void paintEvent(QPaintEvent *e) override;

private:
    QDesignerFormWindowInterface *m_formWindow;
    int m_leftMargin;
    int m_topMargin;
    int m_rightMargin;
    int m_bottomMargin;
};

QT_END_NAMESPACE

#endif // QDESIGNER_WIDGET_H
