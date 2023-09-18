/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the Qt Quick Templates 2 module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL3$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see http://www.qt.io/terms-conditions. For further
** information use the contact form at http://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPLv3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or later as published by the Free
** Software Foundation and appearing in the file LICENSE.GPL included in
** the packaging of this file. Please review the following information to
** ensure the GNU General Public License version 2.0 requirements will be
** met: http://www.gnu.org/licenses/gpl-2.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QQUICKCONTROL_P_P_H
#define QQUICKCONTROL_P_P_H

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

#include <QtQuickTemplates2/private/qquickcontrol_p.h>
#include <QtQuickTemplates2/private/qquickdeferredpointer_p_p.h>
#include <QtQuickTemplates2/private/qquicktheme_p.h>

#include <QtQuick/private/qquickitem_p.h>
#include <QtQuick/private/qquickitemchangelistener_p.h>
#include <QtQml/private/qlazilyallocated_p.h>

#if QT_CONFIG(accessibility)
#include <QtGui/qaccessible.h>
#endif

#include <QtCore/qloggingcategory.h>

QT_BEGIN_NAMESPACE

Q_DECLARE_LOGGING_CATEGORY(lcItemManagement)

class QQuickAccessibleAttached;

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickControlPrivate : public QQuickItemPrivate, public QQuickItemChangeListener
#if QT_CONFIG(accessibility)
    , public QAccessible::ActivationObserver
#endif
{
    Q_DECLARE_PUBLIC(QQuickControl)

public:
    QQuickControlPrivate();
    ~QQuickControlPrivate();

    static QQuickControlPrivate *get(QQuickControl *control)
    {
        return control->d_func();
    }

    void init();

#if QT_CONFIG(quicktemplates2_multitouch)
    virtual bool acceptTouch(const QTouchEvent::TouchPoint &point);
#endif
    virtual void handlePress(const QPointF &point);
    virtual void handleMove(const QPointF &point);
    virtual void handleRelease(const QPointF &point);
    virtual void handleUngrab();

    void mirrorChange() override;

    inline QMarginsF getPadding() const { return QMarginsF(getLeftPadding(), getTopPadding(), getRightPadding(), getBottomPadding()); }
    inline qreal getTopPadding() const { return extra.isAllocated() && extra->hasTopPadding ? extra->topPadding : getVerticalPadding(); }
    inline qreal getLeftPadding() const { return extra.isAllocated() && extra->hasLeftPadding ? extra->leftPadding : getHorizontalPadding(); }
    inline qreal getRightPadding() const { return extra.isAllocated() && extra->hasRightPadding ? extra->rightPadding : getHorizontalPadding(); }
    inline qreal getBottomPadding() const { return extra.isAllocated() && extra->hasBottomPadding ? extra->bottomPadding : getVerticalPadding(); }
    inline qreal getHorizontalPadding() const { return hasHorizontalPadding ? horizontalPadding : padding; }
    inline qreal getVerticalPadding() const { return hasVerticalPadding ? verticalPadding : padding; }

    void setTopPadding(qreal value, bool reset = false);
    void setLeftPadding(qreal value, bool reset = false);
    void setRightPadding(qreal value, bool reset = false);
    void setBottomPadding(qreal value, bool reset = false);
    void setHorizontalPadding(qreal value, bool reset = false);
    void setVerticalPadding(qreal value, bool reset = false);

    inline QMarginsF getInset() const { return QMarginsF(getLeftInset(), getTopInset(), getRightInset(), getBottomInset()); }
    inline qreal getTopInset() const { return extra.isAllocated() ? extra->topInset : 0; }
    inline qreal getLeftInset() const { return extra.isAllocated() ? extra->leftInset : 0; }
    inline qreal getRightInset() const { return extra.isAllocated() ? extra->rightInset : 0; }
    inline qreal getBottomInset() const { return extra.isAllocated() ? extra->bottomInset : 0; }

    void setTopInset(qreal value, bool reset = false);
    void setLeftInset(qreal value, bool reset = false);
    void setRightInset(qreal value, bool reset = false);
    void setBottomInset(qreal value, bool reset = false);

    virtual void resizeBackground();
    virtual void resizeContent();

    virtual QQuickItem *getContentItem();
    void setContentItem_helper(QQuickItem *item, bool notify = true);

#if QT_CONFIG(accessibility)
    void accessibilityActiveChanged(bool active) override;
    QAccessible::Role accessibleRole() const override;
    static QQuickAccessibleAttached *accessibleAttached(const QObject *object);
#endif

    virtual void resolveFont();
    void inheritFont(const QFont &font);
    void updateFont(const QFont &font);
    static void updateFontRecur(QQuickItem *item, const QFont &font);
    inline void setFont_helper(const QFont &font) {
        if (resolvedFont.resolve() == font.resolve() && resolvedFont == font)
            return;
        updateFont(font);
    }
    static QFont parentFont(const QQuickItem *item);

    virtual void resolvePalette();
    void inheritPalette(const QPalette &palette);
    void updatePalette(const QPalette &palette);
    static void updatePaletteRecur(QQuickItem *item, const QPalette &palette);
    inline void setPalette_helper(const QPalette &palette) {
        if (resolvedPalette.resolve() == palette.resolve() && resolvedPalette == palette)
            return;
        updatePalette(palette);
    }
    static QPalette parentPalette(const QQuickItem *item);

    void updateLocale(const QLocale &l, bool e);
    static void updateLocaleRecur(QQuickItem *item, const QLocale &l);
    static QLocale calcLocale(const QQuickItem *item);

#if QT_CONFIG(quicktemplates2_hover)
    void updateHoverEnabled(bool enabled, bool xplicit);
    static void updateHoverEnabledRecur(QQuickItem *item, bool enabled);
    static bool calcHoverEnabled(const QQuickItem *item);
#endif

    virtual void cancelContentItem();
    virtual void executeContentItem(bool complete = false);

    virtual void cancelBackground();
    virtual void executeBackground(bool complete = false);

    static void hideOldItem(QQuickItem *item);

    void updateBaselineOffset();

    static const ChangeTypes ImplicitSizeChanges;

    void addImplicitSizeListener(QQuickItem *item, ChangeTypes changes = ImplicitSizeChanges);
    void removeImplicitSizeListener(QQuickItem *item, ChangeTypes changes = ImplicitSizeChanges);

    static void addImplicitSizeListener(QQuickItem *item, QQuickItemChangeListener *listener, ChangeTypes changes = ImplicitSizeChanges);
    static void removeImplicitSizeListener(QQuickItem *item, QQuickItemChangeListener *listener, ChangeTypes changes = ImplicitSizeChanges);

    void itemImplicitWidthChanged(QQuickItem *item) override;
    void itemImplicitHeightChanged(QQuickItem *item) override;
    void itemGeometryChanged(QQuickItem *item, QQuickGeometryChange change, const QRectF &diff) override;
    void itemDestroyed(QQuickItem *item) override;

    virtual qreal getContentWidth() const;
    virtual qreal getContentHeight() const;

    void updateImplicitContentWidth();
    void updateImplicitContentHeight();
    void updateImplicitContentSize();

    struct ExtraData {
        bool hasTopPadding = false;
        bool hasLeftPadding = false;
        bool hasRightPadding = false;
        bool hasBottomPadding = false;
        bool hasBaselineOffset = false;
        bool hasTopInset = false;
        bool hasLeftInset = false;
        bool hasRightInset = false;
        bool hasBottomInset = false;
        bool hasBackgroundWidth = false;
        bool hasBackgroundHeight = false;
        qreal topPadding = 0;
        qreal leftPadding = 0;
        qreal rightPadding = 0;
        qreal bottomPadding = 0;
        qreal topInset = 0;
        qreal leftInset = 0;
        qreal rightInset = 0;
        qreal bottomInset = 0;
        QFont requestedFont;
        QPalette requestedPalette;
    };
    QLazilyAllocated<ExtraData> extra;

    bool hasHorizontalPadding = false;
    bool hasVerticalPadding = false;
    bool hasLocale = false;
    bool wheelEnabled = false;
#if QT_CONFIG(quicktemplates2_hover)
    bool hovered = false;
    bool explicitHoverEnabled = false;
#endif
    bool resizingBackground = false;
    bool pressWasTouch = false;
    int touchId = -1;
    QPointF previousPressPos;
    qreal padding = 0;
    qreal horizontalPadding = 0;
    qreal verticalPadding = 0;
    qreal implicitContentWidth = 0;
    qreal implicitContentHeight = 0;
    qreal spacing = 0;
    QLocale locale;
    QFont resolvedFont;
    QPalette resolvedPalette;
    Qt::FocusPolicy focusPolicy = Qt::NoFocus;
    Qt::FocusReason focusReason = Qt::OtherFocusReason;
    QQuickDeferredPointer<QQuickItem> background;
    QQuickDeferredPointer<QQuickItem> contentItem;
};

QT_END_NAMESPACE

#endif // QQUICKCONTROL_P_P_H
