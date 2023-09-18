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

#ifndef QQUICKPOPUP_P_H
#define QQUICKPOPUP_P_H

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

#include <QtCore/qobject.h>
#include <QtCore/qmargins.h>
#include <QtGui/qevent.h>
#include <QtCore/qlocale.h>
#include <QtGui/qfont.h>
#include <QtGui/qpalette.h>
#include <QtQuickTemplates2/private/qtquicktemplates2global_p.h>
#include <QtQml/qqml.h>
#include <QtQml/qqmllist.h>
#include <QtQml/qqmlparserstatus.h>
#include <QtQuick/qquickitem.h>

#if QT_CONFIG(accessibility)
#include <QtGui/qaccessible.h>
#endif

QT_BEGIN_NAMESPACE

class QQuickWindow;
class QQuickPopupAnchors;
class QQuickPopupPrivate;
class QQuickTransition;

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickPopup : public QObject, public QQmlParserStatus
{
    Q_OBJECT
    Q_INTERFACES(QQmlParserStatus)
    Q_PROPERTY(qreal x READ x WRITE setX NOTIFY xChanged FINAL)
    Q_PROPERTY(qreal y READ y WRITE setY NOTIFY yChanged FINAL)
    Q_PROPERTY(qreal z READ z WRITE setZ NOTIFY zChanged FINAL)
    Q_PROPERTY(qreal width READ width WRITE setWidth RESET resetWidth NOTIFY widthChanged FINAL)
    Q_PROPERTY(qreal height READ height WRITE setHeight RESET resetHeight NOTIFY heightChanged FINAL)
    Q_PROPERTY(qreal implicitWidth READ implicitWidth WRITE setImplicitWidth NOTIFY implicitWidthChanged FINAL)
    Q_PROPERTY(qreal implicitHeight READ implicitHeight WRITE setImplicitHeight NOTIFY implicitHeightChanged FINAL)
    Q_PROPERTY(qreal contentWidth READ contentWidth WRITE setContentWidth NOTIFY contentWidthChanged FINAL)
    Q_PROPERTY(qreal contentHeight READ contentHeight WRITE setContentHeight NOTIFY contentHeightChanged FINAL)
    Q_PROPERTY(qreal availableWidth READ availableWidth NOTIFY availableWidthChanged FINAL)
    Q_PROPERTY(qreal availableHeight READ availableHeight NOTIFY availableHeightChanged FINAL)
    Q_PROPERTY(qreal margins READ margins WRITE setMargins RESET resetMargins NOTIFY marginsChanged FINAL)
    Q_PROPERTY(qreal topMargin READ topMargin WRITE setTopMargin RESET resetTopMargin NOTIFY topMarginChanged FINAL)
    Q_PROPERTY(qreal leftMargin READ leftMargin WRITE setLeftMargin RESET resetLeftMargin NOTIFY leftMarginChanged FINAL)
    Q_PROPERTY(qreal rightMargin READ rightMargin WRITE setRightMargin RESET resetRightMargin NOTIFY rightMarginChanged FINAL)
    Q_PROPERTY(qreal bottomMargin READ bottomMargin WRITE setBottomMargin RESET resetBottomMargin NOTIFY bottomMarginChanged FINAL)
    Q_PROPERTY(qreal padding READ padding WRITE setPadding RESET resetPadding NOTIFY paddingChanged FINAL)
    Q_PROPERTY(qreal topPadding READ topPadding WRITE setTopPadding RESET resetTopPadding NOTIFY topPaddingChanged FINAL)
    Q_PROPERTY(qreal leftPadding READ leftPadding WRITE setLeftPadding RESET resetLeftPadding NOTIFY leftPaddingChanged FINAL)
    Q_PROPERTY(qreal rightPadding READ rightPadding WRITE setRightPadding RESET resetRightPadding NOTIFY rightPaddingChanged FINAL)
    Q_PROPERTY(qreal bottomPadding READ bottomPadding WRITE setBottomPadding RESET resetBottomPadding NOTIFY bottomPaddingChanged FINAL)
    Q_PROPERTY(QLocale locale READ locale WRITE setLocale RESET resetLocale NOTIFY localeChanged FINAL)
    Q_PROPERTY(QFont font READ font WRITE setFont RESET resetFont NOTIFY fontChanged FINAL)
    Q_PROPERTY(QQuickItem *parent READ parentItem WRITE setParentItem RESET resetParentItem NOTIFY parentChanged FINAL)
    Q_PROPERTY(QQuickItem *background READ background WRITE setBackground NOTIFY backgroundChanged FINAL)
    Q_PROPERTY(QQuickItem *contentItem READ contentItem WRITE setContentItem NOTIFY contentItemChanged FINAL)
    Q_PRIVATE_PROPERTY(QQuickPopup::d_func(), QQmlListProperty<QObject> contentData READ contentData FINAL)
    Q_PRIVATE_PROPERTY(QQuickPopup::d_func(), QQmlListProperty<QQuickItem> contentChildren READ contentChildren NOTIFY contentChildrenChanged FINAL)
    Q_PROPERTY(bool clip READ clip WRITE setClip NOTIFY clipChanged FINAL)
    Q_PROPERTY(bool focus READ hasFocus WRITE setFocus NOTIFY focusChanged FINAL)
    Q_PROPERTY(bool activeFocus READ hasActiveFocus NOTIFY activeFocusChanged FINAL)
    Q_PROPERTY(bool modal READ isModal WRITE setModal NOTIFY modalChanged FINAL)
    Q_PROPERTY(bool dim READ dim WRITE setDim RESET resetDim NOTIFY dimChanged FINAL)
    Q_PROPERTY(bool visible READ isVisible WRITE setVisible NOTIFY visibleChanged FINAL)
    Q_PROPERTY(qreal opacity READ opacity WRITE setOpacity NOTIFY opacityChanged FINAL)
    Q_PROPERTY(qreal scale READ scale WRITE setScale NOTIFY scaleChanged FINAL)
    Q_PROPERTY(ClosePolicy closePolicy READ closePolicy WRITE setClosePolicy RESET resetClosePolicy NOTIFY closePolicyChanged FINAL)
    Q_PROPERTY(TransformOrigin transformOrigin READ transformOrigin WRITE setTransformOrigin FINAL)
    Q_PROPERTY(QQuickTransition *enter READ enter WRITE setEnter NOTIFY enterChanged FINAL)
    Q_PROPERTY(QQuickTransition *exit READ exit WRITE setExit NOTIFY exitChanged FINAL)
    // 2.1 (Qt 5.8)
    Q_PROPERTY(qreal spacing READ spacing WRITE setSpacing RESET resetSpacing NOTIFY spacingChanged FINAL REVISION 1)
    // 2.3 (Qt 5.10)
    Q_PROPERTY(bool opened READ isOpened NOTIFY openedChanged FINAL REVISION 3)
    Q_PROPERTY(bool mirrored READ isMirrored NOTIFY mirroredChanged FINAL REVISION 3)
    Q_PROPERTY(bool enabled READ isEnabled WRITE setEnabled NOTIFY enabledChanged FINAL REVISION 3)
    Q_PROPERTY(QPalette palette READ palette WRITE setPalette RESET resetPalette NOTIFY paletteChanged FINAL REVISION 3)
    // 2.5 (Qt 5.12)
    Q_PROPERTY(qreal horizontalPadding READ horizontalPadding WRITE setHorizontalPadding RESET resetHorizontalPadding NOTIFY horizontalPaddingChanged FINAL)
    Q_PROPERTY(qreal verticalPadding READ verticalPadding WRITE setVerticalPadding RESET resetVerticalPadding NOTIFY verticalPaddingChanged FINAL)
    Q_PRIVATE_PROPERTY(QQuickPopup::d_func(), QQuickPopupAnchors *anchors READ getAnchors DESIGNABLE false CONSTANT FINAL REVISION 5)
    Q_PROPERTY(qreal implicitContentWidth READ implicitContentWidth NOTIFY implicitContentWidthChanged FINAL REVISION 5)
    Q_PROPERTY(qreal implicitContentHeight READ implicitContentHeight NOTIFY implicitContentHeightChanged FINAL REVISION 5)
    Q_PROPERTY(qreal implicitBackgroundWidth READ implicitBackgroundWidth NOTIFY implicitBackgroundWidthChanged FINAL REVISION 5)
    Q_PROPERTY(qreal implicitBackgroundHeight READ implicitBackgroundHeight NOTIFY implicitBackgroundHeightChanged FINAL REVISION 5)
    Q_PROPERTY(qreal topInset READ topInset WRITE setTopInset RESET resetTopInset NOTIFY topInsetChanged FINAL REVISION 5)
    Q_PROPERTY(qreal leftInset READ leftInset WRITE setLeftInset RESET resetLeftInset NOTIFY leftInsetChanged FINAL REVISION 5)
    Q_PROPERTY(qreal rightInset READ rightInset WRITE setRightInset RESET resetRightInset NOTIFY rightInsetChanged FINAL REVISION 5)
    Q_PROPERTY(qreal bottomInset READ bottomInset WRITE setBottomInset RESET resetBottomInset NOTIFY bottomInsetChanged FINAL REVISION 5)
    Q_CLASSINFO("DeferredPropertyNames", "background,contentItem")
    Q_CLASSINFO("DefaultProperty", "contentData")

public:
    explicit QQuickPopup(QObject *parent = nullptr);
    ~QQuickPopup();

    qreal x() const;
    void setX(qreal x);

    qreal y() const;
    void setY(qreal y);

    QPointF position() const;
    void setPosition(const QPointF &pos);

    qreal z() const;
    void setZ(qreal z);

    qreal width() const;
    void setWidth(qreal width);
    void resetWidth();

    qreal height() const;
    void setHeight(qreal height);
    void resetHeight();

    qreal implicitWidth() const;
    void setImplicitWidth(qreal width);

    qreal implicitHeight() const;
    void setImplicitHeight(qreal height);

    qreal contentWidth() const;
    void setContentWidth(qreal width);

    qreal contentHeight() const;
    void setContentHeight(qreal height);

    qreal availableWidth() const;
    qreal availableHeight() const;

    qreal margins() const;
    void setMargins(qreal margins);
    void resetMargins();

    qreal topMargin() const;
    void setTopMargin(qreal margin);
    void resetTopMargin();

    qreal leftMargin() const;
    void setLeftMargin(qreal margin);
    void resetLeftMargin();

    qreal rightMargin() const;
    void setRightMargin(qreal margin);
    void resetRightMargin();

    qreal bottomMargin() const;
    void setBottomMargin(qreal margin);
    void resetBottomMargin();

    qreal padding() const;
    void setPadding(qreal padding);
    void resetPadding();

    qreal topPadding() const;
    void setTopPadding(qreal padding);
    void resetTopPadding();

    qreal leftPadding() const;
    void setLeftPadding(qreal padding);
    void resetLeftPadding();

    qreal rightPadding() const;
    void setRightPadding(qreal padding);
    void resetRightPadding();

    qreal bottomPadding() const;
    void setBottomPadding(qreal padding);
    void resetBottomPadding();

    QLocale locale() const;
    void setLocale(const QLocale &locale);
    void resetLocale();

    QFont font() const;
    void setFont(const QFont &font);
    void resetFont();

    QQuickWindow *window() const;
    QQuickItem *popupItem() const;

    QQuickItem *parentItem() const;
    void setParentItem(QQuickItem *parent);
    void resetParentItem();

    QQuickItem *background() const;
    void setBackground(QQuickItem *background);

    QQuickItem *contentItem() const;
    void setContentItem(QQuickItem *item);

    bool clip() const;
    void setClip(bool clip);

    bool hasFocus() const;
    void setFocus(bool focus);

    bool hasActiveFocus() const;

    bool isModal() const;
    void setModal(bool modal);

    bool dim() const;
    void setDim(bool dim);
    void resetDim();

    bool isVisible() const;
    virtual void setVisible(bool visible);

    qreal opacity() const;
    void setOpacity(qreal opacity);

    qreal scale() const;
    void setScale(qreal scale);

    enum ClosePolicyFlag {
        NoAutoClose = 0x00,
        CloseOnPressOutside = 0x01,
        CloseOnPressOutsideParent = 0x02,
        CloseOnReleaseOutside = 0x04,
        CloseOnReleaseOutsideParent = 0x08,
        CloseOnEscape = 0x10
    };
    Q_DECLARE_FLAGS(ClosePolicy, ClosePolicyFlag)
    Q_FLAG(ClosePolicy)

    ClosePolicy closePolicy() const;
    void setClosePolicy(ClosePolicy policy);
    void resetClosePolicy();

    // keep in sync with Item.TransformOrigin
    enum TransformOrigin {
        TopLeft, Top, TopRight,
        Left, Center, Right,
        BottomLeft, Bottom, BottomRight
    };
    Q_ENUM(TransformOrigin)

    TransformOrigin transformOrigin() const;
    void setTransformOrigin(TransformOrigin);

    QQuickTransition *enter() const;
    void setEnter(QQuickTransition *transition);

    QQuickTransition *exit() const;
    void setExit(QQuickTransition *transition);

    bool filtersChildMouseEvents() const;
    void setFiltersChildMouseEvents(bool filter);

    Q_INVOKABLE void forceActiveFocus(Qt::FocusReason reason = Qt::OtherFocusReason);

    // 2.1 (Qt 5.8)
    qreal spacing() const;
    void setSpacing(qreal spacing);
    void resetSpacing();

    // 2.3 (Qt 5.10)
    bool isOpened() const;
    bool isMirrored() const;

    bool isEnabled() const;
    void setEnabled(bool enabled);

    QPalette palette() const;
    void setPalette(const QPalette &palette);
    void resetPalette();

    // 2.5 (Qt 5.12)
    qreal horizontalPadding() const;
    void setHorizontalPadding(qreal padding);
    void resetHorizontalPadding();

    qreal verticalPadding() const;
    void setVerticalPadding(qreal padding);
    void resetVerticalPadding();

    qreal implicitContentWidth() const;
    qreal implicitContentHeight() const;

    qreal implicitBackgroundWidth() const;
    qreal implicitBackgroundHeight() const;

    qreal topInset() const;
    void setTopInset(qreal inset);
    void resetTopInset();

    qreal leftInset() const;
    void setLeftInset(qreal inset);
    void resetLeftInset();

    qreal rightInset() const;
    void setRightInset(qreal inset);
    void resetRightInset();

    qreal bottomInset() const;
    void setBottomInset(qreal inset);
    void resetBottomInset();

public Q_SLOTS:
    void open();
    void close();

Q_SIGNALS:
    void opened();
    void closed();
    void aboutToShow();
    void aboutToHide();
    void xChanged();
    void yChanged();
    void zChanged();
    void widthChanged();
    void heightChanged();
    void implicitWidthChanged();
    void implicitHeightChanged();
    void contentWidthChanged();
    void contentHeightChanged();
    void availableWidthChanged();
    void availableHeightChanged();
    void marginsChanged();
    void topMarginChanged();
    void leftMarginChanged();
    void rightMarginChanged();
    void bottomMarginChanged();
    void paddingChanged();
    void topPaddingChanged();
    void leftPaddingChanged();
    void rightPaddingChanged();
    void bottomPaddingChanged();
    void fontChanged();
    void localeChanged();
    void parentChanged();
    void backgroundChanged();
    void contentItemChanged();
    void contentChildrenChanged();
    void clipChanged();
    void focusChanged();
    void activeFocusChanged();
    void modalChanged();
    void dimChanged();
    void visibleChanged();
    void opacityChanged();
    void scaleChanged();
    void closePolicyChanged();
    void enterChanged();
    void exitChanged();
    void windowChanged(QQuickWindow *window);
    // 2.1 (Qt 5.8)
    Q_REVISION(1) void spacingChanged();
    // 2.3 (Qt 5.10)
    Q_REVISION(3) void openedChanged();
    Q_REVISION(3) void mirroredChanged();
    Q_REVISION(3) void enabledChanged();
    Q_REVISION(3) void paletteChanged();
    // 2.5 (Qt 5.12)
    Q_REVISION(5) void horizontalPaddingChanged();
    Q_REVISION(5) void verticalPaddingChanged();
    Q_REVISION(5) void implicitContentWidthChanged();
    Q_REVISION(5) void implicitContentHeightChanged();
    Q_REVISION(5) void implicitBackgroundWidthChanged();
    Q_REVISION(5) void implicitBackgroundHeightChanged();
    Q_REVISION(5) void topInsetChanged();
    Q_REVISION(5) void leftInsetChanged();
    Q_REVISION(5) void rightInsetChanged();
    Q_REVISION(5) void bottomInsetChanged();

protected:
    QQuickPopup(QQuickPopupPrivate &dd, QObject *parent);

    void classBegin() override;
    void componentComplete() override;
    bool isComponentComplete() const;

    virtual bool childMouseEventFilter(QQuickItem *child, QEvent *event);
    virtual void focusInEvent(QFocusEvent *event);
    virtual void focusOutEvent(QFocusEvent *event);
    virtual void keyPressEvent(QKeyEvent *event);
    virtual void keyReleaseEvent(QKeyEvent *event);
    virtual void mousePressEvent(QMouseEvent *event);
    virtual void mouseMoveEvent(QMouseEvent *event);
    virtual void mouseReleaseEvent(QMouseEvent *event);
    virtual void mouseDoubleClickEvent(QMouseEvent *event);
    virtual void mouseUngrabEvent();
    virtual bool overlayEvent(QQuickItem *item, QEvent *event);
#if QT_CONFIG(quicktemplates2_multitouch)
    virtual void touchEvent(QTouchEvent *event);
    virtual void touchUngrabEvent();
#endif
#if QT_CONFIG(wheelevent)
    virtual void wheelEvent(QWheelEvent *event);
#endif

    virtual void contentItemChange(QQuickItem *newItem, QQuickItem *oldItem);
    virtual void contentSizeChange(const QSizeF &newSize, const QSizeF &oldSize);
    virtual void fontChange(const QFont &newFont, const QFont &oldFont);
    virtual void geometryChanged(const QRectF &newGeometry, const QRectF &oldGeometry);
    virtual void localeChange(const QLocale &newLocale, const QLocale &oldLocale);
    virtual void itemChange(QQuickItem::ItemChange change, const QQuickItem::ItemChangeData &data);
    virtual void marginsChange(const QMarginsF &newMargins, const QMarginsF &oldMargins);
    virtual void paddingChange(const QMarginsF &newPadding, const QMarginsF &oldPadding);
    virtual void paletteChange(const QPalette &newPalette, const QPalette &oldPalette);
    virtual void spacingChange(qreal newSpacing, qreal oldSpacing);
    virtual void insetChange(const QMarginsF &newInset, const QMarginsF &oldInset);

    virtual QFont defaultFont() const;
    virtual QPalette defaultPalette() const;

#if QT_CONFIG(accessibility)
    virtual QAccessible::Role accessibleRole() const;
    virtual void accessibilityActiveChanged(bool active);
#endif

    QString accessibleName() const;
    void maybeSetAccessibleName(const QString &name);

    QVariant accessibleProperty(const char *propertyName);
    bool setAccessibleProperty(const char *propertyName, const QVariant &value);

private:
    Q_DISABLE_COPY(QQuickPopup)
    Q_DECLARE_PRIVATE(QQuickPopup)
    friend class QQuickPopupItem;
    friend class QQuickOverlay;
    friend class QQuickOverlayPrivate;
};

Q_DECLARE_OPERATORS_FOR_FLAGS(QQuickPopup::ClosePolicy)

QT_END_NAMESPACE

QML_DECLARE_TYPE(QQuickPopup)

#endif // QQUICKPOPUP_P_H
