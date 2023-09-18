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

#ifndef QQUICKCONTROL_P_H
#define QQUICKCONTROL_P_H

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

#include <QtCore/qlocale.h>
#include <QtGui/qpalette.h>
#include <QtQuick/qquickitem.h>
#include <QtQuickTemplates2/private/qtquicktemplates2global_p.h>

QT_BEGIN_NAMESPACE

class QQuickControlPrivate;

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickControl : public QQuickItem
{
    Q_OBJECT
    Q_PROPERTY(QFont font READ font WRITE setFont RESET resetFont NOTIFY fontChanged FINAL)
    Q_PROPERTY(qreal availableWidth READ availableWidth NOTIFY availableWidthChanged FINAL)
    Q_PROPERTY(qreal availableHeight READ availableHeight NOTIFY availableHeightChanged FINAL)
    Q_PROPERTY(qreal padding READ padding WRITE setPadding RESET resetPadding NOTIFY paddingChanged FINAL)
    Q_PROPERTY(qreal topPadding READ topPadding WRITE setTopPadding RESET resetTopPadding NOTIFY topPaddingChanged FINAL)
    Q_PROPERTY(qreal leftPadding READ leftPadding WRITE setLeftPadding RESET resetLeftPadding NOTIFY leftPaddingChanged FINAL)
    Q_PROPERTY(qreal rightPadding READ rightPadding WRITE setRightPadding RESET resetRightPadding NOTIFY rightPaddingChanged FINAL)
    Q_PROPERTY(qreal bottomPadding READ bottomPadding WRITE setBottomPadding RESET resetBottomPadding NOTIFY bottomPaddingChanged FINAL)
    Q_PROPERTY(qreal spacing READ spacing WRITE setSpacing RESET resetSpacing NOTIFY spacingChanged FINAL)
    Q_PROPERTY(QLocale locale READ locale WRITE setLocale RESET resetLocale NOTIFY localeChanged FINAL)
    Q_PROPERTY(bool mirrored READ isMirrored NOTIFY mirroredChanged FINAL)
    Q_PROPERTY(Qt::FocusPolicy focusPolicy READ focusPolicy WRITE setFocusPolicy NOTIFY focusPolicyChanged FINAL)
    Q_PROPERTY(Qt::FocusReason focusReason READ focusReason WRITE setFocusReason NOTIFY focusReasonChanged FINAL)
    Q_PROPERTY(bool visualFocus READ hasVisualFocus NOTIFY visualFocusChanged FINAL)
    Q_PROPERTY(bool hovered READ isHovered NOTIFY hoveredChanged FINAL)
    Q_PROPERTY(bool hoverEnabled READ isHoverEnabled WRITE setHoverEnabled RESET resetHoverEnabled NOTIFY hoverEnabledChanged FINAL)
    Q_PROPERTY(bool wheelEnabled READ isWheelEnabled WRITE setWheelEnabled NOTIFY wheelEnabledChanged FINAL)
    Q_PROPERTY(QQuickItem *background READ background WRITE setBackground NOTIFY backgroundChanged FINAL)
    Q_PROPERTY(QQuickItem *contentItem READ contentItem WRITE setContentItem NOTIFY contentItemChanged FINAL)
    Q_PROPERTY(qreal baselineOffset READ baselineOffset WRITE setBaselineOffset RESET resetBaselineOffset NOTIFY baselineOffsetChanged FINAL)
    // 2.3 (Qt 5.10)
    Q_PROPERTY(QPalette palette READ palette WRITE setPalette RESET resetPalette NOTIFY paletteChanged FINAL REVISION 3)
    // 2.5 (Qt 5.12)
    Q_PROPERTY(qreal horizontalPadding READ horizontalPadding WRITE setHorizontalPadding RESET resetHorizontalPadding NOTIFY horizontalPaddingChanged FINAL REVISION 5)
    Q_PROPERTY(qreal verticalPadding READ verticalPadding WRITE setVerticalPadding RESET resetVerticalPadding NOTIFY verticalPaddingChanged FINAL REVISION 5)
    Q_PROPERTY(qreal implicitContentWidth READ implicitContentWidth NOTIFY implicitContentWidthChanged FINAL REVISION 5)
    Q_PROPERTY(qreal implicitContentHeight READ implicitContentHeight NOTIFY implicitContentHeightChanged FINAL REVISION 5)
    Q_PROPERTY(qreal implicitBackgroundWidth READ implicitBackgroundWidth NOTIFY implicitBackgroundWidthChanged FINAL REVISION 5)
    Q_PROPERTY(qreal implicitBackgroundHeight READ implicitBackgroundHeight NOTIFY implicitBackgroundHeightChanged FINAL REVISION 5)
    Q_PROPERTY(qreal topInset READ topInset WRITE setTopInset RESET resetTopInset NOTIFY topInsetChanged FINAL REVISION 5)
    Q_PROPERTY(qreal leftInset READ leftInset WRITE setLeftInset RESET resetLeftInset NOTIFY leftInsetChanged FINAL REVISION 5)
    Q_PROPERTY(qreal rightInset READ rightInset WRITE setRightInset RESET resetRightInset NOTIFY rightInsetChanged FINAL REVISION 5)
    Q_PROPERTY(qreal bottomInset READ bottomInset WRITE setBottomInset RESET resetBottomInset NOTIFY bottomInsetChanged FINAL REVISION 5)
    Q_CLASSINFO("DeferredPropertyNames", "background,contentItem")

public:
    explicit QQuickControl(QQuickItem *parent = nullptr);
    ~QQuickControl();

    QFont font() const;
    void setFont(const QFont &font);
    void resetFont();

    qreal availableWidth() const;
    qreal availableHeight() const;

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

    qreal spacing() const;
    void setSpacing(qreal spacing);
    void resetSpacing();

    QLocale locale() const;
    void setLocale(const QLocale &locale);
    void resetLocale();

    bool isMirrored() const;

    Qt::FocusPolicy focusPolicy() const;
    void setFocusPolicy(Qt::FocusPolicy policy);

    Qt::FocusReason focusReason() const;
    void setFocusReason(Qt::FocusReason reason);

    bool hasVisualFocus() const;

    bool isHovered() const;
    void setHovered(bool hovered);

    bool isHoverEnabled() const;
    void setHoverEnabled(bool enabled);
    void resetHoverEnabled();

    bool isWheelEnabled() const;
    void setWheelEnabled(bool enabled);

    QQuickItem *background() const;
    void setBackground(QQuickItem *background);

    QQuickItem *contentItem() const;
    void setContentItem(QQuickItem *item);

    qreal baselineOffset() const;
    void setBaselineOffset(qreal offset);
    void resetBaselineOffset();

    // 2.3 (Qt 5.10)
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

Q_SIGNALS:
    void fontChanged();
    void availableWidthChanged();
    void availableHeightChanged();
    void paddingChanged();
    void topPaddingChanged();
    void leftPaddingChanged();
    void rightPaddingChanged();
    void bottomPaddingChanged();
    void spacingChanged();
    void localeChanged();
    void mirroredChanged();
    void focusPolicyChanged();
    void focusReasonChanged();
    void visualFocusChanged();
    void hoveredChanged();
    void hoverEnabledChanged();
    void wheelEnabledChanged();
    void backgroundChanged();
    void contentItemChanged();
    void baselineOffsetChanged();
    // 2.3 (Qt 5.10)
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
    virtual QFont defaultFont() const;
    virtual QPalette defaultPalette() const;

    QQuickControl(QQuickControlPrivate &dd, QQuickItem *parent);

    void classBegin() override;
    void componentComplete() override;

    void itemChange(ItemChange change, const ItemChangeData &value) override;

    void focusInEvent(QFocusEvent *event) override;
    void focusOutEvent(QFocusEvent *event) override;
#if QT_CONFIG(quicktemplates2_hover)
    void hoverEnterEvent(QHoverEvent *event) override;
    void hoverMoveEvent(QHoverEvent *event) override;
    void hoverLeaveEvent(QHoverEvent *event) override;
#endif
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void mouseUngrabEvent() override;
#if QT_CONFIG(quicktemplates2_multitouch)
    void touchEvent(QTouchEvent *event) override;
    void touchUngrabEvent() override;
#endif
#if QT_CONFIG(wheelevent)
    void wheelEvent(QWheelEvent *event) override;
#endif

    void geometryChanged(const QRectF &newGeometry, const QRectF &oldGeometry) override;

    virtual void fontChange(const QFont &newFont, const QFont &oldFont);
#if QT_CONFIG(quicktemplates2_hover)
    virtual void hoverChange();
#endif
    virtual void mirrorChange();
    virtual void spacingChange(qreal newSpacing, qreal oldSpacing);
    virtual void paddingChange(const QMarginsF &newPadding, const QMarginsF &oldPadding);
    virtual void contentItemChange(QQuickItem *newItem, QQuickItem *oldItem);
    virtual void localeChange(const QLocale &newLocale, const QLocale &oldLocale);
    virtual void paletteChange(const QPalette &newPalette, const QPalette &oldPalette);
    virtual void insetChange(const QMarginsF &newInset, const QMarginsF &oldInset);
    virtual void enabledChange();

#if QT_CONFIG(accessibility)
    virtual QAccessible::Role accessibleRole() const;
    virtual void accessibilityActiveChanged(bool active);
#endif

    // helper functions which avoid to check QT_CONFIG(accessibility)
    QString accessibleName() const;
    void maybeSetAccessibleName(const QString &name);

    QVariant accessibleProperty(const char *propertyName);
    bool setAccessibleProperty(const char *propertyName, const QVariant &value);

private:
    Q_DISABLE_COPY(QQuickControl)
    Q_DECLARE_PRIVATE(QQuickControl)
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QQuickControl)

#endif // QQUICKCONTROL_P_H
