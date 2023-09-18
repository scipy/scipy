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

#ifndef QQUICKLABEL_P_P_H
#define QQUICKLABEL_P_P_H

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

#include <QtQml/private/qlazilyallocated_p.h>
#include <QtQuick/private/qquicktext_p_p.h>
#include <QtQuick/private/qquickitemchangelistener_p.h>
#include <QtQuickTemplates2/private/qquickdeferredpointer_p_p.h>

#if QT_CONFIG(accessibility)
#include <QtGui/qaccessible.h>
#endif

QT_BEGIN_NAMESPACE

class QQuickLabelPrivate : public QQuickTextPrivate, public QQuickItemChangeListener
#if QT_CONFIG(accessibility)
    , public QAccessible::ActivationObserver
#endif
{
    Q_DECLARE_PUBLIC(QQuickLabel)

public:
    QQuickLabelPrivate();
    ~QQuickLabelPrivate();

    static QQuickLabelPrivate *get(QQuickLabel *item)
    {
        return static_cast<QQuickLabelPrivate *>(QObjectPrivate::get(item));
    }

    inline QMarginsF getInset() const { return QMarginsF(getLeftInset(), getTopInset(), getRightInset(), getBottomInset()); }
    inline qreal getTopInset() const { return extra.isAllocated() ? extra->topInset : 0; }
    inline qreal getLeftInset() const { return extra.isAllocated() ? extra->leftInset : 0; }
    inline qreal getRightInset() const { return extra.isAllocated() ? extra->rightInset : 0; }
    inline qreal getBottomInset() const { return extra.isAllocated() ? extra->bottomInset : 0; }

    void setTopInset(qreal value, bool reset = false);
    void setLeftInset(qreal value, bool reset = false);
    void setRightInset(qreal value, bool reset = false);
    void setBottomInset(qreal value, bool reset = false);

    void resizeBackground();

    void resolveFont();
    void inheritFont(const QFont &font);
    void updateFont(const QFont &font);
    inline void setFont_helper(const QFont &font) {
        if (sourceFont.resolve() == font.resolve() && sourceFont == font)
            return;
        updateFont(font);
    }

    void resolvePalette();
    void inheritPalette(const QPalette &palette);
    void updatePalette(const QPalette &palette);
    inline void setPalette_helper(const QPalette &palette) {
        if (resolvedPalette.resolve() == palette.resolve() && resolvedPalette == palette)
            return;
        updatePalette(palette);
    }

    void textChanged(const QString &text);

#if QT_CONFIG(accessibility)
    void accessibilityActiveChanged(bool active) override;
    QAccessible::Role accessibleRole() const override;
    void maybeSetAccessibleName(const QString &name);
#endif

    void cancelBackground();
    void executeBackground(bool complete = false);

    void itemGeometryChanged(QQuickItem *item, QQuickGeometryChange change, const QRectF &diff) override;
    void itemImplicitWidthChanged(QQuickItem *item) override;
    void itemImplicitHeightChanged(QQuickItem *item) override;
    void itemDestroyed(QQuickItem *item) override;

    struct ExtraData {
        bool hasTopInset = false;
        bool hasLeftInset = false;
        bool hasRightInset = false;
        bool hasBottomInset = false;
        bool hasBackgroundWidth = false;
        bool hasBackgroundHeight = false;
        qreal topInset = 0;
        qreal leftInset = 0;
        qreal rightInset = 0;
        qreal bottomInset = 0;
        QFont requestedFont;
        QPalette requestedPalette;
    };
    QLazilyAllocated<ExtraData> extra;

    bool resizingBackground = false;
    QPalette resolvedPalette;
    QQuickDeferredPointer<QQuickItem> background;
};

QT_END_NAMESPACE

#endif // QQUICKLABEL_P_P_H
