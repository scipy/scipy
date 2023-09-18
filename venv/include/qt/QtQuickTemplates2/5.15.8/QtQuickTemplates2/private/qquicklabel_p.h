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

#ifndef QQUICKLABEL_P_H
#define QQUICKLABEL_P_H

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

#include <QtGui/qpalette.h>
#include <QtQuick/private/qquicktext_p.h>
#include <QtQuickTemplates2/private/qtquicktemplates2global_p.h>

QT_BEGIN_NAMESPACE

class QQuickLabelPrivate;

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickLabel : public QQuickText
{
    Q_OBJECT
    Q_PROPERTY(QFont font READ font WRITE setFont NOTIFY fontChanged) // override
    Q_PROPERTY(QQuickItem *background READ background WRITE setBackground NOTIFY backgroundChanged FINAL)
    // 2.3 (Qt 5.10)
    Q_PROPERTY(QPalette palette READ palette WRITE setPalette RESET resetPalette NOTIFY paletteChanged FINAL REVISION 3)
    // 2.5 (Qt 5.12)
    Q_PROPERTY(qreal implicitBackgroundWidth READ implicitBackgroundWidth NOTIFY implicitBackgroundWidthChanged FINAL REVISION 5)
    Q_PROPERTY(qreal implicitBackgroundHeight READ implicitBackgroundHeight NOTIFY implicitBackgroundHeightChanged FINAL REVISION 5)
    Q_PROPERTY(qreal topInset READ topInset WRITE setTopInset RESET resetTopInset NOTIFY topInsetChanged FINAL REVISION 5)
    Q_PROPERTY(qreal leftInset READ leftInset WRITE setLeftInset RESET resetLeftInset NOTIFY leftInsetChanged FINAL REVISION 5)
    Q_PROPERTY(qreal rightInset READ rightInset WRITE setRightInset RESET resetRightInset NOTIFY rightInsetChanged FINAL REVISION 5)
    Q_PROPERTY(qreal bottomInset READ bottomInset WRITE setBottomInset RESET resetBottomInset NOTIFY bottomInsetChanged FINAL REVISION 5)
    Q_CLASSINFO("DeferredPropertyNames", "background")

public:
    explicit QQuickLabel(QQuickItem *parent = nullptr);
    ~QQuickLabel();

    QFont font() const;
    void setFont(const QFont &font);

    QQuickItem *background() const;
    void setBackground(QQuickItem *background);

    // 2.3 (Qt 5.10)
    QPalette palette() const;
    void setPalette(const QPalette &palette);
    void resetPalette();

    // 2.5 (Qt 5.12)
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
    void backgroundChanged();
    // 2.3 (Qt 5.10)
    Q_REVISION(3) void paletteChanged();
    // 2.5 (Qt 5.12)
    Q_REVISION(5) void implicitBackgroundWidthChanged();
    Q_REVISION(5) void implicitBackgroundHeightChanged();
    Q_REVISION(5) void topInsetChanged();
    Q_REVISION(5) void leftInsetChanged();
    Q_REVISION(5) void rightInsetChanged();
    Q_REVISION(5) void bottomInsetChanged();

protected:
    void classBegin() override;
    void componentComplete() override;

    void itemChange(ItemChange change, const ItemChangeData &value) override;
    void geometryChanged(const QRectF &newGeometry, const QRectF &oldGeometry) override;
    virtual void insetChange(const QMarginsF &newInset, const QMarginsF &oldInset);

private:
    Q_DISABLE_COPY(QQuickLabel)
    Q_DECLARE_PRIVATE(QQuickLabel)
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QQuickLabel)

#endif // QQUICKLABEL_P_H
