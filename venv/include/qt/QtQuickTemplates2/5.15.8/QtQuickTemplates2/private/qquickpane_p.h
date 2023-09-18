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

#ifndef QQUICKPANE_P_H
#define QQUICKPANE_P_H

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
#include <QtQml/qqmllist.h>

QT_BEGIN_NAMESPACE

class QQuickPanePrivate;

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickPane : public QQuickControl
{
    Q_OBJECT
    Q_PROPERTY(qreal contentWidth READ contentWidth WRITE setContentWidth RESET resetContentWidth NOTIFY contentWidthChanged FINAL)
    Q_PROPERTY(qreal contentHeight READ contentHeight WRITE setContentHeight RESET resetContentHeight NOTIFY contentHeightChanged FINAL)
    Q_PRIVATE_PROPERTY(QQuickPane::d_func(), QQmlListProperty<QObject> contentData READ contentData FINAL)
    Q_PRIVATE_PROPERTY(QQuickPane::d_func(), QQmlListProperty<QQuickItem> contentChildren READ contentChildren NOTIFY contentChildrenChanged FINAL)
    Q_CLASSINFO("DefaultProperty", "contentData")

public:
    explicit QQuickPane(QQuickItem *parent = nullptr);
    ~QQuickPane();

    qreal contentWidth() const;
    void setContentWidth(qreal width);
    void resetContentWidth();

    qreal contentHeight() const;
    void setContentHeight(qreal height);
    void resetContentHeight();

Q_SIGNALS:
    void contentWidthChanged();
    void contentHeightChanged();
    void contentChildrenChanged();

protected:
    QQuickPane(QQuickPanePrivate &dd, QQuickItem *parent);

    void componentComplete() override;

    void contentItemChange(QQuickItem *newItem, QQuickItem *oldItem) override;
    virtual void contentSizeChange(const QSizeF &newSize, const QSizeF &oldSize);

#if QT_CONFIG(accessibility)
    QAccessible::Role accessibleRole() const override;
#endif

private:
    Q_DISABLE_COPY(QQuickPane)
    Q_DECLARE_PRIVATE(QQuickPane)
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QQuickPane)

#endif // QQUICKPANE_P_H
