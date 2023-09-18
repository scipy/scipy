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

#ifndef QQUICKCONTAINER_P_H
#define QQUICKCONTAINER_P_H

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

class QQuickContainerPrivate;

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickContainer : public QQuickControl
{
    Q_OBJECT
    Q_PROPERTY(int count READ count NOTIFY countChanged FINAL)
    Q_PROPERTY(QVariant contentModel READ contentModel CONSTANT FINAL)
    Q_PROPERTY(QQmlListProperty<QObject> contentData READ contentData FINAL)
    Q_PROPERTY(QQmlListProperty<QQuickItem> contentChildren READ contentChildren NOTIFY contentChildrenChanged FINAL)
    Q_PROPERTY(int currentIndex READ currentIndex WRITE setCurrentIndex NOTIFY currentIndexChanged FINAL)
    Q_PROPERTY(QQuickItem *currentItem READ currentItem NOTIFY currentItemChanged FINAL)
    // 2.5 (Qt 5.12)
    Q_PROPERTY(qreal contentWidth READ contentWidth WRITE setContentWidth RESET resetContentWidth NOTIFY contentWidthChanged FINAL REVISION 5)
    Q_PROPERTY(qreal contentHeight READ contentHeight WRITE setContentHeight RESET resetContentHeight NOTIFY contentHeightChanged FINAL REVISION 5)
    Q_CLASSINFO("DefaultProperty", "contentData")

public:
    explicit QQuickContainer(QQuickItem *parent = nullptr);
    ~QQuickContainer();

    int count() const;
    Q_INVOKABLE QQuickItem *itemAt(int index) const;
    Q_INVOKABLE void addItem(QQuickItem *item);
    Q_INVOKABLE void insertItem(int index, QQuickItem *item);
    Q_INVOKABLE void moveItem(int from, int to);
    Q_INVOKABLE void removeItem(const QVariant &item); // ### Qt 6: remove
    void removeItem(QQuickItem *item); // ### Qt 6: Q_INVOKABLE
    // 2.3 (Qt 5.10)
    Q_REVISION(3) Q_INVOKABLE QQuickItem *takeItem(int index);

    QVariant contentModel() const;
    QQmlListProperty<QObject> contentData();
    QQmlListProperty<QQuickItem> contentChildren();

    int currentIndex() const;
    QQuickItem *currentItem() const;

    // 2.5 (Qt 5.12)
    qreal contentWidth() const;
    void setContentWidth(qreal width);
    void resetContentWidth();

    qreal contentHeight() const;
    void setContentHeight(qreal height);
    void resetContentHeight();

public Q_SLOTS:
    void setCurrentIndex(int index);
    // 2.1 (Qt 5.8)
    Q_REVISION(1) void incrementCurrentIndex();
    Q_REVISION(1) void decrementCurrentIndex();

Q_SIGNALS:
    void countChanged();
    void contentChildrenChanged();
    void currentIndexChanged();
    void currentItemChanged();
    // 2.5 (Qt 5.12)
    Q_REVISION(5) void contentWidthChanged();
    Q_REVISION(5) void contentHeightChanged();

protected:
    QQuickContainer(QQuickContainerPrivate &dd, QQuickItem *parent);

    void componentComplete() override;

    void itemChange(ItemChange change, const ItemChangeData &data) override;
    void contentItemChange(QQuickItem *newItem, QQuickItem *oldItem) override;

    virtual bool isContent(QQuickItem *item) const;
    virtual void itemAdded(int index, QQuickItem *item);
    virtual void itemMoved(int index, QQuickItem *item);
    virtual void itemRemoved(int index, QQuickItem *item);

private:
    Q_DISABLE_COPY(QQuickContainer)
    Q_DECLARE_PRIVATE(QQuickContainer)
    Q_PRIVATE_SLOT(d_func(), void _q_currentIndexChanged())
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QQuickContainer)

#endif // QQUICKCONTAINER_P_H
