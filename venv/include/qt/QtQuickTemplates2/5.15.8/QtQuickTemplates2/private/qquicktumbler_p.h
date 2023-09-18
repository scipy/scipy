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

#ifndef QQUICKTUMBLER_P_H
#define QQUICKTUMBLER_P_H

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

#include <QtCore/qvariant.h>
#include <QtQml/qqmlcomponent.h>
#include <QtQuickTemplates2/private/qquickcontrol_p.h>

QT_BEGIN_NAMESPACE

class QQuickTumblerAttached;
class QQuickTumblerPrivate;

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickTumbler : public QQuickControl
{
    Q_OBJECT
    Q_PROPERTY(QVariant model READ model WRITE setModel NOTIFY modelChanged FINAL)
    Q_PROPERTY(int count READ count NOTIFY countChanged FINAL)
    Q_PROPERTY(int currentIndex READ currentIndex WRITE setCurrentIndex NOTIFY currentIndexChanged FINAL)
    Q_PROPERTY(QQuickItem *currentItem READ currentItem NOTIFY currentItemChanged FINAL)
    Q_PROPERTY(QQmlComponent *delegate READ delegate WRITE setDelegate NOTIFY delegateChanged FINAL)
    Q_PROPERTY(int visibleItemCount READ visibleItemCount WRITE setVisibleItemCount NOTIFY visibleItemCountChanged FINAL)
    // 2.1 (Qt 5.8)
    Q_PROPERTY(bool wrap READ wrap WRITE setWrap RESET resetWrap NOTIFY wrapChanged FINAL REVISION 1)
    // 2.2 (Qt 5.9)
    Q_PROPERTY(bool moving READ isMoving NOTIFY movingChanged FINAL REVISION 2)

public:
    explicit QQuickTumbler(QQuickItem *parent = nullptr);
    ~QQuickTumbler();

    QVariant model() const;
    void setModel(const QVariant &model);

    int count() const;

    int currentIndex() const;
    void setCurrentIndex(int currentIndex);
    QQuickItem *currentItem() const;

    QQmlComponent *delegate() const;
    void setDelegate(QQmlComponent *delegate);

    int visibleItemCount() const;
    void setVisibleItemCount(int visibleItemCount);

    static QQuickTumblerAttached *qmlAttachedProperties(QObject *object);

    // 2.1 (Qt 5.8)
    bool wrap() const;
    void setWrap(bool wrap);
    void resetWrap();

    // 2.2 (Qt 5.9)
    bool isMoving() const;

    enum PositionMode {
        Beginning,
        Center,
        End,
        Visible, // ListView-only
        Contain,
        SnapPosition
    };
    Q_ENUM(PositionMode)

    // 2.5 (Qt 5.12)
    Q_REVISION(5) Q_INVOKABLE void positionViewAtIndex(int index, PositionMode mode);

Q_SIGNALS:
    void modelChanged();
    void countChanged();
    void currentIndexChanged();
    void currentItemChanged();
    void delegateChanged();
    void visibleItemCountChanged();
    // 2.1 (Qt 5.8)
    Q_REVISION(1) void wrapChanged();
    // 2.2 (Qt 5.9)
    Q_REVISION(2) void movingChanged();

protected:
    void geometryChanged(const QRectF &newGeometry, const QRectF &oldGeometry) override;
    void componentComplete() override;
    void contentItemChange(QQuickItem *newItem, QQuickItem *oldItem) override;
    void keyPressEvent(QKeyEvent *event) override;
    void updatePolish() override;

    QFont defaultFont() const override;
    QPalette defaultPalette() const override;

private:
    Q_DISABLE_COPY(QQuickTumbler)
    Q_DECLARE_PRIVATE(QQuickTumbler)

    Q_PRIVATE_SLOT(d_func(), void _q_updateItemWidths())
    Q_PRIVATE_SLOT(d_func(), void _q_updateItemHeights())
    Q_PRIVATE_SLOT(d_func(), void _q_onViewCurrentIndexChanged())
    Q_PRIVATE_SLOT(d_func(), void _q_onViewCountChanged())
    Q_PRIVATE_SLOT(d_func(), void _q_onViewOffsetChanged())
    Q_PRIVATE_SLOT(d_func(), void _q_onViewContentYChanged())
};

class QQuickTumblerAttachedPrivate;

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickTumblerAttached : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QQuickTumbler *tumbler READ tumbler CONSTANT FINAL)
    Q_PROPERTY(qreal displacement READ displacement NOTIFY displacementChanged FINAL)

public:
    explicit QQuickTumblerAttached(QObject *parent = nullptr);

    QQuickTumbler *tumbler() const;
    qreal displacement() const;

Q_SIGNALS:
    void displacementChanged();

private:
    Q_DISABLE_COPY(QQuickTumblerAttached)
    Q_DECLARE_PRIVATE(QQuickTumblerAttached)
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QQuickTumbler)
QML_DECLARE_TYPEINFO(QQuickTumbler, QML_HAS_ATTACHED_PROPERTIES)

#endif // QQUICKTUMBLER_P_H
