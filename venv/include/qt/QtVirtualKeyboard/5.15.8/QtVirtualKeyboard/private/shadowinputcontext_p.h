/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Virtual Keyboard module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:GPL$
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
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef SHADOWINPUTCONTEXT_P_H
#define SHADOWINPUTCONTEXT_P_H

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

#include <QObject>
#include <QPointer>
#include <QMetaType>
#include <QRectF>
#include <QtVirtualKeyboard/qvirtualkeyboard_global.h>

QT_BEGIN_NAMESPACE

class QVirtualKeyboardInputContext;
class QVirtualKeyboardInputContextPrivate;

namespace QtVirtualKeyboard {

class ShadowInputContextPrivate;

class QVIRTUALKEYBOARD_EXPORT ShadowInputContext : public QObject
{
    Q_OBJECT
    Q_DISABLE_COPY(ShadowInputContext)
    Q_DECLARE_PRIVATE(ShadowInputContext)
    Q_PROPERTY(QObject *inputItem READ inputItem WRITE setInputItem NOTIFY inputItemChanged)
    Q_PROPERTY(QRectF anchorRectangle READ anchorRectangle NOTIFY anchorRectangleChanged)
    Q_PROPERTY(QRectF cursorRectangle READ cursorRectangle NOTIFY cursorRectangleChanged)
    Q_PROPERTY(bool anchorRectIntersectsClipRect READ anchorRectIntersectsClipRect NOTIFY anchorRectIntersectsClipRectChanged)
    Q_PROPERTY(bool cursorRectIntersectsClipRect READ cursorRectIntersectsClipRect NOTIFY cursorRectIntersectsClipRectChanged)
    Q_PROPERTY(bool selectionControlVisible READ selectionControlVisible NOTIFY selectionControlVisibleChanged)

    explicit ShadowInputContext(QObject *parent = nullptr);

    void setInputContext(QVirtualKeyboardInputContext *inputContext);

public:
    QObject *inputItem() const;
    void setInputItem(QObject *inputItem);
    QRectF anchorRectangle() const;
    QRectF cursorRectangle() const;
    bool anchorRectIntersectsClipRect() const;
    bool cursorRectIntersectsClipRect() const;
    bool selectionControlVisible() const;

    Q_INVOKABLE void setSelectionOnFocusObject(const QPointF &anchorPos, const QPointF &cursorPos);
    Q_INVOKABLE void updateSelectionProperties();

signals:
    void inputItemChanged();
    void anchorRectangleChanged();
    void cursorRectangleChanged();
    void anchorRectIntersectsClipRectChanged();
    void cursorRectIntersectsClipRectChanged();
    void selectionControlVisibleChanged();

private:
    void update(Qt::InputMethodQueries queries);
    QVariant queryFocusObject(Qt::InputMethodQuery query, QVariant argument);

private:
    friend class ::QVirtualKeyboardInputContextPrivate;
    friend class ::QVirtualKeyboardInputContext;
};

} // namespace QtVirtualKeyboard
QT_END_NAMESPACE

#endif // SHADOWINPUTCONTEXT_P_H
