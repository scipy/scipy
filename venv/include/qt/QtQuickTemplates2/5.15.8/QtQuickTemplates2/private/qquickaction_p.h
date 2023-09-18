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

#ifndef QQUICKACTION_P_H
#define QQUICKACTION_P_H

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

#include <QtQuickTemplates2/private/qtquicktemplates2global_p.h>
#include <QtQuickTemplates2/private/qquickicon_p.h>
#include <QtCore/qobject.h>
#include <QtQml/qqml.h>

QT_BEGIN_NAMESPACE

class QQuickIcon;
class QQuickActionPrivate;

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickAction : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QString text READ text WRITE setText NOTIFY textChanged FINAL)
    Q_PROPERTY(QQuickIcon icon READ icon WRITE setIcon NOTIFY iconChanged FINAL)
    Q_PROPERTY(bool enabled READ isEnabled WRITE setEnabled NOTIFY enabledChanged RESET resetEnabled FINAL)
    Q_PROPERTY(bool checked READ isChecked WRITE setChecked NOTIFY checkedChanged FINAL)
    Q_PROPERTY(bool checkable READ isCheckable WRITE setCheckable NOTIFY checkableChanged FINAL)
#if QT_CONFIG(shortcut)
    Q_PRIVATE_PROPERTY(QQuickAction::d_func(), QVariant shortcut READ shortcut WRITE setShortcut NOTIFY shortcutChanged FINAL)
#endif

public:
    explicit QQuickAction(QObject *parent = nullptr);
    ~QQuickAction();

    QString text() const;
    void setText(const QString &text);

    QQuickIcon icon() const;
    void setIcon(const QQuickIcon &icon);

    bool isEnabled() const;
    void setEnabled(bool enabled);
    void resetEnabled();

    bool isChecked() const;
    void setChecked(bool checked);

    bool isCheckable() const;
    void setCheckable(bool checkable);

#if QT_CONFIG(shortcut)
    QKeySequence shortcut() const;
    void setShortcut(const QKeySequence &shortcut);
#endif

public Q_SLOTS:
    void toggle(QObject *source = nullptr);
    void trigger(QObject *source = nullptr);

Q_SIGNALS:
    void textChanged(const QString &text);
    void iconChanged(const QQuickIcon &icon);
    void enabledChanged(bool enabled);
    void checkedChanged(bool checked);
    void checkableChanged(bool checkable);
#if QT_CONFIG(shortcut)
    void shortcutChanged(const QKeySequence &shortcut);
#endif

    void toggled(QObject *source = nullptr);
    void triggered(QObject *source = nullptr);

protected:
    bool event(QEvent *event) override;
    bool eventFilter(QObject *object, QEvent *event) override;

private:
    Q_DISABLE_COPY(QQuickAction)
    Q_DECLARE_PRIVATE(QQuickAction)
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QQuickAction)

#endif // QQUICKACTION_P_H
