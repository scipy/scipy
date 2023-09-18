/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
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

#ifndef QVIRTUALKEYBOARDTRACE_H
#define QVIRTUALKEYBOARDTRACE_H

#include <QObject>
#include <QVariant>
#include <QPointF>
#include <QtVirtualKeyboard/qvirtualkeyboard_global.h>

QT_BEGIN_NAMESPACE

class QVirtualKeyboardTracePrivate;

class QVIRTUALKEYBOARD_EXPORT QVirtualKeyboardTrace : public QObject
{
    Q_OBJECT
    Q_DECLARE_PRIVATE(QVirtualKeyboardTrace)
    Q_PROPERTY(int traceId READ traceId WRITE setTraceId NOTIFY traceIdChanged)
    Q_PROPERTY(QStringList channels READ channels WRITE setChannels NOTIFY channelsChanged)
    Q_PROPERTY(int length READ length NOTIFY lengthChanged)
    Q_PROPERTY(bool final READ isFinal WRITE setFinal NOTIFY finalChanged)
    Q_PROPERTY(bool canceled READ isCanceled WRITE setCanceled NOTIFY canceledChanged)
    Q_PROPERTY(qreal opacity READ opacity WRITE setOpacity NOTIFY opacityChanged)
public:
    explicit QVirtualKeyboardTrace(QObject *parent = nullptr);
    ~QVirtualKeyboardTrace();

    int traceId() const;
    void setTraceId(int id);

    QStringList channels() const;
    void setChannels(const QStringList &channels);

    int length() const;

    Q_INVOKABLE QVariantList points(int pos = 0, int count = -1) const;
    Q_INVOKABLE int addPoint(const QPointF &point);

    Q_INVOKABLE void setChannelData(const QString &channel, int index, const QVariant &data);
    Q_INVOKABLE QVariantList channelData(const QString &channel, int pos = 0, int count = -1) const;

    bool isFinal() const;
    void setFinal(bool final);

    bool isCanceled() const;
    void setCanceled(bool canceled);

    qreal opacity() const;
    void setOpacity(qreal opacity);

Q_SIGNALS:
    void traceIdChanged(int traceId);
    void channelsChanged();
    void lengthChanged(int length);
    void finalChanged(bool isFinal);
    void canceledChanged(bool isCanceled);
    void opacityChanged(qreal opacity);
};

QT_END_NAMESPACE

#endif // QVIRTUALKEYBOARDTRACE_H
