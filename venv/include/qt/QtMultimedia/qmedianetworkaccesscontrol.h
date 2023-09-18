/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl-3.0.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or (at your option) the GNU General
** Public license version 3 or any later version approved by the KDE Free
** Qt Foundation. The licenses are as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL2 and LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-2.0.html and
** https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/


#ifndef QMEDIANETWORKACCESSCONTROL_H
#define QMEDIANETWORKACCESSCONTROL_H

#if 0
#pragma qt_class(QMediaNetworkAccessControl)
#endif

#include <QtMultimedia/qmediacontrol.h>

#include <QtCore/qlist.h>
#include <QtNetwork/qnetworkconfiguration.h>

#ifndef QT_NO_BEARERMANAGEMENT

QT_BEGIN_NAMESPACE

QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED

// Required for QDoc workaround
class QString;
class QT_DEPRECATED_VERSION_5_15 Q_MULTIMEDIA_EXPORT QMediaNetworkAccessControl : public QMediaControl
{
    Q_OBJECT
public:

    virtual ~QMediaNetworkAccessControl();

    virtual void setConfigurations(const QList<QNetworkConfiguration> &configuration) = 0;
    virtual QNetworkConfiguration currentConfiguration() const = 0;

Q_SIGNALS:
    void configurationChanged(const QNetworkConfiguration& configuration);

protected:
    explicit QMediaNetworkAccessControl(QObject *parent = nullptr);
};

#define QMediaNetworkAccessControl_iid "org.qt-project.qt.medianetworkaccesscontrol/5.0"
Q_MEDIA_DECLARE_CONTROL(QMediaNetworkAccessControl, QMediaNetworkAccessControl_iid)

QT_WARNING_POP

QT_END_NAMESPACE

#endif

#endif
