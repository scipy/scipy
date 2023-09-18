/****************************************************************************
**
** Copyright (C) 2017 QNX Software Systems. All rights reserved.
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

#ifndef QCUSTOMAUDIOROLECONTROL_H
#define QCUSTOMAUDIOROLECONTROL_H

#include <QtMultimedia/qmediacontrol.h>
#include <QtMultimedia/qaudio.h>

QT_BEGIN_NAMESPACE

// Class forward declaration required for QDoc bug
class QString;
class QStringList;

class Q_MULTIMEDIA_EXPORT QCustomAudioRoleControl : public QMediaControl
{
    Q_OBJECT

public:
    virtual ~QCustomAudioRoleControl();

    virtual QString customAudioRole() const = 0;
    virtual void setCustomAudioRole(const QString &role) = 0;

    virtual QStringList supportedCustomAudioRoles() const = 0;

Q_SIGNALS:
    void customAudioRoleChanged(const QString &role);

protected:
    explicit QCustomAudioRoleControl(QObject *parent = nullptr);
};

#define QCustomAudioRoleControl_iid "org.qt-project.qt.customaudiorolecontrol/5.11"
Q_MEDIA_DECLARE_CONTROL(QCustomAudioRoleControl, QCustomAudioRoleControl_iid)

QT_END_NAMESPACE

#endif // QCUSTOMAUDIOROLECONTROL_H
