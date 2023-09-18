/****************************************************************************
**
** Copyright (C) 2017 Klaralvdalens Datakonsult AB (KDAB).
** Contact: http://www.qt-project.org/legal
**
** This file is part of the Qt3D module of the Qt Toolkit.
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

#ifndef QT3DANIMATION_QANIMATIONCLIPLOADER_H
#define QT3DANIMATION_QANIMATIONCLIPLOADER_H

#include <Qt3DAnimation/qt3danimation_global.h>
#include <Qt3DAnimation/qabstractanimationclip.h>
#include <QtCore/qurl.h>

QT_BEGIN_NAMESPACE

namespace Qt3DAnimation {

class QAnimationClipLoaderPrivate;

class Q_3DANIMATIONSHARED_EXPORT QAnimationClipLoader : public QAbstractAnimationClip
{
    Q_OBJECT
    Q_PROPERTY(QUrl source READ source WRITE setSource NOTIFY sourceChanged)
    Q_PROPERTY(Status status READ status NOTIFY statusChanged)

public:
    explicit QAnimationClipLoader(Qt3DCore::QNode *parent = nullptr);
    explicit QAnimationClipLoader(const QUrl &source,
                                  Qt3DCore::QNode *parent = nullptr);
    ~QAnimationClipLoader();

    enum Status {
        NotReady = 0,
        Ready,
        Error
    };
    Q_ENUM(Status) // LCOV_EXCL_LINE

    QUrl source() const;
    Status status() const;

public Q_SLOTS:
    void setSource(const QUrl &source);

Q_SIGNALS:
    void sourceChanged(const QUrl &source);
    void statusChanged(Status status);

protected:
    explicit QAnimationClipLoader(QAnimationClipLoaderPrivate &dd, Qt3DCore::QNode *parent = nullptr);
    // TODO Unused remove in Qt6
    void sceneChangeEvent(const Qt3DCore::QSceneChangePtr &change) override;

private:
    Q_DECLARE_PRIVATE(QAnimationClipLoader)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

} // namespace Qt3DAnimation

QT_END_NAMESPACE

#endif // QT3DANIMATION_QANIMATIONCLIPLOADER_H
