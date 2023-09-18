/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
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

#ifndef QRENDERCAPTURE_H
#define QRENDERCAPTURE_H

#include <Qt3DRender/QFrameGraphNode>
#include <QtGui/QImage>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QRenderCapturePrivate;
class QRenderCaptureReplyPrivate;

class Q_3DRENDERSHARED_EXPORT QRenderCaptureReply : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QImage image READ image CONSTANT)
    Q_PROPERTY(int captureId READ captureId CONSTANT)
    Q_PROPERTY(bool complete READ isComplete NOTIFY completed)

public:

    QImage image() const;
    Q_DECL_DEPRECATED int captureId() const;
    bool isComplete() const;

    Q_INVOKABLE bool saveImage(const QString &fileName) const;
#if QT_DEPRECATED_SINCE(5, 9)
    // ### Qt 6: remove this
    Q_DECL_DEPRECATED_X("Use saveImage instead") Q_INVOKABLE void saveToFile(const QString &fileName) const;
#endif

Q_SIGNALS:
    Q_DECL_DEPRECATED_X("Use completed instead") void completeChanged(bool isComplete);
    void completed();

private:
    Q_DECLARE_PRIVATE(QRenderCaptureReply)

    QRenderCaptureReply(QObject *parent = nullptr);

    friend class QRenderCapturePrivate;
};

class Q_3DRENDERSHARED_EXPORT QRenderCapture : public QFrameGraphNode
{
    Q_OBJECT
public:
    explicit QRenderCapture(Qt3DCore::QNode *parent = nullptr);

    Q_INVOKABLE Q_DECL_DEPRECATED_X("Use the overload with no id parameter")
    Qt3DRender::QRenderCaptureReply *requestCapture(int captureId);
    Q_REVISION(9) Q_INVOKABLE Qt3DRender::QRenderCaptureReply *requestCapture();
    Q_REVISION(10) Q_INVOKABLE Qt3DRender::QRenderCaptureReply *requestCapture(const QRect &rect);

protected:
    void sceneChangeEvent(const Qt3DCore::QSceneChangePtr &change) override;

private:
    Q_DECLARE_PRIVATE(QRenderCapture)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

} // Qt3DRender

QT_END_NAMESPACE

#endif // QRENDERCAPTURE_H
