/****************************************************************************
**
** Copyright (C) 2017 Klaralvdalens Datakonsult AB (KDAB).
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** BSD License Usage
** Alternatively, you may use this file under the terms of the BSD license
** as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of The Qt Company Ltd nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QT3DEXTRAS_QEXTRUDEDTEXTGEOMETRY_H
#define QT3DEXTRAS_QEXTRUDEDTEXTGEOMETRY_H

#include <Qt3DExtras/qt3dextras_global.h>
#include <Qt3DRender/qgeometry.h>
#include <QtCore/QString>
#include <QtGui/QFont>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QAttribute;

} // namespace Qt3DRender

namespace Qt3DExtras {

class QExtrudedTextGeometryPrivate;

class Q_3DEXTRASSHARED_EXPORT QExtrudedTextGeometry : public Qt3DRender::QGeometry
{
    Q_OBJECT
    Q_PROPERTY(QString text READ text WRITE setText NOTIFY textChanged)
    Q_PROPERTY(QFont font READ font WRITE setFont NOTIFY fontChanged)
    Q_PROPERTY(float extrusionLength READ extrusionLength WRITE setDepth NOTIFY depthChanged)
    Q_PROPERTY(Qt3DRender::QAttribute *positionAttribute READ positionAttribute CONSTANT)
    Q_PROPERTY(Qt3DRender::QAttribute *normalAttribute READ normalAttribute CONSTANT)
    Q_PROPERTY(Qt3DRender::QAttribute *indexAttribute READ indexAttribute CONSTANT)

public:
    explicit QExtrudedTextGeometry(Qt3DCore::QNode *parent = nullptr);
    ~QExtrudedTextGeometry();

    Qt3DRender::QAttribute *positionAttribute() const;
    Qt3DRender::QAttribute *normalAttribute() const;
    Qt3DRender::QAttribute *indexAttribute() const;
    QString text() const;
    QFont font() const;
    float extrusionLength() const;

public Q_SLOTS:
    void setText(const QString &text);
    void setFont(const QFont &font);
    void setDepth(float extrusionLength);

Q_SIGNALS:
    void textChanged(const QString &text);
    void fontChanged(const QFont &font);
    void depthChanged(float extrusionLength);

protected:
    explicit QExtrudedTextGeometry(QExtrudedTextGeometryPrivate &dd, QNode *parent = nullptr);

private:
    Q_DECLARE_PRIVATE(QExtrudedTextGeometry)
};

} // namespace Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DEXTRAS_QEXTRUDEDTEXTGEOMETRY_H
