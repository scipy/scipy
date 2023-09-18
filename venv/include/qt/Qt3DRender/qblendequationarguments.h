/****************************************************************************
**
** Copyright (C) 2014 Klaralvdalens Datakonsult AB (KDAB).
** Copyright (C) 2016 The Qt Company Ltd and/or its subsidiary(-ies).
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
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

#ifndef QT3DRENDER_QBLENDEQUATIONARGUMENTS_H
#define QT3DRENDER_QBLENDEQUATIONARGUMENTS_H

#include <Qt3DRender/qrenderstate.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QBlendEquationArgumentsPrivate;

class Q_3DRENDERSHARED_EXPORT QBlendEquationArguments : public QRenderState
{
    Q_OBJECT
    Q_PROPERTY(Blending sourceRgb READ sourceRgb WRITE setSourceRgb NOTIFY sourceRgbChanged)
    Q_PROPERTY(Blending sourceAlpha READ sourceAlpha WRITE setSourceAlpha NOTIFY sourceAlphaChanged)
    Q_PROPERTY(Blending destinationRgb READ destinationRgb WRITE setDestinationRgb NOTIFY destinationRgbChanged)
    Q_PROPERTY(Blending destinationAlpha READ destinationAlpha WRITE setDestinationAlpha NOTIFY destinationAlphaChanged)
    Q_PROPERTY(int bufferIndex READ bufferIndex WRITE setBufferIndex NOTIFY bufferIndexChanged)

public:

    enum Blending
    {
        Zero = 0,
        One = 1,
        SourceColor = 0x0300,
        SourceAlpha = 0x0302,
        Source1Alpha, // ### Qt 6: Fix -> has same value as OneMinusSourceAlpha
        Source1Color, // ### Qt 6: Fix -> has same value as DestinationAlpha
        DestinationColor = 0x0306,
        DestinationAlpha = 0x0304,
        SourceAlphaSaturate = 0x0308,
        ConstantColor = 0x8001,
        ConstantAlpha = 0x8003,
        OneMinusSourceColor = 0x0301,
        OneMinusSourceAlpha = 0x0303,
        OneMinusDestinationAlpha = 0x0305,
        OneMinusDestinationColor = 0x0307,
        OneMinusConstantColor = 0x8002,
        OneMinusConstantAlpha = 0x8004,
        OneMinusSource1Alpha,
        OneMinusSource1Color,
        OneMinusSource1Color0 = OneMinusSource1Color // ### Qt 6: Remove
    };
    Q_ENUM(Blending) // LCOV_EXCL_LINE

    explicit QBlendEquationArguments(Qt3DCore::QNode *parent = nullptr);
    ~QBlendEquationArguments();

    Blending sourceRgb() const;
    Blending destinationRgb() const;
    Blending sourceAlpha() const;
    Blending destinationAlpha() const;
    int bufferIndex() const;

public Q_SLOTS:
    void setSourceRgb(Blending sourceRgb);
    void setDestinationRgb(Blending destinationRgb);
    void setSourceAlpha(Blending sourceAlpha);
    void setDestinationAlpha(Blending destinationAlpha);
    void setSourceRgba(Blending sourceRgba);
    void setDestinationRgba(Blending destinationRgba);
    void setBufferIndex(int index);

Q_SIGNALS:
    void sourceRgbChanged(Blending sourceRgb);
    void sourceAlphaChanged(Blending sourceAlpha);
    void destinationRgbChanged(Blending destinationRgb);
    void destinationAlphaChanged(Blending destinationAlpha);
    void sourceRgbaChanged(Blending sourceRgba);
    void destinationRgbaChanged(Blending destinationRgba);
    void bufferIndexChanged(int index);

protected:
    explicit QBlendEquationArguments(QBlendEquationArgumentsPrivate &dd, Qt3DCore::QNode *parent = nullptr);

private:
    Q_DECLARE_PRIVATE(QBlendEquationArguments)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_QBLENDEQUATIONARGUMENTS_H
