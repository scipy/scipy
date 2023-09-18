/****************************************************************************
**
** Copyright (C) 2015 Paul Lemire
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

#ifndef QT3DRENDER_QSTENCILTESTARGUMENTS_H
#define QT3DRENDER_QSTENCILTESTARGUMENTS_H

#include <QtCore/QObject>
#include <Qt3DRender/qt3drender_global.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QStencilTestPrivate;
class QStencilTestArgumentsPrivate;

class Q_3DRENDERSHARED_EXPORT QStencilTestArguments : public QObject
{
    Q_OBJECT
    Q_PROPERTY(StencilFaceMode faceMode READ faceMode NOTIFY faceModeChanged)
    Q_PROPERTY(uint comparisonMask READ comparisonMask WRITE setComparisonMask NOTIFY comparisonMaskChanged)
    Q_PROPERTY(int referenceValue READ referenceValue WRITE setReferenceValue NOTIFY referenceValueChanged)
    Q_PROPERTY(StencilFunction stencilFunction READ stencilFunction WRITE setStencilFunction NOTIFY stencilFunctionChanged)

public:
    enum StencilFaceMode
    {
        Front = 0x0404,
        Back = 0x0405,
        FrontAndBack = 0x0408
    };
    Q_ENUM(StencilFaceMode) // LCOV_EXCL_LINE

    enum StencilFunction
    {
        Never = 0x0200,
        Always = 0x0207,
        Less = 0x0201,
        LessOrEqual = 0x0203,
        Equal = 0x0202,
        GreaterOrEqual = 0x0206,
        Greater = 0x0204,
        NotEqual = 0x0205
    };
    Q_ENUM(StencilFunction) // LCOV_EXCL_LINE

    ~QStencilTestArguments();

    uint comparisonMask() const;
    int referenceValue() const;
    StencilFunction stencilFunction() const;

    StencilFaceMode faceMode() const;

public Q_SLOTS:
    void setComparisonMask(uint comparisonMask);
    void setReferenceValue(int referenceValue);
    void setStencilFunction(StencilFunction stencilFunction);

Q_SIGNALS:
    void comparisonMaskChanged(uint comparisonMask);
    void stencilFunctionChanged(StencilFunction stencilFunction);
    void referenceValueChanged(int referenceValue);
    void faceModeChanged(StencilFaceMode faceMode);

private:
    explicit QStencilTestArguments(StencilFaceMode face, QObject *parent = nullptr);

    friend class QStencilTestPrivate;

    Q_DECLARE_PRIVATE(QStencilTestArguments)
};

} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_QSTENCILTESTARGUMENTS_H
