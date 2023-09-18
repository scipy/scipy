/****************************************************************************
**
** Copyright (C) 2014 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DRENDER_QSHADERPROGRAM_H
#define QT3DRENDER_QSHADERPROGRAM_H

#include <Qt3DCore/qnode.h>
#include <Qt3DRender/qt3drender_global.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QShaderProgramPrivate;

class Q_3DRENDERSHARED_EXPORT QShaderProgram : public Qt3DCore::QNode
{
    Q_OBJECT
    Q_PROPERTY(QByteArray vertexShaderCode READ vertexShaderCode WRITE setVertexShaderCode NOTIFY vertexShaderCodeChanged)
    Q_PROPERTY(QByteArray tessellationControlShaderCode READ tessellationControlShaderCode WRITE setTessellationControlShaderCode NOTIFY tessellationControlShaderCodeChanged)
    Q_PROPERTY(QByteArray tessellationEvaluationShaderCode READ tessellationEvaluationShaderCode WRITE setTessellationEvaluationShaderCode NOTIFY tessellationEvaluationShaderCodeChanged)
    Q_PROPERTY(QByteArray geometryShaderCode READ geometryShaderCode WRITE setGeometryShaderCode NOTIFY geometryShaderCodeChanged)
    Q_PROPERTY(QByteArray fragmentShaderCode READ fragmentShaderCode WRITE setFragmentShaderCode NOTIFY fragmentShaderCodeChanged)
    Q_PROPERTY(QByteArray computeShaderCode READ computeShaderCode WRITE setComputeShaderCode NOTIFY computeShaderCodeChanged)
    Q_PROPERTY(QString log READ log NOTIFY logChanged REVISION 9)
    Q_PROPERTY(Status status READ status NOTIFY statusChanged REVISION 9)
    Q_PROPERTY(Format format READ format WRITE setFormat NOTIFY formatChanged REVISION 15)

public:
    explicit QShaderProgram(Qt3DCore::QNode *parent = nullptr);
    ~QShaderProgram();

    enum ShaderType {
        Vertex = 0,
        Fragment,
        TessellationControl,
        TessellationEvaluation,
        Geometry,
        Compute
    };
    Q_ENUM(ShaderType) // LCOV_EXCL_LINE

    enum Status {
        NotReady = 0,
        Ready,
        Error
    };
    Q_ENUM(Status) // LCOV_EXCL_LINE

    enum Format {
        GLSL = 0,
        SPIRV
    };
    Q_ENUM(Format) // LCOV_EXCL_LINE

    // Source code in-line
    QByteArray vertexShaderCode() const;
    QByteArray tessellationControlShaderCode() const;
    QByteArray tessellationEvaluationShaderCode() const;
    QByteArray geometryShaderCode() const;
    QByteArray fragmentShaderCode() const;
    QByteArray computeShaderCode() const;

    void setShaderCode(ShaderType type, const QByteArray &shaderCode);
    QByteArray shaderCode(ShaderType type) const;

    QString log() const;
    Status status() const;

    void setFormat(Format format);
    Format format() const;

    Q_INVOKABLE static QByteArray loadSource(const QUrl &sourceUrl);

public Q_SLOTS:
    void setVertexShaderCode(const QByteArray &vertexShaderCode);
    void setTessellationControlShaderCode(const QByteArray &tessellationControlShaderCode);
    void setTessellationEvaluationShaderCode(const QByteArray &tessellationEvaluationShaderCode);
    void setGeometryShaderCode(const QByteArray &geometryShaderCode);
    void setFragmentShaderCode(const QByteArray &fragmentShaderCode);
    void setComputeShaderCode(const QByteArray &computeShaderCode);

Q_SIGNALS:
    void vertexShaderCodeChanged(const QByteArray &vertexShaderCode);
    void tessellationControlShaderCodeChanged(const QByteArray &tessellationControlShaderCode);
    void tessellationEvaluationShaderCodeChanged(const QByteArray &tessellationEvaluationShaderCode);
    void geometryShaderCodeChanged(const QByteArray &geometryShaderCode);
    void fragmentShaderCodeChanged(const QByteArray &fragmentShaderCode);
    void computeShaderCodeChanged(const QByteArray &computeShaderCode);
    void logChanged(const QString &log);
    void statusChanged(Status status);
    void formatChanged(Format format);

protected:
    explicit QShaderProgram(QShaderProgramPrivate &dd, Qt3DCore::QNode *parent = nullptr);
    void sceneChangeEvent(const Qt3DCore::QSceneChangePtr &change) override;

private:
    Q_DECLARE_PRIVATE(QShaderProgram)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

}

QT_END_NAMESPACE

#endif // QT3DRENDER_QSHADERPROGRAM_H
