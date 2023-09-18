/****************************************************************************
**
** Copyright (C) 2017 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DRENDER_QSHADERPROGRAMBUILDER_H
#define QT3DRENDER_QSHADERPROGRAMBUILDER_H

#include <Qt3DCore/qnode.h>
#include <Qt3DRender/qt3drender_global.h>

#include <QtCore/qurl.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QShaderProgram;
class QShaderProgramBuilderPrivate;

class Q_3DRENDERSHARED_EXPORT QShaderProgramBuilder : public Qt3DCore::QNode
{
    Q_OBJECT
    Q_PROPERTY(Qt3DRender::QShaderProgram* shaderProgram READ shaderProgram WRITE setShaderProgram NOTIFY shaderProgramChanged)
    Q_PROPERTY(QStringList enabledLayers READ enabledLayers WRITE setEnabledLayers NOTIFY enabledLayersChanged)
    Q_PROPERTY(QUrl vertexShaderGraph READ vertexShaderGraph WRITE setVertexShaderGraph NOTIFY vertexShaderGraphChanged)
    Q_PROPERTY(QUrl tessellationControlShaderGraph READ tessellationControlShaderGraph WRITE setTessellationControlShaderGraph NOTIFY tessellationControlShaderGraphChanged)
    Q_PROPERTY(QUrl tessellationEvaluationShaderGraph READ tessellationEvaluationShaderGraph WRITE setTessellationEvaluationShaderGraph NOTIFY tessellationEvaluationShaderGraphChanged)
    Q_PROPERTY(QUrl geometryShaderGraph READ geometryShaderGraph WRITE setGeometryShaderGraph NOTIFY geometryShaderGraphChanged)
    Q_PROPERTY(QUrl fragmentShaderGraph READ fragmentShaderGraph WRITE setFragmentShaderGraph NOTIFY fragmentShaderGraphChanged)
    Q_PROPERTY(QUrl computeShaderGraph READ computeShaderGraph WRITE setComputeShaderGraph NOTIFY computeShaderGraphChanged)
    Q_PROPERTY(QByteArray vertexShaderCode READ vertexShaderCode NOTIFY vertexShaderCodeChanged REVISION 13)
    Q_PROPERTY(QByteArray tessellationControlShaderCode READ tessellationControlShaderCode NOTIFY tessellationControlShaderCodeChanged REVISION 13)
    Q_PROPERTY(QByteArray tessellationEvaluationShaderCode READ tessellationEvaluationShaderCode  NOTIFY tessellationEvaluationShaderCodeChanged REVISION 13)
    Q_PROPERTY(QByteArray geometryShaderCode READ geometryShaderCode NOTIFY geometryShaderCodeChanged REVISION 13)
    Q_PROPERTY(QByteArray fragmentShaderCode READ fragmentShaderCode NOTIFY fragmentShaderCodeChanged REVISION 13)
    Q_PROPERTY(QByteArray computeShaderCode READ computeShaderCode NOTIFY computeShaderCodeChanged REVISION 13)

public:
    explicit QShaderProgramBuilder(Qt3DCore::QNode *parent = nullptr);
    ~QShaderProgramBuilder();

    QShaderProgram *shaderProgram() const;
    QStringList enabledLayers() const;
    QUrl vertexShaderGraph() const;
    QUrl tessellationControlShaderGraph() const;
    QUrl tessellationEvaluationShaderGraph() const;
    QUrl geometryShaderGraph() const;
    QUrl fragmentShaderGraph() const;
    QUrl computeShaderGraph() const;
    QByteArray vertexShaderCode() const;
    QByteArray tessellationControlShaderCode() const;
    QByteArray tessellationEvaluationShaderCode() const;
    QByteArray geometryShaderCode() const;
    QByteArray fragmentShaderCode() const;
    QByteArray computeShaderCode() const;

public Q_SLOTS:
    void setShaderProgram(Qt3DRender::QShaderProgram *program);
    void setEnabledLayers(const QStringList &layers);
    void setVertexShaderGraph(const QUrl &vertexShaderGraph);
    void setTessellationControlShaderGraph(const QUrl &tessellationControlShaderGraph);
    void setTessellationEvaluationShaderGraph(const QUrl &tessellationEvaluationShaderGraph);
    void setGeometryShaderGraph(const QUrl &geometryShaderGraph);
    void setFragmentShaderGraph(const QUrl &fragmentShaderGraph);
    void setComputeShaderGraph(const QUrl &computeShaderGraph);

Q_SIGNALS:
    void shaderProgramChanged(Qt3DRender::QShaderProgram *shaderProgram);
    void enabledLayersChanged(const QStringList &layers);
    void vertexShaderGraphChanged(const QUrl &vertexShaderGraph);
    void tessellationControlShaderGraphChanged(const QUrl &tessellationControlShaderGraph);
    void tessellationEvaluationShaderGraphChanged(const QUrl &tessellationEvaluationShaderGraph);
    void geometryShaderGraphChanged(const QUrl &geometryShaderGraph);
    void fragmentShaderGraphChanged(const QUrl &fragmentShaderGraph);
    void computeShaderGraphChanged(const QUrl &computeShaderGraph);
    Q_REVISION(13) void vertexShaderCodeChanged(const QByteArray &vertexShaderCode);
    Q_REVISION(13) void tessellationControlShaderCodeChanged(const QByteArray &tessellationControlShaderCode);
    Q_REVISION(13) void tessellationEvaluationShaderCodeChanged(const QByteArray &tessellationEvaluationShaderCode);
    Q_REVISION(13) void geometryShaderCodeChanged(const QByteArray &geometryShaderCode);
    Q_REVISION(13) void fragmentShaderCodeChanged(const QByteArray &fragmentShaderCode);
    Q_REVISION(13) void computeShaderCodeChanged(const QByteArray &computeShaderCode);

protected:
    explicit QShaderProgramBuilder(QShaderProgramBuilderPrivate &dd, Qt3DCore::QNode *parent = nullptr);
    // TODO Unused remove in Qt6
    void sceneChangeEvent(const Qt3DCore::QSceneChangePtr &change) override;

private:
    Q_DECLARE_PRIVATE(QShaderProgramBuilder)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

}

QT_END_NAMESPACE

#endif // QT3DRENDER_QSHADERPROGRAMBUILDER_H
