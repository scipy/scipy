/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Charts module of the Qt Toolkit.
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

//  W A R N I N G
//  -------------
//
// This file is not part of the Qt Chart API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#ifndef GLWIDGET_H
#define GLWIDGET_H

#ifndef QT_NO_OPENGL

#include <QtWidgets/QOpenGLWidget>
#include <QtWidgets/QGraphicsView>
#include <QtGui/QOpenGLFunctions>
#include <QtGui/QOpenGLVertexArrayObject>
#include <QtGui/QOpenGLBuffer>
#include <QtGui/QOpenGLFramebufferObject>
#include <QtCore/QHash>
#include <QtCharts/QAbstractSeries>
#include <QtCharts/QXYSeries>
#include <QtCharts/QChart>

QT_FORWARD_DECLARE_CLASS(QOpenGLShaderProgram)

QT_CHARTS_BEGIN_NAMESPACE

class GLXYSeriesDataManager;

class GLWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    GLWidget(GLXYSeriesDataManager *xyDataManager, QChart *chart, QGraphicsView *parent = 0);
    ~GLWidget();

    bool needsReset() const;

public Q_SLOTS:
    void cleanup();
    void cleanXYSeriesResources(const QXYSeries *series);

protected:
    void initializeGL() override;
    void paintGL() override;
    void resizeGL(int width, int height) override;
    void mouseDoubleClickEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;

private:
    QXYSeries *findSeriesAtEvent(QMouseEvent *event);
    void render(bool selection);
    void recreateSelectionFbo();
    QXYSeries *chartSeries(const QXYSeries *cSeries);

    QOpenGLShaderProgram *m_program;
    int m_shaderAttribLoc;
    int m_colorUniformLoc;
    int m_minUniformLoc;
    int m_deltaUniformLoc;
    int m_pointSizeUniformLoc;
    int m_matrixUniformLoc;
    QOpenGLVertexArrayObject m_vao;

    QHash<const QAbstractSeries *, QOpenGLBuffer *> m_seriesBufferMap;
    GLXYSeriesDataManager *m_xyDataManager;
    bool m_antiAlias;
    QGraphicsView *m_view;
    QOpenGLFramebufferObject *m_selectionFbo;
    QSize m_fboSize;
    QVector<const QXYSeries *> m_selectionVector;
    QChart *m_chart;
    bool m_recreateSelectionFbo;
    bool m_selectionRenderNeeded;
    QPoint m_mousePressPos;
    bool m_mousePressed;
    QXYSeries *m_lastPressSeries;
    QXYSeries *m_lastHoverSeries;
};

QT_CHARTS_END_NAMESPACE
#endif
#endif
