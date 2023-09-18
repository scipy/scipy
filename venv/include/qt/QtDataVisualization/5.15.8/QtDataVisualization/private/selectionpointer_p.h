/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Data Visualization module of the Qt Toolkit.
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

//
//  W A R N I N G
//  -------------
//
// This file is not part of the QtDataVisualization API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#ifndef SELECTIONPOINTER_P_H
#define SELECTIONPOINTER_P_H

#include "datavisualizationglobal_p.h"
#include "surface3dcontroller_p.h"

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class ShaderHelper;
class ObjectHelper;
class TextureHelper;
class Drawer;

class SelectionPointer : public QObject, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    explicit SelectionPointer(Drawer *drawer);
    ~SelectionPointer();

    void renderSelectionPointer(GLuint defaultFboHandle = 0, bool useOrtho = false);
    void renderSelectionLabel(GLuint defaultFboHandle = 0, bool useOrtho = false);
    void setPosition(const QVector3D &position);
    void setLabel(const QString &label, bool themeChange = false);
    void setPointerObject(ObjectHelper *object);
    void setLabelObject(ObjectHelper *object);
    void handleDrawerChange();
    void updateBoundingRect(const QRect &rect);
    void updateScene(Q3DScene *scene);
    void updateSliceData(bool sliceActivated, GLfloat autoScaleAdjustment);
    void setHighlightColor(const QVector4D &colorVector);
    void setRotation(const QQuaternion &rotation);

private:
    void initializeOpenGL();
    void initShaders();

private:
    ShaderHelper *m_labelShader;
    ShaderHelper *m_pointShader;
    ObjectHelper *m_labelObj; // Not owned
    ObjectHelper *m_pointObj; // Not owned
    TextureHelper *m_textureHelper;
    Q3DTheme *m_cachedTheme;
    bool m_labelBackground;
    LabelItem m_labelItem;
    Drawer *m_drawer;
    QRect m_mainViewPort;
    QVector3D m_position;
    Q3DScene *m_cachedScene;
    QString m_label;
    bool m_cachedIsSlicingActivated;
    GLfloat m_autoScaleAdjustment;
    QVector4D m_highlightColor;
    QQuaternion m_rotation;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
