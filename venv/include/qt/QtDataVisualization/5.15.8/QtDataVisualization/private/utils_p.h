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

#ifndef UTILS_P_H
#define UTILS_P_H

#include "datavisualizationglobal_p.h"

QT_FORWARD_DECLARE_CLASS(QLinearGradient)

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class Utils
{
public:
    enum ParamType {
        ParamTypeUnknown = 0,
        ParamTypeInt,
        ParamTypeUInt,
        ParamTypeReal
    };

    static GLuint getNearestPowerOfTwo(GLuint value);
    static QVector4D vectorFromColor(const QColor &color);
    static QColor colorFromVector(const QVector3D &colorVector);
    static QColor colorFromVector(const QVector4D &colorVector);
    static QImage printTextToImage(const QFont &font,
                                   const QString &text,
                                   const QColor &bgrColor,
                                   const QColor &txtColor,
                                   bool labelBackground,
                                   bool borders = false,
                                   int maxLabelWidth = 0);
    static QVector4D getSelection(QPoint mousepos, int height);
    static QImage getGradientImage(QLinearGradient &gradient);

    static ParamType preParseFormat(const QString &format, QString &preStr, QString &postStr,
                                    int &precision, char &formatSpec);
    static QString formatLabelSprintf(const QByteArray &format, ParamType paramType, qreal value);
    static QString formatLabelLocalized(ParamType paramType, qreal value,
                               const QLocale &locale, const QString &preStr, const QString &postStr,
                               int precision, char formatSpec, const QByteArray &format);
    static QString defaultLabelFormat();

    static float wrapValue(float value, float min, float max);
    static QQuaternion calculateRotation(const QVector3D &xyzRotations);
    static bool isOpenGLES();
    static void resolveStatics();

private:
    static ParamType mapFormatCharToParamType(char formatSpec);
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
