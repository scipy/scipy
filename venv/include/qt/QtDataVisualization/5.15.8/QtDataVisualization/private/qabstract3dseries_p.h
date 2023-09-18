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

#ifndef QABSTRACT3DSERIES_P_H
#define QABSTRACT3DSERIES_P_H

#include "datavisualizationglobal_p.h"
#include "qabstract3dseries.h"

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class QAbstractDataProxy;
class Abstract3DController;

struct QAbstract3DSeriesChangeBitField {
    bool meshChanged                    : 1;
    bool meshSmoothChanged              : 1;
    bool meshRotationChanged            : 1;
    bool userDefinedMeshChanged         : 1;
    bool colorStyleChanged              : 1;
    bool baseColorChanged               : 1;
    bool baseGradientChanged            : 1;
    bool singleHighlightColorChanged    : 1;
    bool singleHighlightGradientChanged : 1;
    bool multiHighlightColorChanged     : 1;
    bool multiHighlightGradientChanged  : 1;
    bool nameChanged                    : 1;
    bool itemLabelChanged               : 1;
    bool itemLabelVisibilityChanged     : 1;
    bool visibilityChanged              : 1;

    QAbstract3DSeriesChangeBitField()
        : meshChanged(true),
          meshSmoothChanged(true),
          meshRotationChanged(true),
          userDefinedMeshChanged(true),
          colorStyleChanged(true),
          baseColorChanged(true),
          baseGradientChanged(true),
          singleHighlightColorChanged(true),
          singleHighlightGradientChanged(true),
          multiHighlightColorChanged(true),
          multiHighlightGradientChanged(true),
          nameChanged(true),
          itemLabelChanged(true),
          itemLabelVisibilityChanged(true),
          visibilityChanged(true)
    {
    }
};

struct QAbstract3DSeriesThemeOverrideBitField {
    bool colorStyleOverride              : 1;
    bool baseColorOverride               : 1;
    bool baseGradientOverride            : 1;
    bool singleHighlightColorOverride    : 1;
    bool singleHighlightGradientOverride : 1;
    bool multiHighlightColorOverride     : 1;
    bool multiHighlightGradientOverride  : 1;

    QAbstract3DSeriesThemeOverrideBitField()
        : colorStyleOverride(false),
          baseColorOverride(false),
          baseGradientOverride(false),
          singleHighlightColorOverride(false),
          singleHighlightGradientOverride(false),
          multiHighlightColorOverride(false),
          multiHighlightGradientOverride(false)
    {
    }
};

class QAbstract3DSeriesPrivate : public QObject
{
    Q_OBJECT
public:
    QAbstract3DSeriesPrivate(QAbstract3DSeries *q, QAbstract3DSeries::SeriesType type);
    virtual ~QAbstract3DSeriesPrivate();

    QAbstractDataProxy *dataProxy() const;
    virtual void setDataProxy(QAbstractDataProxy *proxy);
    virtual void setController(Abstract3DController *controller);
    virtual void connectControllerAndProxy(Abstract3DController *newController) = 0;
    virtual void createItemLabel() = 0;

    void setItemLabelFormat(const QString &format);
    void setVisible(bool visible);
    void setMesh(QAbstract3DSeries::Mesh mesh);
    void setMeshSmooth(bool enable);
    void setMeshRotation(const QQuaternion &rotation);
    void setUserDefinedMesh(const QString &meshFile);

    void setColorStyle(Q3DTheme::ColorStyle style);
    void setBaseColor(const QColor &color);
    void setBaseGradient(const QLinearGradient &gradient);
    void setSingleHighlightColor(const QColor &color);
    void setSingleHighlightGradient(const QLinearGradient &gradient);
    void setMultiHighlightColor(const QColor &color);
    void setMultiHighlightGradient(const QLinearGradient &gradient);
    void setName(const QString &name);

    void resetToTheme(const Q3DTheme &theme, int seriesIndex, bool force);
    QString itemLabel();
    void markItemLabelDirty();
    inline bool itemLabelDirty() const { return m_itemLabelDirty; }
    void setItemLabelVisible(bool visible);

    QAbstract3DSeriesChangeBitField m_changeTracker;
    QAbstract3DSeriesThemeOverrideBitField m_themeTracker;
    QAbstract3DSeries *q_ptr;
    QAbstract3DSeries::SeriesType m_type;
    QString m_itemLabelFormat;
    QAbstractDataProxy *m_dataProxy;
    bool m_visible;
    Abstract3DController *m_controller;
    QAbstract3DSeries::Mesh m_mesh;
    bool m_meshSmooth;
    QQuaternion m_meshRotation;
    QString m_userDefinedMesh;

    Q3DTheme::ColorStyle m_colorStyle;
    QColor m_baseColor;
    QLinearGradient m_baseGradient;
    QColor m_singleHighlightColor;
    QLinearGradient m_singleHighlightGradient;
    QColor m_multiHighlightColor;
    QLinearGradient m_multiHighlightGradient;

    QString m_name;
    QString m_itemLabel;
    bool m_itemLabelDirty;
    bool m_itemLabelVisible;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
