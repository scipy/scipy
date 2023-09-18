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

#ifndef QABSTRACT3DSERIES_H
#define QABSTRACT3DSERIES_H

#include <QtDataVisualization/q3dtheme.h>
#include <QtCore/QObject>
#include <QtCore/QScopedPointer>
#include <QtCore/QString>
#include <QtGui/QLinearGradient>
#include <QtGui/QQuaternion>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class QAbstract3DSeriesPrivate;

class QT_DATAVISUALIZATION_EXPORT QAbstract3DSeries : public QObject
{
    Q_OBJECT
    Q_ENUMS(SeriesType)
    Q_ENUMS(Mesh)
    Q_PROPERTY(SeriesType type READ type CONSTANT)
    Q_PROPERTY(QString itemLabelFormat READ itemLabelFormat WRITE setItemLabelFormat NOTIFY itemLabelFormatChanged)
    Q_PROPERTY(bool visible READ isVisible WRITE setVisible NOTIFY visibilityChanged)
    Q_PROPERTY(Mesh mesh READ mesh WRITE setMesh NOTIFY meshChanged)
    Q_PROPERTY(bool meshSmooth READ isMeshSmooth WRITE setMeshSmooth NOTIFY meshSmoothChanged)
    Q_PROPERTY(QQuaternion meshRotation READ meshRotation WRITE setMeshRotation NOTIFY meshRotationChanged)
    Q_PROPERTY(QString userDefinedMesh READ userDefinedMesh WRITE setUserDefinedMesh NOTIFY userDefinedMeshChanged)
    Q_PROPERTY(QtDataVisualization::Q3DTheme::ColorStyle colorStyle READ colorStyle WRITE setColorStyle NOTIFY colorStyleChanged)
    Q_PROPERTY(QColor baseColor READ baseColor WRITE setBaseColor NOTIFY baseColorChanged)
    Q_PROPERTY(QLinearGradient baseGradient READ baseGradient WRITE setBaseGradient NOTIFY baseGradientChanged)
    Q_PROPERTY(QColor singleHighlightColor READ singleHighlightColor WRITE setSingleHighlightColor NOTIFY singleHighlightColorChanged)
    Q_PROPERTY(QLinearGradient singleHighlightGradient READ singleHighlightGradient WRITE setSingleHighlightGradient NOTIFY singleHighlightGradientChanged)
    Q_PROPERTY(QColor multiHighlightColor READ multiHighlightColor WRITE setMultiHighlightColor NOTIFY multiHighlightColorChanged)
    Q_PROPERTY(QLinearGradient multiHighlightGradient READ multiHighlightGradient WRITE setMultiHighlightGradient NOTIFY multiHighlightGradientChanged)
    Q_PROPERTY(QString name READ name WRITE setName NOTIFY nameChanged)
    Q_PROPERTY(QString itemLabel READ itemLabel NOTIFY itemLabelChanged REVISION 1)
    Q_PROPERTY(bool itemLabelVisible READ isItemLabelVisible WRITE setItemLabelVisible NOTIFY itemLabelVisibilityChanged REVISION 1)

public:
    enum SeriesType {
        SeriesTypeNone = 0,
        SeriesTypeBar = 1,
        SeriesTypeScatter = 2,
        SeriesTypeSurface = 4
    };

    enum Mesh {
        MeshUserDefined = 0,
        MeshBar,
        MeshCube,
        MeshPyramid,
        MeshCone,
        MeshCylinder,
        MeshBevelBar,
        MeshBevelCube,
        MeshSphere,
        MeshMinimal,
        MeshArrow,
        MeshPoint
    };

protected:
    explicit QAbstract3DSeries(QAbstract3DSeriesPrivate *d, QObject *parent = nullptr);

public:
    virtual ~QAbstract3DSeries();

    SeriesType type() const;

    void setItemLabelFormat(const QString &format);
    QString itemLabelFormat() const;

    void setVisible(bool visible);
    bool isVisible() const;

    void setMesh(Mesh mesh);
    Mesh mesh() const;

    void setMeshSmooth(bool enable);
    bool isMeshSmooth() const;

    void setMeshRotation(const QQuaternion &rotation);
    QQuaternion meshRotation() const;
    Q_INVOKABLE void setMeshAxisAndAngle(const QVector3D &axis, float angle);

    void setUserDefinedMesh(const QString &fileName);
    QString userDefinedMesh() const;

    void setColorStyle(Q3DTheme::ColorStyle style);
    Q3DTheme::ColorStyle colorStyle() const;
    void setBaseColor(const QColor &color);
    QColor baseColor() const;
    void setBaseGradient(const QLinearGradient &gradient);
    QLinearGradient baseGradient() const;
    void setSingleHighlightColor(const QColor &color);
    QColor singleHighlightColor() const;
    void setSingleHighlightGradient(const QLinearGradient &gradient);
    QLinearGradient singleHighlightGradient() const;
    void setMultiHighlightColor(const QColor &color);
    QColor multiHighlightColor() const;
    void setMultiHighlightGradient(const QLinearGradient &gradient);
    QLinearGradient multiHighlightGradient() const;

    void setName(const QString &name);
    QString name() const;

    QString itemLabel() const;
    void setItemLabelVisible(bool visible);
    bool isItemLabelVisible() const;

Q_SIGNALS:
    void itemLabelFormatChanged(const QString &format);
    void visibilityChanged(bool visible);
    void meshChanged(QAbstract3DSeries::Mesh mesh);
    void meshSmoothChanged(bool enabled);
    void meshRotationChanged(const QQuaternion &rotation);
    void userDefinedMeshChanged(const QString &fileName);
    void colorStyleChanged(Q3DTheme::ColorStyle style);
    void baseColorChanged(const QColor &color);
    void baseGradientChanged(const QLinearGradient &gradient);
    void singleHighlightColorChanged(const QColor &color);
    void singleHighlightGradientChanged(const QLinearGradient &gradient);
    void multiHighlightColorChanged(const QColor &color);
    void multiHighlightGradientChanged(const QLinearGradient &gradient);
    void nameChanged(const QString &name);
    Q_REVISION(1) void itemLabelChanged(const QString &label);
    Q_REVISION(1) void itemLabelVisibilityChanged(bool visible);

protected:
    QScopedPointer<QAbstract3DSeriesPrivate> d_ptr;

private:
    Q_DISABLE_COPY(QAbstract3DSeries)

    friend class Abstract3DController;
    friend class Bars3DController;
    friend class Surface3DController;
    friend class Surface3DRenderer;
    friend class Scatter3DRenderer;
    friend class Scatter3DController;
    friend class QBar3DSeries;
    friend class SeriesRenderCache;
    friend class Abstract3DRenderer;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
