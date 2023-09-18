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

#ifndef QABSTRACTAXIS_H
#define QABSTRACTAXIS_H

#include <QtCharts/QChartGlobal>
#include <QtGui/QPen>
#include <QtGui/QFont>
#include <QtCore/QVariant>

QT_CHARTS_BEGIN_NAMESPACE

class QAbstractAxisPrivate;

class Q_CHARTS_EXPORT QAbstractAxis : public QObject
{
    Q_OBJECT
    //visibility
    Q_PROPERTY(bool visible READ isVisible WRITE setVisible NOTIFY visibleChanged)
    //arrow
    Q_PROPERTY(bool lineVisible READ isLineVisible WRITE setLineVisible NOTIFY lineVisibleChanged)
    Q_PROPERTY(QPen linePen READ linePen WRITE setLinePen NOTIFY linePenChanged)
    Q_PROPERTY(QColor color READ linePenColor WRITE setLinePenColor NOTIFY colorChanged)
    //labels
    Q_PROPERTY(bool labelsVisible READ labelsVisible WRITE setLabelsVisible NOTIFY labelsVisibleChanged)
    Q_PROPERTY(QBrush labelsBrush READ labelsBrush WRITE setLabelsBrush NOTIFY labelsBrushChanged)
    Q_PROPERTY(int labelsAngle READ labelsAngle WRITE setLabelsAngle NOTIFY labelsAngleChanged)
    Q_PROPERTY(QFont labelsFont READ labelsFont WRITE setLabelsFont NOTIFY labelsFontChanged)
    Q_PROPERTY(QColor labelsColor READ labelsColor WRITE setLabelsColor NOTIFY labelsColorChanged)
    //grid
    Q_PROPERTY(bool gridVisible READ isGridLineVisible WRITE setGridLineVisible NOTIFY gridVisibleChanged)
    Q_PROPERTY(QPen gridLinePen READ gridLinePen WRITE setGridLinePen NOTIFY gridLinePenChanged)
    Q_PROPERTY(bool minorGridVisible READ isMinorGridLineVisible WRITE setMinorGridLineVisible NOTIFY minorGridVisibleChanged)
    Q_PROPERTY(QPen minorGridLinePen READ minorGridLinePen WRITE setMinorGridLinePen NOTIFY minorGridLinePenChanged)
    Q_PROPERTY(QColor gridLineColor READ gridLineColor WRITE setGridLineColor NOTIFY gridLineColorChanged)
    Q_PROPERTY(QColor minorGridLineColor READ minorGridLineColor WRITE setMinorGridLineColor NOTIFY minorGridLineColorChanged)
    //shades
    Q_PROPERTY(bool shadesVisible READ shadesVisible WRITE setShadesVisible NOTIFY shadesVisibleChanged)
    Q_PROPERTY(QColor shadesColor READ shadesColor WRITE setShadesColor NOTIFY shadesColorChanged)
    Q_PROPERTY(QColor shadesBorderColor READ shadesBorderColor WRITE setShadesBorderColor NOTIFY shadesBorderColorChanged)
    Q_PROPERTY(QPen shadesPen READ shadesPen WRITE setShadesPen NOTIFY shadesPenChanged)
    Q_PROPERTY(QBrush shadesBrush READ shadesBrush WRITE setShadesBrush NOTIFY shadesBrushChanged)
    //title
    Q_PROPERTY(QString titleText READ titleText WRITE setTitleText NOTIFY titleTextChanged)
    Q_PROPERTY(QBrush titleBrush READ titleBrush WRITE setTitleBrush NOTIFY titleBrushChanged)
    Q_PROPERTY(bool titleVisible READ isTitleVisible WRITE setTitleVisible NOTIFY titleVisibleChanged)
    Q_PROPERTY(QFont titleFont READ titleFont WRITE setTitleFont NOTIFY titleFontChanged)
    //orientation
    Q_PROPERTY(Qt::Orientation orientation READ orientation)
    //aligment
    Q_PROPERTY(Qt::Alignment alignment READ alignment)
    Q_PROPERTY(bool reverse READ isReverse WRITE setReverse NOTIFY reverseChanged)

public:

    enum AxisType {
        AxisTypeNoAxis = 0x0,
        AxisTypeValue = 0x1,
        AxisTypeBarCategory = 0x2,
        AxisTypeCategory = 0x4,
        AxisTypeDateTime = 0x8,
        AxisTypeLogValue = 0x10
    };

    Q_DECLARE_FLAGS(AxisTypes, AxisType)

protected:
    explicit QAbstractAxis(QAbstractAxisPrivate &d, QObject *parent = nullptr);

public:
    ~QAbstractAxis();

    virtual AxisType type() const = 0;

    //visibility handling
    bool isVisible() const;
    void setVisible(bool visible = true);
    void show();
    void hide();

    //arrow handling
    bool isLineVisible() const;
    void setLineVisible(bool visible = true);
    void setLinePen(const QPen &pen);
    QPen linePen() const;
    void setLinePenColor(QColor color);
    QColor linePenColor() const;

    //grid handling
    bool isGridLineVisible() const;
    void setGridLineVisible(bool visible = true);
    void setGridLinePen(const QPen &pen);
    QPen gridLinePen() const;
    bool isMinorGridLineVisible() const;
    void setMinorGridLineVisible(bool visible = true);
    void setMinorGridLinePen(const QPen &pen);
    QPen minorGridLinePen() const;
    void setGridLineColor(const QColor &color);
    QColor gridLineColor();
    void setMinorGridLineColor(const QColor &color);
    QColor minorGridLineColor();

    //labels handling
    bool labelsVisible() const;
    void setLabelsVisible(bool visible = true);
    void setLabelsBrush(const QBrush &brush);
    QBrush labelsBrush() const;
    void setLabelsFont(const QFont &font);
    QFont labelsFont() const;
    void setLabelsAngle(int angle);
    int labelsAngle() const;
    void setLabelsColor(QColor color);
    QColor labelsColor() const;

    //title handling
    bool isTitleVisible() const;
    void setTitleVisible(bool visible = true);
    void setTitleBrush(const QBrush &brush);
    QBrush titleBrush() const;
    void setTitleFont(const QFont &font);
    QFont titleFont() const;
    void setTitleText(const QString &title);
    QString titleText() const;

    //shades handling
    bool shadesVisible() const;
    void setShadesVisible(bool visible = true);
    void setShadesPen(const QPen &pen);
    QPen shadesPen() const;
    void setShadesBrush(const QBrush &brush);
    QBrush shadesBrush() const;
    void setShadesColor(QColor color);
    QColor shadesColor() const;
    void setShadesBorderColor(QColor color);
    QColor shadesBorderColor() const;

    Qt::Orientation orientation() const;
    Qt::Alignment alignment() const;

    //range handling
    void setMin(const QVariant &min);
    void setMax(const QVariant &max);
    void setRange(const QVariant &min, const QVariant &max);

    //reverse handling
    void setReverse(bool reverse = true);
    bool isReverse() const;

    //label editable handling
    void setLabelsEditable(bool editable = true);
    bool labelsEditable() const;

Q_SIGNALS:
    void visibleChanged(bool visible);
    void linePenChanged(const QPen &pen);
    void lineVisibleChanged(bool visible);
    void labelsVisibleChanged(bool visible);
    void labelsBrushChanged(const QBrush &brush);
    void labelsFontChanged(const QFont &pen);
    void labelsAngleChanged(int angle);
    void gridLinePenChanged(const QPen &pen);
    void gridVisibleChanged(bool visible);
    void minorGridVisibleChanged(bool visible);
    void minorGridLinePenChanged(const QPen &pen);
    void gridLineColorChanged(const QColor &color);
    void minorGridLineColorChanged(const QColor &color);
    void colorChanged(QColor color);
    void labelsColorChanged(QColor color);
    void titleTextChanged(const QString &title);
    void titleBrushChanged(const QBrush &brush);
    void titleVisibleChanged(bool visible);
    void titleFontChanged(const QFont &font);
    void shadesVisibleChanged(bool visible);
    void shadesColorChanged(QColor color);
    void shadesBorderColorChanged(QColor color);
    void shadesPenChanged(const QPen &pen);
    void shadesBrushChanged(const QBrush &brush);
    void reverseChanged(bool reverse);
    void labelsEditableChanged(bool editable);

protected:
    QScopedPointer<QAbstractAxisPrivate> d_ptr;
    friend class ChartDataSet;
    friend class ChartPresenter;
    friend class ChartThemeManager;
    friend class AbstractDomain;
    friend class ChartAxisElement;
    friend class XYChart;

private:
    Q_DISABLE_COPY(QAbstractAxis)
};

QT_CHARTS_END_NAMESPACE

#endif // QABSTRACTAXIS_H
