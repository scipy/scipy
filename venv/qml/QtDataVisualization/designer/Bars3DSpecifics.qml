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

import QtQuick 2.0
import HelperWidgets 2.0
import QtQuick.Layouts 1.0
import QtQuick.Controls 1.1 as Controls

Column {
    anchors.left: parent.left
    anchors.right: parent.right

    Section {
        anchors.left: parent.left
        anchors.right: parent.right
        caption: qsTr("Bars3D")

        SectionLayout {
            Label {
                text: qsTr("multiSeriesUniform")
                tooltip: qsTr("Multiseries Uniform")
                Layout.fillWidth: true
            }
            SecondColumnLayout {
                CheckBox {
                    backendValue: backendValues.multiSeriesUniform
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("barThickness")
                tooltip: qsTr("Bar Thickness Ratio")
                Layout.fillWidth: true
            }
            SecondColumnLayout {
                SpinBox {
                    backendValue: backendValues.barThickness
                    minimumValue: 0.01
                    maximumValue: 100.0
                    stepSize: 0.01
                    decimals: 2
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("barSpacing")
                tooltip: qsTr("Bar Spacing")
                Layout.fillWidth: true
            }
            SecondColumnLayout {
                SpinBox {
                    prefix: "col: "
                    backendValue: backendValues.barSpacing_width
                    minimumValue: 0.0
                    maximumValue: 10.0
                    stepSize: 0.01
                    decimals: 2
                    Layout.fillWidth: true
                }
                SpinBox {
                    prefix: "row: "
                    backendValue: backendValues.barSpacing_height
                    minimumValue: 0.0
                    maximumValue: 10.0
                    stepSize: 0.01
                    decimals: 2
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("barSpacingRelative")
                tooltip: qsTr("Bar Spacing Relative")
                Layout.fillWidth: true
            }
            SecondColumnLayout {
                CheckBox {
                    backendValue: backendValues.barSpacingRelative
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("renderingMode")
                tooltip: qsTr("Rendering Mode")
                Layout.fillWidth: true
            }
            SecondColumnLayout {
                ComboBox {
                    backendValue: backendValues.renderingMode
                    model: ["RenderIndirect", "RenderDirectToBackground",
                        "RenderDirectToBackground_NoClear"]
                    Layout.fillWidth: true
                    scope: "AbstractGraph3D"
                }
            }
            Label {
                text: qsTr("msaaSamples")
                tooltip: qsTr("MSAA Sample Count")
                Layout.fillWidth: true
            }
            SpinBox {
                suffix: " x MSAA"
                backendValue: backendValues.msaaSamples
                minimumValue: 0
                maximumValue: 16
                Layout.fillWidth: true
            }
            Label {
                text: qsTr("shadowQuality")
                tooltip: qsTr("Shadow Quality")
                Layout.fillWidth: true
            }
            SecondColumnLayout {
                ComboBox {
                    backendValue: backendValues.shadowQuality
                    model: ["ShadowQualityNone", "ShadowQualityLow", "ShadowQualityMedium",
                        "ShadowQualityHigh", "ShadowQualitySoftLow", "ShadowQualitySoftMedium",
                        "ShadowQualitySoftHigh"]
                    Layout.fillWidth: true
                    scope: "AbstractGraph3D"
                }
            }
            Label {
                text: qsTr("selectionMode")
                tooltip: qsTr("Selection Mode")
                Layout.fillWidth: true
            }
            SecondColumnLayout {
                id: selectionLayout
                property bool isInModel: backendValue.isInModel;
                property bool isInSubState: backendValue.isInSubState;
                property bool selectionChangedFlag: selectionChanged
                property variant backendValue: backendValues.selectionMode
                property variant valueFromBackend: backendValue.value
                property string enumScope: "AbstractGraph3D"
                property string enumSeparator: " | "
                property int checkedCount: 0
                property bool selectionItem: false
                property bool selectionRow: false
                property bool selectionColumn: false
                property bool selectionSlice: false
                property bool selectionMulti: false

                function checkValue(checkedVariable, variableText, expressionBase) {
                    var expressionStr = expressionBase
                    if (checkedVariable) {
                        if (expressionStr !== "") {
                            expressionStr += enumSeparator
                        }
                        expressionStr += enumScope
                        expressionStr += "."
                        expressionStr += variableText
                        checkedCount++
                    }
                    return expressionStr
                }

                function composeSelectionMode() {
                    var expressionStr = ""
                    checkedCount = 0
                    expressionStr = checkValue(selectionItem, "SelectionItem", expressionStr)
                    expressionStr = checkValue(selectionRow, "SelectionRow", expressionStr)
                    expressionStr = checkValue(selectionColumn, "SelectionColumn", expressionStr)
                    expressionStr = checkValue(selectionSlice, "SelectionSlice", expressionStr)
                    expressionStr = checkValue(selectionMulti, "SelectionMultiSeries", expressionStr)

                    if (checkedCount === 0)
                        backendValue.expression = enumScope + ".SelectionNone"
                    else
                        backendValue.expression = expressionStr
                }

                function evaluate() {
                    if (backendValue.value === undefined)
                        return

                    selectionItem = (backendValue.expression.indexOf("SelectionItem") !== -1)
                    selectionRow = (backendValue.expression.indexOf("SelectionRow") !== -1)
                    selectionColumn = (backendValue.expression.indexOf("SelectionColumn") !== -1)
                    selectionSlice = (backendValue.expression.indexOf("SelectionSlice") !== -1)
                    selectionMulti = (backendValue.expression.indexOf("SelectionMultiSeries") !== -1)

                    selectionItemBox.checked = selectionItem
                    selectionRowBox.checked = selectionRow
                    selectionColumnBox.checked = selectionColumn
                    selectionSliceBox.checked = selectionSlice
                    selectionMultiSeriesBox.checked = selectionMulti
                }

                onSelectionChangedFlagChanged: evaluate()

                onIsInModelChanged: evaluate()

                onIsInSubStateChanged: evaluate()

                onBackendValueChanged: evaluate()

                onValueFromBackendChanged: evaluate()

                ColumnLayout {
                    anchors.fill: parent

                    Controls.CheckBox {
                        id: selectionItemBox
                        style: checkBox.style
                        text: "SelectionItem"
                        Layout.fillWidth: true
                        onClicked: {
                            selectionLayout.selectionItem = checked
                            selectionLayout.composeSelectionMode()
                        }
                    }
                    Controls.CheckBox {
                        id: selectionRowBox
                        style: checkBox.style
                        text: "SelectionRow"
                        Layout.fillWidth: true
                        onClicked: {
                            selectionLayout.selectionRow = checked
                            selectionLayout.composeSelectionMode()
                        }
                    }
                    Controls.CheckBox {
                        id: selectionColumnBox
                        style: checkBox.style
                        text: "SelectionColumn"
                        Layout.fillWidth: true
                        onClicked: {
                            selectionLayout.selectionColumn = checked
                            selectionLayout.composeSelectionMode()
                        }
                    }
                    Controls.CheckBox {
                        id: selectionSliceBox
                        style: checkBox.style
                        text: "SelectionSlice"
                        Layout.fillWidth: true
                        onClicked: {
                            selectionLayout.selectionSlice = checked
                            selectionLayout.composeSelectionMode()
                        }
                    }
                    Controls.CheckBox {
                        id: selectionMultiSeriesBox
                        style: checkBox.style
                        text: "SelectionMultiSeries"
                        Layout.fillWidth: true
                        onClicked: {
                            selectionLayout.selectionMulti = checked
                            selectionLayout.composeSelectionMode()
                        }
                    }
                }
            }
            Label {
                text: qsTr("measureFps")
                tooltip: qsTr("Measure Frames Per Second")
                Layout.fillWidth: true
            }
            SecondColumnLayout {
                CheckBox {
                    backendValue: backendValues.measureFps
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("orthoProjection")
                tooltip: qsTr("Use Orthographic Projection")
                Layout.fillWidth: true
            }
            SecondColumnLayout {
                CheckBox {
                    backendValue: backendValues.orthoProjection
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("aspectRatio")
                tooltip: qsTr("Aspect Ratio")
                Layout.fillWidth: true
            }
            SecondColumnLayout {
                SpinBox {
                    backendValue: backendValues.aspectRatio
                    minimumValue: 0.01
                    maximumValue: 100.0
                    stepSize: 0.01
                    decimals: 2
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("floorLevel")
                tooltip: qsTr("Floor Level")
                Layout.fillWidth: true
            }
            SecondColumnLayout {
                LineEdit {
                    backendValue: backendValues.floorLevel
                    inputMethodHints: Qt.ImhFormattedNumbersOnly
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("horizontalAspectRatio")
                tooltip: qsTr("Horizontal Aspect Ratio")
                Layout.fillWidth: true
            }
            SecondColumnLayout {
                SpinBox {
                    backendValue: backendValues.horizontalAspectRatio
                    minimumValue: 0.0
                    maximumValue: 100.0
                    stepSize: 0.01
                    decimals: 2
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("reflection")
                tooltip: qsTr("Reflection")
                Layout.fillWidth: true
            }
            SecondColumnLayout {
                CheckBox {
                    id: reflectionCheckbox
                    backendValue: backendValues.reflection
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("reflectivity")
                tooltip: qsTr("Reflectivity")
                Layout.fillWidth: true
                visible: reflectionCheckbox.checked
            }
            SecondColumnLayout {
                visible: reflectionCheckbox.checked
                SpinBox {
                    backendValue: backendValues.reflectivity
                    minimumValue: 0.0
                    maximumValue: 1.0
                    stepSize: 0.01
                    decimals: 1
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("margin")
                tooltip: qsTr("Graph Margin")
                Layout.fillWidth: true
            }
            SecondColumnLayout {
                SpinBox {
                    backendValue: backendValues.margin
                    minimumValue: -1.0
                    maximumValue: 100.0
                    stepSize: 0.1
                    decimals: 1
                    Layout.fillWidth: true
                }
            }

            // Kept for debugging
            Label { }
            SecondColumnLayout {
                TextEdit {
                    id: debugLabel
                    Layout.fillWidth: true
                    wrapMode: TextEdit.WordWrap
                    textFormat: TextEdit.RichText
                    width: 400
                    visible: false
                }
            }
            Controls.CheckBox {
                property color textColor: colorLogic.textColor
                id: checkBox
                style: CustomCheckBoxStyle {}
                visible: false
                ColorLogic {
                    id: colorLogic
                    backendValue: backendValues.selectionMode
                }
            }
        }
    }
}
