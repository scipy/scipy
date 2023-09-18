/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtScxml module of the Qt Toolkit.
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

#ifndef QSCXMLTABLEDATA_P_H
#define QSCXMLTABLEDATA_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <QtScxml/qscxmltabledata.h>
#include <QtCore/qsharedpointer.h>
#include <QtCore/qhash.h>
#include <QtCore/qvector.h>

#include <functional>

QT_BEGIN_NAMESPACE
class QTextStream;
class QScxmlInvokableServiceFactory;

namespace DocumentModel {
struct ScxmlDocument;
}

namespace QScxmlInternal {
class Q_SCXML_EXPORT GeneratedTableData: public QScxmlTableData
{
public:
    typedef std::function<
        int(const QScxmlExecutableContent::InvokeInfo &invokeInfo,
            const QVector<QScxmlExecutableContent::StringId> &namelist,
            const QVector<QScxmlExecutableContent::ParameterInfo> &params,
            QSharedPointer<DocumentModel::ScxmlDocument> content)
    > CreateFactoryId;

    struct MetaDataInfo {
        QStringList stateNames;
    };

    struct DataModelInfo {
        QHash<QScxmlExecutableContent::EvaluatorId, QString> stringEvaluators;
        QHash<QScxmlExecutableContent::EvaluatorId, QString> boolEvaluators;
        QHash<QScxmlExecutableContent::EvaluatorId, QString> variantEvaluators;
        QHash<QScxmlExecutableContent::EvaluatorId, QString> voidEvaluators;
    };

public:
    static void build(DocumentModel::ScxmlDocument *doc, GeneratedTableData *table,
                      MetaDataInfo *metaDataInfo, DataModelInfo *dataModelInfo,
                      CreateFactoryId func);
    static QString toString(const int *stateMachineTable);

public:
    QString string(QScxmlExecutableContent::StringId id) const override final;
    QScxmlExecutableContent::InstructionId *instructions() const override final;
    QScxmlExecutableContent::EvaluatorInfo evaluatorInfo(
            QScxmlExecutableContent::EvaluatorId evaluatorId) const override final;
    QScxmlExecutableContent::AssignmentInfo assignmentInfo(
            QScxmlExecutableContent::EvaluatorId assignmentId) const override final;
    QScxmlExecutableContent::ForeachInfo foreachInfo(
            QScxmlExecutableContent::EvaluatorId foreachId) const override final;
    QScxmlExecutableContent::StringId *dataNames(int *count) const override final;
    QScxmlExecutableContent::ContainerId initialSetup() const override final;
    QString name() const override final;
    const qint32 *stateMachineTable() const override final;
    QScxmlInvokableServiceFactory *serviceFactory(int id) const override;

public:
    QVector<qint32> theStateMachineTable;
    QStringList theStrings;
    QVector<qint32> theInstructions;
    QVector<QScxmlExecutableContent::EvaluatorInfo> theEvaluators;
    QVector<QScxmlExecutableContent::AssignmentInfo> theAssignments;
    QVector<QScxmlExecutableContent::ForeachInfo> theForeaches;
    QVector<QScxmlExecutableContent::StringId> theDataNameIds;
    QScxmlExecutableContent::ContainerId theInitialSetup;
    int theName;
};
} // QScxmlInternal namespace

QT_END_NAMESPACE

#endif // QSCXMLTABLEDATA_P_H
