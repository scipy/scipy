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

#ifndef QSCXMLEXECUTABLECONTENT_H
#define QSCXMLEXECUTABLECONTENT_H

#include <QtScxml/qscxmlglobals.h>

QT_BEGIN_NAMESPACE

namespace QScxmlExecutableContent {

typedef qint32 ContainerId;
enum { NoContainer = -1 };
typedef qint32 StringId;
enum { NoString = -1 };
typedef qint32 InstructionId;
enum { NoInstruction = -1 };
typedef qint32 EvaluatorId;
enum { NoEvaluator = -1 };

#if defined(Q_CC_MSVC) || defined(Q_CC_GNU)
#pragma pack(push, 4) // 4 == sizeof(qint32)
#endif
struct EvaluatorInfo {
    StringId expr;
    StringId context;
};

struct AssignmentInfo {
    StringId dest;
    StringId expr;
    StringId context;
};

struct ForeachInfo {
    StringId array;
    StringId item;
    StringId index;
    StringId context;
};

struct ParameterInfo {
    StringId name;
    EvaluatorId expr;
    StringId location;
};

struct InvokeInfo {
    StringId id;
    StringId prefix;
    StringId location;
    StringId context;
    EvaluatorId expr;
    ContainerId finalize;
    bool autoforward;
};
#if defined(Q_CC_MSVC) || defined(Q_CC_GNU)
#pragma pack(pop)
#endif

} // QScxmlExecutableContent namespace

QT_END_NAMESPACE

#endif // QSCXMLEXECUTABLECONTENT_H
