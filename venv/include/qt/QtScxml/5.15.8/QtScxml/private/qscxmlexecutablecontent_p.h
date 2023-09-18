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

#ifndef QSCXMLEXECUTABLECONTENT_P_H
#define QSCXMLEXECUTABLECONTENT_P_H

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

#include <QtScxml/qscxmlexecutablecontent.h>
#include <QtScxml/private/qscxmltabledata_p.h>
#include <QtScxml/private/qscxmlcompiler_p.h>
#include <QtCore/qtextstream.h>

#ifndef BUILD_QSCXMLC
#include <QtScxml/qscxmldatamodel.h>
#include <QtScxml/qscxmlstatemachine.h>
#endif // BUILD_QSCXMLC

QT_BEGIN_NAMESPACE

namespace QScxmlExecutableContent {

static inline bool operator<(const EvaluatorInfo &ei1, const EvaluatorInfo &ei2)
{
    if (ei1.expr != ei2.expr)
        return ei1.expr < ei2.expr;
    else
        return ei1.context < ei2.context;
}

static inline bool operator<(const AssignmentInfo &ai1, const AssignmentInfo &ai2)
{
    if (ai1.dest != ai2.dest)
        return ai1.dest < ai2.dest;
    else if (ai1.expr != ai2.expr)
        return ai1.expr < ai2.expr;
    else
        return ai1.context < ai2.context;
}

static inline bool operator<(const ForeachInfo &fi1, const ForeachInfo &fi2)
{
    if (fi1.array != fi2.array) return fi1.array < fi2.array;
    if (fi1.item != fi2.item) return fi1.item < fi2.item;
    if (fi1.index != fi2.index) return fi1.index < fi2.index;
    return fi1.context < fi2.context;
}

#if defined(Q_CC_MSVC) || defined(Q_CC_GNU)
#pragma pack(push, 4) // 4 == sizeof(qint32)
#endif

template <typename T>
struct Array
{
    qint32 count;
    // T[] data;
    T *data() { return const_cast<T *>(const_data()); }
    const T *const_data() const { return reinterpret_cast<const T *>(reinterpret_cast<const char *>(this) + sizeof(Array<T>)); }

    const T &at(int pos) const { return *(const_data() + pos); }
    int dataSize() const { return count * sizeof(T) / sizeof(qint32); }
    int size() const { return sizeof(Array<T>) / sizeof(qint32) + dataSize(); }
};

struct Q_SCXML_EXPORT Instruction
{
    enum InstructionType: qint32 {
        Sequence = 1,
        Sequences,
        Send,
        Raise,
        Log,
        JavaScript,
        Assign,
        Initialize,
        If,
        Foreach,
        Cancel,
        DoneData
    } instructionType;
};

struct Q_SCXML_EXPORT DoneData: Instruction
{
    StringId location;
    StringId contents;
    EvaluatorId expr;
    Array<ParameterInfo> params;

    static InstructionType kind() { return Instruction::DoneData; }
};

struct Q_SCXML_EXPORT InstructionSequence: Instruction
{
    qint32 entryCount; // the amount of qint32's that the instructions take up
    // Instruction[] instructions;

    static InstructionType kind() { return Instruction::Sequence; }
    const InstructionId *instructions() const
    {
        return reinterpret_cast<const InstructionId *>(this)
                + sizeof(InstructionSequence) / sizeof(qint32);
    }
    int size() const { return sizeof(InstructionSequence) / sizeof(qint32) + entryCount; }
};

struct Q_SCXML_EXPORT InstructionSequences: Instruction
{
    qint32 sequenceCount;
    qint32 entryCount; // the amount of qint32's that the sequences take up
    // InstructionSequence[] sequences;

    static InstructionType kind() { return Instruction::Sequences; }
    const InstructionSequence *sequences() const {
        return reinterpret_cast<const InstructionSequence *>(
                    reinterpret_cast<const InstructionId *>(this)
                    + sizeof(InstructionSequences) / sizeof(qint32));
    }
    int size() const { return sizeof(InstructionSequences)/sizeof(qint32) + entryCount; }
    const InstructionId *at(int pos) const
    {
        const InstructionId *seq = reinterpret_cast<const InstructionId *>(sequences());
        while (pos--) {
            seq += reinterpret_cast<const InstructionSequence *>(seq)->size();
        }
        return seq;
    }
};

struct Q_SCXML_EXPORT Send: Instruction
{
    StringId instructionLocation;
    StringId event;
    EvaluatorId eventexpr;
    StringId type;
    EvaluatorId typeexpr;
    StringId target;
    EvaluatorId targetexpr;
    StringId id;
    StringId idLocation;
    StringId delay;
    EvaluatorId delayexpr;
    StringId content;
    EvaluatorId contentexpr;
    Array<StringId> namelist;
//    Array<Param> params;

    static InstructionType kind() { return Instruction::Send; }

    int paramsOffset() const
    {
        return sizeof(Send) / sizeof(qint32) + namelist.dataSize();
    }

    int size() const
    {
        return paramsOffset() + params()->size();
    }

    const Array<ParameterInfo> *params() const {
        return reinterpret_cast<const Array<ParameterInfo> *>(
                    reinterpret_cast<const InstructionId *>(this) + paramsOffset());
    }

    Array<ParameterInfo> *params() {
        return reinterpret_cast<Array<ParameterInfo> *>(
                    reinterpret_cast<InstructionId *>(this) + paramsOffset());
    }

    static int calculateExtraSize(int paramCount, int nameCount) {
        return 1 + paramCount * sizeof(ParameterInfo) / sizeof(qint32)
                + nameCount * sizeof(StringId) / sizeof(qint32);
    }
};

struct Q_SCXML_EXPORT Raise: Instruction
{
    StringId event;

    static InstructionType kind() { return Instruction::Raise; }
    int size() const { return sizeof(Raise) / sizeof(qint32); }
};

struct Q_SCXML_EXPORT Log: Instruction
{
    StringId label;
    EvaluatorId expr;

    static InstructionType kind() { return Instruction::Log; }
    int size() const { return sizeof(Log) / sizeof(qint32); }
};

struct Q_SCXML_EXPORT JavaScript: Instruction
{
    EvaluatorId go;

    static InstructionType kind() { return Instruction::JavaScript; }
    int size() const { return sizeof(JavaScript) / sizeof(qint32); }
};

struct Q_SCXML_EXPORT Assign: Instruction
{
    EvaluatorId expression;

    static InstructionType kind() { return Instruction::Assign; }
    int size() const { return sizeof(Assign) / sizeof(qint32); }
};

struct Q_SCXML_EXPORT Initialize: Instruction
{
    EvaluatorId expression;

    static InstructionType kind() { return Instruction::Initialize; }
    int size() const { return sizeof(Initialize) / sizeof(qint32); }
};

struct Q_SCXML_EXPORT If: Instruction
{
    Array<EvaluatorId> conditions;
    // InstructionSequences blocks;
    const InstructionSequences *blocks() const {
        return reinterpret_cast<const InstructionSequences *>(
                    reinterpret_cast<const InstructionId *>(this) + sizeof(If) / sizeof(qint32)
                    + conditions.dataSize());
    }

    static InstructionType kind() { return Instruction::If; }
    int size() const
    {
        return sizeof(If) / sizeof(qint32) + blocks()->size() + conditions.dataSize();
    }
};

struct Q_SCXML_EXPORT Foreach: Instruction
{
    EvaluatorId doIt;
    InstructionSequence block;

    static InstructionType kind() { return Instruction::Foreach; }
    int size() const { return sizeof(Foreach) / sizeof(qint32) + block.entryCount; }
    const InstructionId *blockstart() const
    {
        return reinterpret_cast<const InstructionId *>(&block);
    }
};

struct Q_SCXML_EXPORT Cancel: Instruction
{
    StringId sendid;
    EvaluatorId sendidexpr;

    static InstructionType kind() { return Instruction::Cancel; }
    int size() const { return sizeof(Cancel) / sizeof(qint32); }
};

struct StateTable {
    int version;
    int name;
    enum: int {
        InvalidDataModel = -1,
        NullDataModel = 0,
        EcmaScriptDataModel = 1,
        CppDataModel = 2
    } dataModel;
    int childStates; // offset into offsets
    int initialTransition;
    int initialSetup;
    enum: int { InvalidBinding = -1, EarlyBinding = 0, LateBinding = 1 } binding;
    int maxServiceId;
    int stateOffset, stateCount;
    int transitionOffset, transitionCount;
    int arrayOffset, arraySize;

    enum { terminator = 0xc0ff33 };
    enum { InvalidIndex = -1 };

    struct State {
        int name;
        int parent;
        enum: int {
            Invalid = -1,
            Normal = 0,
            Parallel = 1,
            Final = 2,
            ShallowHistory = 3,
            DeepHistory = 4
        } type;
        int initialTransition;
        int initInstructions;
        int entryInstructions;
        int exitInstructions;
        int doneData;
        int childStates; // offset into arrays
        int transitions; // offset into arrays
        int serviceFactoryIds; // offset into arrays

        State()
            : name(InvalidIndex)
            , parent(InvalidIndex)
            , type(Invalid)
            , initialTransition(InvalidIndex)
            , initInstructions(InvalidIndex)
            , entryInstructions(InvalidIndex)
            , exitInstructions(InvalidIndex)
            , doneData(InvalidIndex)
            , childStates(InvalidIndex)
            , transitions(InvalidIndex)
            , serviceFactoryIds(InvalidIndex)
        {}

        bool isAtomic() const
        { return childStates == InvalidIndex; }

        bool isCompound() const
        { return type == Normal && childStates != InvalidIndex; }

        bool parentIsScxmlElement() const
        { return parent == InvalidIndex; }

        bool isHistoryState() const
        { return type == ShallowHistory || type == DeepHistory; }

        bool isParallel() const
        { return type == Parallel; }
    };

    struct Transition {
        int events; // offset into offsets
        int condition;
        enum: int {
            Invalid = -1,
            Internal = 0,
            External = 1,
            Synthetic = 2
        } type;
        int source;
        int targets; // offset into offsets
        int transitionInstructions;

        Transition()
            : events(InvalidIndex)
            , condition(InvalidIndex)
            , type(Invalid)
            , source(InvalidIndex)
            , targets(InvalidIndex)
            , transitionInstructions(InvalidIndex)
        {}
    };

    struct Array {
        Array(const int *start): start(start) {}
        int size() const { return *start; }
        bool isValid() const { return start != nullptr; }

        int operator[](int idx) const {
            Q_ASSERT(idx >= 0);
            Q_ASSERT(idx < size());
            return *(start + idx + 1);
        }

        struct const_iterator: public std::iterator<std::forward_iterator_tag, int, ptrdiff_t,
                                                    const int *, const int &>
        {
            const_iterator(const Array &a, int pos): a(a), pos(pos) {}

            const_iterator &operator++() {
                if (pos < a.size()) ++pos;
                return *this;
            }

            bool operator==(const const_iterator &other) const
            { return &other.a == &a && other.pos == pos; }

            bool operator!=(const StateTable::Array::const_iterator &other)
            { return !this->operator==(other); }

            int operator*() const {
                if (pos < a.size())
                    return a[pos];
                else
                    return -1;
            }

        private:
            const Array &a;
            int pos;
        };

        const_iterator begin() const
        { return const_iterator(*this, 0); }

        const_iterator end() const
        { return const_iterator(*this, size()); }

    private:
        const int *start;
    };

    StateTable()
        : version(InvalidIndex)
        , name(InvalidIndex)
        , dataModel(InvalidDataModel)
        , childStates(InvalidIndex)
        , initialTransition(InvalidIndex)
        , initialSetup(InvalidIndex)
        , binding(InvalidBinding)
        , maxServiceId(InvalidIndex)
        , stateOffset(InvalidIndex), stateCount(InvalidIndex)
        , transitionOffset(InvalidIndex), transitionCount(InvalidIndex)
        , arrayOffset(InvalidIndex), arraySize(InvalidIndex)
    {}

    const State &state(int idx) const
    {
        Q_ASSERT(idx >= 0);
        Q_ASSERT(idx < stateCount);
        return reinterpret_cast<const State *>(
                    reinterpret_cast<const int *>(this) + stateOffset)[idx];
    }

    const Transition &transition(int idx) const
    {
        Q_ASSERT(idx >= 0);
        Q_ASSERT(idx < transitionCount);
        return reinterpret_cast<const Transition *>(
                    reinterpret_cast<const int *>(this) + transitionOffset)[idx];
    }

    const Array array(int idx) const
    {
        Q_ASSERT(idx < arraySize);
        if (idx >= 0) {
            const int *start = reinterpret_cast<const int *>(this) + arrayOffset + idx;
            Q_ASSERT(*start + idx < arraySize);
            return Array(start);
        } else {
            return Array(nullptr);
        }
    }
};

#if defined(Q_CC_MSVC) || defined(Q_CC_GNU)
#pragma pack(pop)
#endif

} // QScxmlExecutableContent namespace

class QScxmlExecutionEngine
{
    Q_DISABLE_COPY(QScxmlExecutionEngine)

public:
    QScxmlExecutionEngine(QScxmlStateMachine *stateMachine);

    bool execute(QScxmlExecutableContent::ContainerId ip, const QVariant &extraData = QVariant());

private:
    const QScxmlExecutableContent::InstructionId *step(
            const QScxmlExecutableContent::InstructionId *ip, bool *ok);

    QScxmlStateMachine *stateMachine;
    QVariant extraData;
};

QT_END_NAMESPACE

#endif // QSCXMLEXECUTABLECONTENT_P_H
