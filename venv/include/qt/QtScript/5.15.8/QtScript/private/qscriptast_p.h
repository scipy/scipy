/****************************************************************************
**
** Copyright (C) 2015 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the QtScript module of the Qt Toolkit.
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

#ifndef QSCRIPTAST_P_H
#define QSCRIPTAST_P_H

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

#include <QtCore/qstring.h>

#include "qscriptastvisitor_p.h"

QT_BEGIN_NAMESPACE

class QScriptNameIdImpl;

namespace QSOperator // ### rename
{

enum Op {
    Add,
    And,
    InplaceAnd,
    Assign,
    BitAnd,
    BitOr,
    BitXor,
    InplaceSub,
    Div,
    InplaceDiv,
    Equal,
    Ge,
    Gt,
    In,
    InplaceAdd,
    InstanceOf,
    Le,
    LShift,
    InplaceLeftShift,
    Lt,
    Mod,
    InplaceMod,
    Mul,
    InplaceMul,
    NotEqual,
    Or,
    InplaceOr,
    RShift,
    InplaceRightShift,
    StrictEqual,
    StrictNotEqual,
    Sub,
    URShift,
    InplaceURightShift,
    InplaceXor
};

} // namespace QSOperator

namespace QScript { namespace AST {

class Node
{
public:
    enum Kind {
        Kind_Node,
        Kind_ExpressionNode,
        Kind_Statement,
        Kind_ThisExpression,
        Kind_IdentifierExpression,
        Kind_NullExpression,
        Kind_TrueLiteral,
        Kind_FalseLiteral,
        Kind_NumericLiteral,
        Kind_StringLiteral,
        Kind_RegExpLiteral,
        Kind_ArrayLiteral,
        Kind_ObjectLiteral,
        Kind_ElementList,
        Kind_Elision,
        Kind_PropertyNameAndValueList,
        Kind_PropertyName,
        Kind_IdentifierPropertyName,
        Kind_StringLiteralPropertyName,
        Kind_NumericLiteralPropertyName,
        Kind_ArrayMemberExpression,
        Kind_FieldMemberExpression,
        Kind_NewMemberExpression,
        Kind_NewExpression,
        Kind_CallExpression,
        Kind_ArgumentList,
        Kind_PostIncrementExpression,
        Kind_PostDecrementExpression,
        Kind_DeleteExpression,
        Kind_VoidExpression,
        Kind_TypeOfExpression,
        Kind_PreIncrementExpression,
        Kind_PreDecrementExpression,
        Kind_UnaryPlusExpression,
        Kind_UnaryMinusExpression,
        Kind_TildeExpression,
        Kind_NotExpression,
        Kind_BinaryExpression,
        Kind_ConditionalExpression,
        Kind_Expression,
        Kind_Block,
        Kind_StatementList,
        Kind_VariableStatement,
        Kind_VariableDeclarationList,
        Kind_VariableDeclaration,
        Kind_EmptyStatement,
        Kind_ExpressionStatement,
        Kind_IfStatement,
        Kind_DoWhileStatement,
        Kind_WhileStatement,
        Kind_ForStatement,
        Kind_LocalForStatement,
        Kind_ForEachStatement,
        Kind_LocalForEachStatement,
        Kind_ContinueStatement,
        Kind_BreakStatement,
        Kind_ReturnStatement,
        Kind_WithStatement,
        Kind_SwitchStatement,
        Kind_CaseBlock,
        Kind_CaseClauses,
        Kind_CaseClause,
        Kind_DefaultClause,
        Kind_LabelledStatement,
        Kind_ThrowStatement,
        Kind_TryStatement,
        Kind_Catch,
        Kind_Finally,
        Kind_FunctionDeclaration,
        Kind_FunctionExpression,
        Kind_FormalParameterList,
        Kind_FunctionBody,
        Kind_Program,
        Kind_SourceElements,
        Kind_SourceElement,
        Kind_FunctionSourceElement,
        Kind_StatementSourceElement,
        Kind_DebuggerStatement
    };

    inline Node():
        startLine(0), startColumn(0),
        endLine(0), endColumn(0), kind(Kind_Node) {}

    virtual ~Node() {}

    virtual ExpressionNode *expressionCast();
    virtual BinaryExpression *binaryExpressionCast();
    virtual Statement *statementCast();

    inline void accept(Visitor *visitor)
    {
        if (visitor->preVisit(this)) {
            accept0(visitor);
            visitor->postVisit(this);
        }
    }

    static void acceptChild(Node *node, Visitor *visitor)
    {
        if (node)
            node->accept(visitor);
    }

    virtual void accept0(Visitor *visitor) = 0;

    int startLine;
    int startColumn;
    int endLine;
    int endColumn;
    Kind kind;
};

class ExpressionNode: public Node
{
public:
    ExpressionNode() { kind = Kind_ExpressionNode; }
    virtual ~ExpressionNode() {}

    virtual ExpressionNode *expressionCast();
};

class Statement: public Node
{
public:
    Statement() { kind = Kind_Statement; }
    virtual ~Statement() {}

    virtual Statement *statementCast();
};

class ThisExpression: public ExpressionNode
{
public:
    ThisExpression() { kind = Kind_ThisExpression; }
    virtual ~ThisExpression() {}

    virtual void accept0(Visitor *visitor);
};

class IdentifierExpression: public ExpressionNode
{
public:
    IdentifierExpression(QScriptNameIdImpl *n):
        name (n) { kind = Kind_IdentifierExpression; }

    virtual ~IdentifierExpression() {}

    virtual void accept0(Visitor *visitor);

// attributes
    QScriptNameIdImpl *name;
};

class NullExpression: public ExpressionNode
{
public:
    NullExpression() { kind = Kind_NullExpression; }
    virtual ~NullExpression() {}

    virtual void accept0(Visitor *visitor);
};

class TrueLiteral: public ExpressionNode
{
public:
    TrueLiteral() { kind = Kind_TrueLiteral; }
    virtual ~TrueLiteral() {}

    virtual void accept0(Visitor *visitor);
};

class FalseLiteral: public ExpressionNode
{
public:
    FalseLiteral() { kind = Kind_FalseLiteral; }
    virtual ~FalseLiteral() {}

    virtual void accept0(Visitor *visitor);
};

class NumericLiteral: public ExpressionNode
{
public:
    NumericLiteral(double v):
        value (v) { kind = Kind_NumericLiteral; }
    virtual ~NumericLiteral() {}

    virtual void accept0(Visitor *visitor);

// attributes:
    double value;
};

class StringLiteral: public ExpressionNode
{
public:
    StringLiteral(QScriptNameIdImpl *v):
        value (v) { kind = Kind_StringLiteral; }

    virtual ~StringLiteral() {}

    virtual void accept0(Visitor *visitor);

// attributes:
    QScriptNameIdImpl *value;
};

class RegExpLiteral: public ExpressionNode
{
public:
    RegExpLiteral(QScriptNameIdImpl *p, int f):
        pattern (p), flags (f) { kind = Kind_RegExpLiteral; }

    virtual ~RegExpLiteral() {}

    virtual void accept0(Visitor *visitor);

// attributes:
    QScriptNameIdImpl *pattern;
    int flags;
};

class ArrayLiteral: public ExpressionNode
{
public:
    ArrayLiteral(Elision *e):
        elements (0), elision (e)
        { kind = Kind_ArrayLiteral; }

    ArrayLiteral(ElementList *elts):
        elements (elts), elision (0)
        { kind = Kind_ArrayLiteral; }

    ArrayLiteral(ElementList *elts, Elision *e):
        elements (elts), elision (e)
        { kind = Kind_ArrayLiteral; }

    virtual ~ArrayLiteral() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ElementList *elements;
    Elision *elision;
};

class ObjectLiteral: public ExpressionNode
{
public:
    ObjectLiteral():
        properties (0) { kind = Kind_ObjectLiteral; }

    ObjectLiteral(PropertyNameAndValueList *plist):
        properties (plist) { kind = Kind_ObjectLiteral; }

    virtual ~ObjectLiteral() {}

    virtual void accept0(Visitor *visitor);

// attributes
    PropertyNameAndValueList *properties;
};

class ElementList: public Node
{
public:
    ElementList(Elision *e, ExpressionNode *expr):
        elision (e), expression (expr), next (this)
        { kind = Kind_ElementList; }

    ElementList(ElementList *previous, Elision *e, ExpressionNode *expr):
        elision (e), expression (expr)
    {
        kind = Kind_ElementList;
        next = previous->next;
        previous->next = this;
    }

    virtual ~ElementList() {}

    inline ElementList *finish ()
    {
        ElementList *front = next;
        next = 0;
        return front;
    }

    virtual void accept0(Visitor *visitor);

// attributes
    Elision *elision;
    ExpressionNode *expression;
    ElementList *next;
};

class Elision: public Node
{
public:
    Elision():
        next (this) { kind = Kind_Elision; }

    Elision(Elision *previous)
    {
        kind = Kind_Elision;
        next = previous->next;
        previous->next = this;
    }

    virtual ~Elision() {}

    virtual void accept0(Visitor *visitor);

    inline Elision *finish ()
    {
        Elision *front = next;
        next = 0;
        return front;
    }

// attributes
    Elision *next;
};

class PropertyNameAndValueList: public Node
{
public:
    PropertyNameAndValueList(PropertyName *n, ExpressionNode *v):
        name (n), value (v), next (this)
        { kind = Kind_PropertyNameAndValueList; }

    PropertyNameAndValueList(PropertyNameAndValueList *previous, PropertyName *n, ExpressionNode *v):
        name (n), value (v)
    {
        kind = Kind_PropertyNameAndValueList;
        next = previous->next;
        previous->next = this;
    }

    virtual ~PropertyNameAndValueList() {}

    virtual void accept0(Visitor *visitor);

    inline PropertyNameAndValueList *finish ()
    {
        PropertyNameAndValueList *front = next;
        next = 0;
        return front;
    }

// attributes
    PropertyName *name;
    ExpressionNode *value;
    PropertyNameAndValueList *next;
};

class PropertyName: public Node
{
public:
    PropertyName() { kind = Kind_PropertyName; }
    virtual ~PropertyName() {}
};

class IdentifierPropertyName: public PropertyName
{
public:
    IdentifierPropertyName(QScriptNameIdImpl *n):
        id (n) { kind = Kind_IdentifierPropertyName; }

    virtual ~IdentifierPropertyName() {}

    virtual void accept0(Visitor *visitor);

// attributes
    QScriptNameIdImpl *id;
};

class StringLiteralPropertyName: public PropertyName
{
public:
    StringLiteralPropertyName(QScriptNameIdImpl *n):
        id (n) { kind = Kind_StringLiteralPropertyName; }
    virtual ~StringLiteralPropertyName() {}

    virtual void accept0(Visitor *visitor);

// attributes
    QScriptNameIdImpl *id;
};

class NumericLiteralPropertyName: public PropertyName
{
public:
    NumericLiteralPropertyName(double n):
        id (n) { kind = Kind_NumericLiteralPropertyName; }
    virtual ~NumericLiteralPropertyName() {}

    virtual void accept0(Visitor *visitor);

// attributes
    double id;
};

class ArrayMemberExpression: public ExpressionNode
{
public:
    ArrayMemberExpression(ExpressionNode *b, ExpressionNode *e):
        base (b), expression (e)
        { kind = Kind_ArrayMemberExpression; }

    virtual ~ArrayMemberExpression() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *base;
    ExpressionNode *expression;
};

class FieldMemberExpression: public ExpressionNode
{
public:
    FieldMemberExpression(ExpressionNode *b, QScriptNameIdImpl *n):
        base (b), name (n)
        { kind = Kind_FieldMemberExpression; }

    virtual ~FieldMemberExpression() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *base;
    QScriptNameIdImpl *name;
};

class NewMemberExpression: public ExpressionNode
{
public:
    NewMemberExpression(ExpressionNode *b, ArgumentList *a):
        base (b), arguments (a)
        { kind = Kind_NewMemberExpression; }

    virtual ~NewMemberExpression() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *base;
    ArgumentList *arguments;
};

class NewExpression: public ExpressionNode
{
public:
    NewExpression(ExpressionNode *e):
        expression (e) { kind = Kind_NewExpression; }

    virtual ~NewExpression() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *expression;
};

class CallExpression: public ExpressionNode
{
public:
    CallExpression(ExpressionNode *b, ArgumentList *a):
        base (b), arguments (a)
        { kind = Kind_CallExpression; }

    virtual ~CallExpression() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *base;
    ArgumentList *arguments;
};

class ArgumentList: public Node
{
public:
    ArgumentList(ExpressionNode *e):
        expression (e), next (this)
        { kind = Kind_ArgumentList; }

    ArgumentList(ArgumentList *previous, ExpressionNode *e):
        expression (e)
    {
        kind = Kind_ArgumentList;
        next = previous->next;
        previous->next = this;
    }

    virtual ~ArgumentList() {}

    virtual void accept0(Visitor *visitor);

    inline ArgumentList *finish ()
    {
        ArgumentList *front = next;
        next = 0;
        return front;
    }

// attributes
    ExpressionNode *expression;
    ArgumentList *next;
};

class PostIncrementExpression: public ExpressionNode
{
public:
    PostIncrementExpression(ExpressionNode *b):
        base (b) { kind = Kind_PostIncrementExpression; }

    virtual ~PostIncrementExpression() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *base;
};

class PostDecrementExpression: public ExpressionNode
{
public:
    PostDecrementExpression(ExpressionNode *b):
        base (b) { kind = Kind_PostDecrementExpression; }

    virtual ~PostDecrementExpression() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *base;
};

class DeleteExpression: public ExpressionNode
{
public:
    DeleteExpression(ExpressionNode *e):
        expression (e) { kind = Kind_DeleteExpression; }
    virtual ~DeleteExpression() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *expression;
};

class VoidExpression: public ExpressionNode
{
public:
    VoidExpression(ExpressionNode *e):
        expression (e) { kind = Kind_VoidExpression; }

    virtual ~VoidExpression() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *expression;
};

class TypeOfExpression: public ExpressionNode
{
public:
    TypeOfExpression(ExpressionNode *e):
        expression (e) { kind = Kind_TypeOfExpression; }

    virtual ~TypeOfExpression() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *expression;
};

class PreIncrementExpression: public ExpressionNode
{
public:
    PreIncrementExpression(ExpressionNode *e):
        expression (e) { kind = Kind_PreIncrementExpression; }

    virtual ~PreIncrementExpression() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *expression;
};

class PreDecrementExpression: public ExpressionNode
{
public:
    PreDecrementExpression(ExpressionNode *e):
        expression (e) { kind = Kind_PreDecrementExpression; }

    virtual ~PreDecrementExpression() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *expression;
};

class UnaryPlusExpression: public ExpressionNode
{
public:
    UnaryPlusExpression(ExpressionNode *e):
        expression (e) { kind = Kind_UnaryPlusExpression; }

    virtual ~UnaryPlusExpression() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *expression;
};

class UnaryMinusExpression: public ExpressionNode
{
public:
    UnaryMinusExpression(ExpressionNode *e):
        expression (e) { kind = Kind_UnaryMinusExpression; }

    virtual ~UnaryMinusExpression() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *expression;
};

class TildeExpression: public ExpressionNode
{
public:
    TildeExpression(ExpressionNode *e):
        expression (e) { kind = Kind_TildeExpression; }

    virtual ~TildeExpression() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *expression;
};

class NotExpression: public ExpressionNode
{
public:
    NotExpression(ExpressionNode *e):
        expression (e) { kind = Kind_NotExpression; }

    virtual ~NotExpression() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *expression;
};

class BinaryExpression: public ExpressionNode
{
public:
    BinaryExpression(ExpressionNode *l, int o, ExpressionNode *r):
        left (l), op (o), right (r)
        { kind = Kind_BinaryExpression; }

    virtual ~BinaryExpression() {}

    virtual BinaryExpression *binaryExpressionCast();

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *left;
    int op;
    ExpressionNode *right;
};

class ConditionalExpression: public ExpressionNode
{
public:
    ConditionalExpression(ExpressionNode *e, ExpressionNode *t, ExpressionNode *f):
        expression (e), ok (t), ko (f)
        { kind = Kind_ConditionalExpression; }

    virtual ~ConditionalExpression() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *expression;
    ExpressionNode *ok;
    ExpressionNode *ko;
};

class Expression: public ExpressionNode // ### rename
{
public:
    Expression(ExpressionNode *l, ExpressionNode *r):
        left (l), right (r) { kind = Kind_Expression; }

    virtual ~Expression() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *left;
    ExpressionNode *right;
};

class Block: public Statement
{
public:
    Block(StatementList *slist):
        statements (slist) { kind = Kind_Block; }

    virtual ~Block() {}

    virtual void accept0(Visitor *visitor);

// attributes
    StatementList *statements;
};

class StatementList: public Node
{
public:
    StatementList(Statement *stmt):
        statement (stmt), next (this)
        { kind = Kind_StatementList; }

    StatementList(StatementList *previous, Statement *stmt):
        statement (stmt)
    {
        kind = Kind_StatementList;
        next = previous->next;
        previous->next = this;
    }

    virtual ~StatementList() {}

    virtual void accept0(Visitor *visitor);

    inline StatementList *finish ()
    {
        StatementList *front = next;
        next = 0;
        return front;
    }

// attributes
    Statement *statement;
    StatementList *next;
};

class VariableStatement: public Statement
{
public:
    VariableStatement(VariableDeclarationList *vlist):
        declarations (vlist)
        { kind = Kind_VariableStatement; }

    virtual ~VariableStatement() {}

    virtual void accept0(Visitor *visitor);

// attributes
    VariableDeclarationList *declarations;
};

class VariableDeclaration: public Node
{
public:
    VariableDeclaration(QScriptNameIdImpl *n, ExpressionNode *e):
        name (n), expression (e), readOnly(false)
        { kind = Kind_VariableDeclaration; }

    virtual ~VariableDeclaration() {}

    virtual void accept0(Visitor *visitor);

// attributes
    QScriptNameIdImpl *name;
    ExpressionNode *expression;
    bool readOnly;
};

class VariableDeclarationList: public Node
{
public:
    VariableDeclarationList(VariableDeclaration *decl):
        declaration (decl), next (this)
        { kind = Kind_VariableDeclarationList; }

    VariableDeclarationList(VariableDeclarationList *previous, VariableDeclaration *decl):
        declaration (decl)
    {
        kind = Kind_VariableDeclarationList;
        next = previous->next;
        previous->next = this;
    }

    virtual ~VariableDeclarationList() {}

    virtual void accept0(Visitor *visitor);

    inline VariableDeclarationList *finish (bool readOnly)
    {
        VariableDeclarationList *front = next;
        next = 0;
        if (readOnly) {
            VariableDeclarationList *vdl;
            for (vdl = front; vdl != 0; vdl = vdl->next)
                vdl->declaration->readOnly = true;
        }
        return front;
    }

// attributes
    VariableDeclaration *declaration;
    VariableDeclarationList *next;
};

class EmptyStatement: public Statement
{
public:
    EmptyStatement() { kind = Kind_EmptyStatement; }
    virtual ~EmptyStatement() {}

    virtual void accept0(Visitor *visitor);
};

class ExpressionStatement: public Statement
{
public:
    ExpressionStatement(ExpressionNode *e):
        expression (e) { kind = Kind_ExpressionStatement; }

    virtual ~ExpressionStatement() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *expression;
};

class IfStatement: public Statement
{
public:
    IfStatement(ExpressionNode *e, Statement *t, Statement *f = 0):
        expression (e), ok (t), ko (f)
        { kind = Kind_IfStatement; }

    virtual ~IfStatement() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *expression;
    Statement *ok;
    Statement *ko;
};

class DoWhileStatement: public Statement
{
public:
    DoWhileStatement(Statement *stmt, ExpressionNode *e):
        statement (stmt), expression (e)
        { kind = Kind_DoWhileStatement; }

    virtual ~DoWhileStatement() {}

    virtual void accept0(Visitor *visitor);

// attributes
    Statement *statement;
    ExpressionNode *expression;
};

class WhileStatement: public Statement
{
public:
    WhileStatement(ExpressionNode *e, Statement *stmt):
        expression (e), statement (stmt)
        { kind = Kind_WhileStatement; }

    virtual ~WhileStatement() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *expression;
    Statement *statement;
};

class ForStatement: public Statement
{
public:
    ForStatement(ExpressionNode *i, ExpressionNode *c, ExpressionNode *e, Statement *stmt):
        initialiser (i), condition (c), expression (e), statement (stmt)
        { kind = Kind_ForStatement; }

    virtual ~ForStatement() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *initialiser;
    ExpressionNode *condition;
    ExpressionNode *expression;
    Statement *statement;
};

class LocalForStatement: public Statement
{
public:
    LocalForStatement(VariableDeclarationList *vlist, ExpressionNode *c, ExpressionNode *e, Statement *stmt):
        declarations (vlist), condition (c), expression (e), statement (stmt)
        { kind = Kind_LocalForStatement; }

    virtual ~LocalForStatement() {}

    virtual void accept0(Visitor *visitor);

// attributes
    VariableDeclarationList *declarations;
    ExpressionNode *condition;
    ExpressionNode *expression;
    Statement *statement;
};

class ForEachStatement: public Statement
{
public:
    ForEachStatement(ExpressionNode *i, ExpressionNode *e, Statement *stmt):
        initialiser (i), expression (e), statement (stmt)
        { kind = Kind_ForEachStatement; }

    virtual ~ForEachStatement() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *initialiser;
    ExpressionNode *expression;
    Statement *statement;
};

class LocalForEachStatement: public Statement
{
public:
    LocalForEachStatement(VariableDeclaration *v, ExpressionNode *e, Statement *stmt):
        declaration (v), expression (e), statement (stmt)
        { kind = Kind_LocalForEachStatement; }

    virtual ~LocalForEachStatement() {}

    virtual void accept0(Visitor *visitor);

// attributes
    VariableDeclaration *declaration;
    ExpressionNode *expression;
    Statement *statement;
};

class ContinueStatement: public Statement
{
public:
    ContinueStatement(QScriptNameIdImpl *l = 0):
        label (l) { kind = Kind_ContinueStatement; }

    virtual ~ContinueStatement() {}

    virtual void accept0(Visitor *visitor);

// attributes
    QScriptNameIdImpl *label;
};

class BreakStatement: public Statement
{
public:
    BreakStatement(QScriptNameIdImpl *l = 0):
        label (l) { kind = Kind_BreakStatement; }

    virtual ~BreakStatement() {}

    virtual void accept0(Visitor *visitor);

// attributes
    QScriptNameIdImpl *label;
};

class ReturnStatement: public Statement
{
public:
    ReturnStatement(ExpressionNode *e):
        expression (e) { kind = Kind_ReturnStatement; }

    virtual ~ReturnStatement() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *expression;
};

class WithStatement: public Statement
{
public:
    WithStatement(ExpressionNode *e, Statement *stmt):
        expression (e), statement (stmt)
        { kind = Kind_WithStatement; }

    virtual ~WithStatement() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *expression;
    Statement *statement;
};

class SwitchStatement: public Statement
{
public:
    SwitchStatement(ExpressionNode *e, CaseBlock *b):
        expression (e), block (b)
        { kind = Kind_SwitchStatement; }

    virtual ~SwitchStatement() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *expression;
    CaseBlock *block;
};

class CaseBlock: public Node
{
public:
    CaseBlock(CaseClauses *c, DefaultClause *d = 0, CaseClauses *r = 0):
        clauses (c), defaultClause (d), moreClauses (r)
        { kind = Kind_CaseBlock; }

    virtual ~CaseBlock() {}

    virtual void accept0(Visitor *visitor);

// attributes
    CaseClauses *clauses;
    DefaultClause *defaultClause;
    CaseClauses *moreClauses;
};

class CaseClauses: public Node
{
public:
    CaseClauses(CaseClause *c):
        clause (c), next (this)
        { kind = Kind_CaseClauses; }

    CaseClauses(CaseClauses *previous, CaseClause *c):
        clause (c)
    {
        kind = Kind_CaseClauses;
        next = previous->next;
        previous->next = this;
    }

    virtual ~CaseClauses() {}

    virtual void accept0(Visitor *visitor);

    inline CaseClauses *finish ()
    {
        CaseClauses *front = next;
        next = 0;
        return front;
    }

//attributes
    CaseClause *clause;
    CaseClauses *next;
};

class CaseClause: public Node
{
public:
    CaseClause(ExpressionNode *e, StatementList *slist):
        expression (e), statements (slist)
        { kind = Kind_CaseClause; }

    virtual ~CaseClause() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *expression;
    StatementList *statements;
};

class DefaultClause: public Node
{
public:
    DefaultClause(StatementList *slist):
        statements (slist)
        { kind = Kind_DefaultClause; }

    virtual ~DefaultClause() {}

    virtual void accept0(Visitor *visitor);

// attributes
    StatementList *statements;
};

class LabelledStatement: public Statement
{
public:
    LabelledStatement(QScriptNameIdImpl *l, Statement *stmt):
        label (l), statement (stmt)
        { kind = Kind_LabelledStatement; }

    virtual ~LabelledStatement() {}

    virtual void accept0(Visitor *visitor);

// attributes
    QScriptNameIdImpl *label;
    Statement *statement;
};

class ThrowStatement: public Statement
{
public:
    ThrowStatement(ExpressionNode *e):
        expression (e) { kind = Kind_ThrowStatement; }

    virtual ~ThrowStatement() {}

    virtual void accept0(Visitor *visitor);

// attributes
    ExpressionNode *expression;
};

class TryStatement: public Statement
{
public:
    TryStatement(Statement *stmt, Catch *c, Finally *f):
        statement (stmt), catchExpression (c), finallyExpression (f)
        { kind = Kind_TryStatement; }

    TryStatement(Statement *stmt, Finally *f):
        statement (stmt), catchExpression (0), finallyExpression (f)
        { kind = Kind_TryStatement; }

    TryStatement(Statement *stmt, Catch *c):
        statement (stmt), catchExpression (c), finallyExpression (0)
        { kind = Kind_TryStatement; }

    virtual ~TryStatement() {}

    virtual void accept0(Visitor *visitor);

// attributes
    Statement *statement;
    Catch *catchExpression;
    Finally *finallyExpression;
};

class Catch: public Node
{
public:
    Catch(QScriptNameIdImpl *n, Statement *stmt):
        name (n), statement (stmt)
        { kind = Kind_Catch; }

    virtual ~Catch() {}

    virtual void accept0(Visitor *visitor);

// attributes
    QScriptNameIdImpl *name;
    Statement *statement;
};

class Finally: public Node
{
public:
    Finally(Statement *stmt):
        statement (stmt)
        { kind = Kind_Finally; }

    virtual ~Finally() {}

    virtual void accept0(Visitor *visitor);

// attributes
    Statement *statement;
};

class FunctionExpression: public ExpressionNode
{
public:
    FunctionExpression(QScriptNameIdImpl *n, FormalParameterList *f, FunctionBody *b):
        name (n), formals (f), body (b)
        { kind = Kind_FunctionExpression; }

    virtual ~FunctionExpression() {}

    virtual void accept0(Visitor *visitor);

// attributes
    QScriptNameIdImpl *name;
    FormalParameterList *formals;
    FunctionBody *body;
};

class FunctionDeclaration: public FunctionExpression
{
public:
    FunctionDeclaration(QScriptNameIdImpl *n, FormalParameterList *f, FunctionBody *b):
        FunctionExpression(n, f, b)
        { kind = Kind_FunctionDeclaration; }

    virtual ~FunctionDeclaration() {}

    virtual void accept0(Visitor *visitor);
};

class FormalParameterList: public Node
{
public:
    FormalParameterList(QScriptNameIdImpl *n):
        name (n), next (this)
        { kind = Kind_FormalParameterList; }

    FormalParameterList(FormalParameterList *previous, QScriptNameIdImpl *n):
        name (n)
    {
        kind = Kind_FormalParameterList;
        next = previous->next;
        previous->next = this;
    }

    virtual ~FormalParameterList() {}

    virtual void accept0(Visitor *visitor);

    inline FormalParameterList *finish ()
    {
        FormalParameterList *front = next;
        next = 0;
        return front;
    }

// attributes
    QScriptNameIdImpl *name;
    FormalParameterList *next;
};

class FunctionBody: public Node
{
public:
    FunctionBody(SourceElements *elts):
        elements (elts)
        { kind = Kind_FunctionBody; }

    virtual ~FunctionBody() {}

    virtual void accept0(Visitor *visitor);

// attributes
    SourceElements *elements;
};

class Program: public Node
{
public:
    Program(SourceElements *elts):
        elements (elts)
        { kind = Kind_Program; }

    virtual ~Program() {}

    virtual void accept0(Visitor *visitor);

// attributes
    SourceElements *elements;
};

class SourceElements: public Node
{
public:
    SourceElements(SourceElement *elt):
        element (elt), next (this)
        { kind = Kind_SourceElements; }

    SourceElements(SourceElements *previous, SourceElement *elt):
        element (elt)
    {
        kind = Kind_SourceElements;
        next = previous->next;
        previous->next = this;
    }

    virtual ~SourceElements() {}

    virtual void accept0(Visitor *visitor);

    inline SourceElements *finish ()
    {
        SourceElements *front = next;
        next = 0;
        return front;
    }

// attributes
    SourceElement *element;
    SourceElements *next;
};

class SourceElement: public Node
{
public:
    inline SourceElement()
        { kind = Kind_SourceElement; }

    virtual ~SourceElement() {}
};

class FunctionSourceElement: public SourceElement
{
public:
    FunctionSourceElement(FunctionDeclaration *f):
        declaration (f)
        { kind = Kind_FunctionSourceElement; }

    virtual ~FunctionSourceElement() {}

    virtual void accept0(Visitor *visitor);

// attributes
    FunctionDeclaration *declaration;
};

class StatementSourceElement: public SourceElement
{
public:
    StatementSourceElement(Statement *stmt):
        statement (stmt)
        { kind = Kind_StatementSourceElement; }

    virtual ~StatementSourceElement() {}

    virtual void accept0(Visitor *visitor);

// attributes
    Statement *statement;
};

class DebuggerStatement: public Statement
{
public:
    DebuggerStatement()
        { kind = Kind_DebuggerStatement; }

    virtual ~DebuggerStatement() {}

    virtual void accept0(Visitor *visitor);
};

} } // namespace AST

QT_END_NAMESPACE

#endif
