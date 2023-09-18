----------------------------------------------------------------------------
--
-- Copyright (C) 2014-2015 Ford Motor Company.
-- Contact: https://www.qt.io/licensing/
--
-- This file is part of the QtRemoteObjects module of the Qt Toolkit.
--
-- $QT_BEGIN_LICENSE:LGPL$
-- Commercial License Usage
-- Licensees holding valid commercial Qt licenses may use this file in
-- accordance with the commercial license agreement provided with the
-- Software or, alternatively, in accordance with the terms contained in
-- a written agreement between you and The Qt Company. For licensing terms
-- and conditions see https://www.qt.io/terms-conditions. For further
-- information use the contact form at https://www.qt.io/contact-us.
--
-- GNU Lesser General Public License Usage
-- Alternatively, this file may be used under the terms of the GNU Lesser
-- General Public License version 3 as published by the Free Software
-- Foundation and appearing in the file LICENSE.LGPL3 included in the
-- packaging of this file. Please review the following information to
-- ensure the GNU Lesser General Public License version 3 requirements
-- will be met: https://www.gnu.org/licenses/lgpl-3.0.html.
--
-- GNU General Public License Usage
-- Alternatively, this file may be used under the terms of the GNU
-- General Public License version 2.0 or (at your option) the GNU General
-- Public license version 3 or any later version approved by the KDE Free
-- Qt Foundation. The licenses are as published by the Free Software
-- Foundation and appearing in the file LICENSE.GPL2 and LICENSE.GPL3
-- included in the packaging of this file. Please review the following
-- information to ensure the GNU General Public License requirements will
-- be met: https://www.gnu.org/licenses/gpl-2.0.html and
-- https://www.gnu.org/licenses/gpl-3.0.html.
--
-- $QT_END_LICENSE$
--
----------------------------------------------------------------------------

%parser rep_grammar
%decl repparser.h
%impl repparser.cpp

%token_prefix Token_
%token semicolon "[semicolon];"
%token class "[class]class[ \\t]+(?<name>[A-Za-z_][A-Za-z0-9_]+)[ \\t]*"
%token pod "[pod]POD[ \\t]*(?<name>[A-Za-z_][A-Za-z0-9_]+)[ \\t]*\\((?<types>[^\\)]*)\\);?[ \\t]*"
%token enum "[enum][ \\t]*ENUM[ \\t]+(?<name>[A-Za-z_][A-Za-z0-9_]*)[ \\t]*"
%token enum_param "[enum_param][ \\t]*(?<name>[A-Za-z_][A-Za-z0-9_]*)[ \\t]*(=[ \\t]*(?<value>-\\d+|0[xX][0-9A-Fa-f]+|\\d+))?[ \\t]*"
%token prop "[prop][ \\t]*PROP[ \\t]*\\((?<args>[^\\)]+)\\);?[ \\t]*"
%token use_enum "[use_enum]USE_ENUM[ \\t]*\\((?<name>[^\\)]*)\\);?[ \\t]*"
%token signal "[signal][ \\t]*SIGNAL[ \\t]*\\([ \\t]*(?<name>\\S+)[ \\t]*\\((?<args>[^\\)]*)\\)[ \\t]*\\);?[ \\t]*"
%token slot "[slot][ \\t]*SLOT[ \\t]*\\((?<type>[^\\(]*)\\((?<args>[^\\)]*)\\)[ \\t]*\\);?[ \\t]*"
%token model "[model][ \\t]*MODEL[ \\t]+(?<name>[A-Za-z_][A-Za-z0-9_]+)\\((?<args>[^\\)]+)\\)[ \\t]*;?[ \\t]*"
%token childrep "[childrep][ \\t]*CLASS[ \\t]+(?<name>[A-Za-z_][A-Za-z0-9_]+)\\((?<type>[^\\)]+)\\)[ \\t]*;?[ \\t]*"
%token start "[start][ \\t]*\\{[ \\t]*"
%token stop "[stop][ \\t]*\\};?[ \\t]*"
%token comma "[comma],"
%token comment "[comment](?<comment>[ \\t]*//[^\\n]*\\n)"
%token preprocessor_directive "[preprocessor_directive](?<preprocessor_directive>#[ \\t]*[^\\n]*\\n)"
%token newline "[newline](\\r)?\\n"

%start TopLevel

/:
#ifndef REPPARSER_H
#define REPPARSER_H

#include <rep_grammar_p.h>
#include <qregexparser.h>
#include <QStringList>
#include <QVector>
#include <QRegExp>

QT_BEGIN_NAMESPACE
class QIODevice;

/// A property of a Class declaration
struct ASTProperty
{
    enum Modifier
    {
        Constant,
        ReadOnly,
        ReadPush,
        ReadWrite,
        SourceOnlySetter
    };

    ASTProperty();
    ASTProperty(const QString &type, const QString &name, const QString &defaultValue, Modifier modifier, bool persisted,
                bool isPointer=false);

    QString type;
    QString name;
    QString defaultValue;
    Modifier modifier;
    bool persisted;
    bool isPointer;
};
Q_DECLARE_TYPEINFO(ASTProperty, Q_MOVABLE_TYPE);

struct ASTDeclaration
{
    enum VariableType {
        None = 0,
        Constant = 1,
        Reference = 2,
    };
    Q_DECLARE_FLAGS(VariableTypes, VariableType)

    ASTDeclaration(const QString &declarationType = QString(), const QString &declarationName = QString(), VariableTypes declarationVariableType = None)
        : type(declarationType),
          name(declarationName),
          variableType(declarationVariableType)
    {
    }

    QString asString(bool withName) const;

    QString type;
    QString name;
    VariableTypes variableType;
};
Q_DECLARE_TYPEINFO(ASTDeclaration, Q_MOVABLE_TYPE);

struct ASTFunction
{
    enum ParamsAsStringFormat {
        Default,
        Normalized
    };

    explicit ASTFunction(const QString &name = QString(), const QString &returnType = QLatin1String("void"));

    QString paramsAsString(ParamsAsStringFormat format = Default) const;
    QStringList paramNames() const;

    QString returnType;
    QString name;
    QVector<ASTDeclaration> params;
};
Q_DECLARE_TYPEINFO(ASTFunction, Q_MOVABLE_TYPE);

struct ASTEnumParam
{
    ASTEnumParam(const QString &paramName = QString(), int paramValue = 0)
        : name(paramName),
          value(paramValue)
    {
    }

    QString asString() const;

    QString name;
    int value;
};
Q_DECLARE_TYPEINFO(ASTEnumParam, Q_MOVABLE_TYPE);

struct ASTEnum
{
    explicit ASTEnum(const QString &name = QString());

    QString name;
    QVector<ASTEnumParam> params;
    bool isSigned;
    int max;
};
Q_DECLARE_TYPEINFO(ASTEnum, Q_MOVABLE_TYPE);

struct ASTModelRole
{
    ASTModelRole(const QString &roleName = QString())
        : name(roleName)
    {
    }

    QString name;
};
Q_DECLARE_TYPEINFO(ASTModelRole, Q_MOVABLE_TYPE);

struct ASTModel
{
    ASTModel(int index = -1) : propertyIndex(index) {}

    QVector<ASTModelRole> roles;
    int propertyIndex;
};
Q_DECLARE_TYPEINFO(ASTModel, Q_MOVABLE_TYPE);

/// A Class declaration
struct ASTClass
{
    explicit ASTClass(const QString& name = QString());

    bool isValid() const;
    bool hasPointerObjects() const;

    QString name;
    QVector<ASTProperty> properties;
    QVector<ASTFunction> signalsList;
    QVector<ASTFunction> slotsList;
    QVector<ASTEnum> enums;
    bool hasPersisted;
    QVector<ASTModel> modelMetadata;
    QVector<int> subClassPropertyIndices;
};
Q_DECLARE_TYPEINFO(ASTClass, Q_MOVABLE_TYPE);

// The attribute of a POD
struct PODAttribute
{
    explicit PODAttribute(const QString &type_ = QString(), const QString &name_ = QString())
        : type(type_),
          name(name_)
    {}
    QString type;
    QString name;
};
Q_DECLARE_TYPEINFO(PODAttribute, Q_MOVABLE_TYPE);

// A POD declaration
struct POD
{
    QString name;
    QVector<PODAttribute> attributes;
};
Q_DECLARE_TYPEINFO(POD, Q_MOVABLE_TYPE);

// The AST representation of a .rep file
struct AST
{
    QVector<ASTClass> classes;
    QVector<POD> pods;
    QVector<ASTEnum> enums;
    QVector<QString> enumUses;
    QStringList preprocessorDirectives;
};
Q_DECLARE_TYPEINFO(AST, Q_MOVABLE_TYPE);

class RepParser: public QRegexParser<RepParser, $table>
{
public:
    explicit RepParser(QIODevice &outputDevice);
    virtual ~RepParser() {}

    bool parse() override { return QRegexParser<RepParser, $table>::parse(); }

    void reset() override;
    int nextToken();
    bool consumeRule(int ruleno);

    AST ast() const;

private:
    struct TypeParser
    {
        void parseArguments(const QString &arguments);
        void appendParams(ASTFunction &slot);
        void appendPods(POD &pods);
        void generateFunctionParameter(QString variableName, const QString &propertyType, int &variableNameIndex, ASTDeclaration::VariableTypes variableType);
        //Type, Variable
        QList<ASTDeclaration> arguments;
    };

    bool parseProperty(ASTClass &astClass, const QString &propertyDeclaration);
    /// A helper function to parse modifier flag of property declaration
    bool parseModifierFlag(const QString &flag, ASTProperty::Modifier &modifier, bool &persisted);

    bool parseRoles(ASTModel &astModel, const QString &modelRoles);

    AST m_ast;

    ASTClass m_astClass;
    ASTEnum m_astEnum;
    int m_astEnumValue;
};
QT_END_NAMESPACE
#endif
:/


/.
#include "repparser.h"

#include <QDebug>
#include <QTextStream>

// for normalizeTypeInternal
#include <private/qmetaobject_p.h>
#include <private/qmetaobject_moc_p.h>

// Code copied from moc.cpp
// We cannot depend on QMetaObject::normalizedSignature,
// since repc is linked against Qt5Bootstrap (which doesn't offer QMetaObject) when cross-compiling
// Thus, just use internal API which is exported in private headers, as moc does
static QByteArray normalizeType(const QByteArray &ba, bool fixScope = false)
{
    const char *s = ba.constData();
    int len = ba.size();
    char stackbuf[64];
    char *buf = (len >= 64 ? new char[len + 1] : stackbuf);
    char *d = buf;
    char last = 0;
    while (*s && is_space(*s))
        s++;
    while (*s) {
        while (*s && !is_space(*s))
            last = *d++ = *s++;
        while (*s && is_space(*s))
            s++;
        if (*s && ((is_ident_char(*s) && is_ident_char(last))
                   || ((*s == ':') && (last == '<')))) {
            last = *d++ = ' ';
        }
    }
    *d = '\0';
    QByteArray result = normalizeTypeInternal(buf, d, fixScope);
    if (buf != stackbuf)
        delete [] buf;
    return result;
}

ASTProperty::ASTProperty()
    : modifier(ReadPush), persisted(false), isPointer(false)
{
}

ASTProperty::ASTProperty(const QString &type, const QString &name, const QString &defaultValue, Modifier modifier, bool persisted, bool isPointer)
    : type(type), name(name), defaultValue(defaultValue), modifier(modifier), persisted(persisted), isPointer(isPointer)
{
}

QString ASTDeclaration::asString(bool withName) const
{
    QString str;
    if (variableType & ASTDeclaration::Constant)
        str += QLatin1String("const ");
    str += type;
    if (variableType & ASTDeclaration::Reference)
        str += QLatin1String(" &");
    if (withName)
        str += QString::fromLatin1(" %1").arg(name);
    return str;
}

ASTFunction::ASTFunction(const QString &name, const QString &returnType)
    : returnType(returnType), name(name)
{
}

QString ASTFunction::paramsAsString(ParamsAsStringFormat format) const
{
    QString str;
    for (const ASTDeclaration &param : params) {
        QString paramStr = param.asString(format != Normalized);
        if (format == Normalized) {
            paramStr = QString::fromLatin1(::normalizeType(paramStr.toLatin1().constData()));
            str += paramStr + QLatin1Char(',');
        } else {
            str += paramStr + QLatin1String(", ");
        }
    }

    str.chop((format == Normalized ? 1 : 2)); // chop trailing ',' or ', '

    return str;
}

QStringList ASTFunction::paramNames() const
{
    QStringList names;
    names.reserve(params.size());
    for (const ASTDeclaration &param : params)
        names << param.name;
    return names;
}

ASTEnum::ASTEnum(const QString &name)
    : name(name), isSigned(false), max(0)
{
}

ASTClass::ASTClass(const QString &name)
    : name(name), hasPersisted(false)
{
}

bool ASTClass::isValid() const
{
    return !name.isEmpty();
}

bool ASTClass::hasPointerObjects() const
{
    int count = modelMetadata.size() + subClassPropertyIndices.size();
    return count > 0;
}

RepParser::RepParser(QIODevice &outputDevice)
    : QRegexParser(), m_astEnumValue(-1)
{
    setBufferFromDevice(&outputDevice);
}

void RepParser::reset()
{
    m_ast = AST();
    m_astClass = ASTClass();
    m_astEnum = ASTEnum();
    //setDebug();
}

bool RepParser::parseModifierFlag(const QString &flag, ASTProperty::Modifier &modifier, bool &persisted)
{
    QRegExp regex(QStringLiteral("\\s*,\\s*"));
    QStringList flags = flag.split(regex);
    persisted = flags.removeAll(QStringLiteral("PERSISTED")) > 0;
    if (flags.length() == 0)
        return true;
    if (flags.length() > 1) {
        // Only valid combination is "READONLY" and "CONSTANT"
        if (flags.length() == 2 && flags.contains(QStringLiteral("READONLY")) &&
            flags.contains(QStringLiteral("CONSTANT"))) {
            // If we have READONLY and CONSTANT that means CONSTANT
            modifier = ASTProperty::Constant;
            return true;
        } else {
            setErrorString(QLatin1String("Invalid property declaration: combination not allowed (%1)").arg(flag));
            return false;
        }
    }
    const QString &f = flags.at(0);
    if (f == QLatin1String("READONLY"))
        modifier = ASTProperty::ReadOnly;
    else if (f == QLatin1String("CONSTANT"))
        modifier = ASTProperty::Constant;
    else if (f == QLatin1String("READPUSH"))
        modifier = ASTProperty::ReadPush;
    else if (f == QLatin1String("READWRITE"))
        modifier = ASTProperty::ReadWrite;
    else if (f == QLatin1String("SOURCEONLYSETTER"))
        modifier = ASTProperty::SourceOnlySetter;
    else {
        setErrorString(QLatin1String("Invalid property declaration: flag %1 is unknown").arg(flag));
        return false;
    }

    return true;
}

QString stripArgs(const QString &arguments)
{
    // This repc parser searches for the longest possible matches, which can be multiline.
    // This method "cleans" the string input, removing comments and converting to a single
    // line for subsequent parsing.
    QStringList lines = arguments.split(QRegExp(QStringLiteral("\r?\n")));
    for (auto & line : lines)
        line.replace(QRegExp(QStringLiteral("//.*")),QString());
    return lines.join(QString());
}

bool RepParser::parseProperty(ASTClass &astClass, const QString &propertyDeclaration)
{
    QString input = stripArgs(propertyDeclaration).trimmed();
    const QRegExp whitespace(QStringLiteral("\\s"));

    QString propertyType;
    QString propertyName;
    QString propertyDefaultValue;
    ASTProperty::Modifier propertyModifier = ASTProperty::ReadPush;
    bool persisted = false;

    // parse type declaration which could be a nested template as well
    bool inTemplate = false;
    int templateDepth = 0;
    int nameIndex = -1;

    for (int i = 0; i < input.size(); ++i) {
        const QChar inputChar(input.at(i));
        if (inputChar == QLatin1Char('<')) {
            propertyType += inputChar;
            inTemplate = true;
            ++templateDepth;
        } else if (inputChar == QLatin1Char('>')) {
            propertyType += inputChar;
            --templateDepth;
            if (templateDepth == 0)
                inTemplate = false;
        } else if (inputChar.isSpace()) {
            if (!inTemplate) {
                nameIndex = i;
                break;
            } else {
                propertyType += inputChar;
            }
        } else {
            propertyType += inputChar;
        }
    }

    if (nameIndex == -1) {
        setErrorString(QLatin1String("PROP: Invalid property declaration: %1").arg(propertyDeclaration));
        return false;
    }

    // parse the name of the property
    input = input.mid(nameIndex).trimmed();

    const int equalSignIndex = input.indexOf(QLatin1Char('='));
    if (equalSignIndex != -1) { // we have a default value
        propertyName = input.left(equalSignIndex).trimmed();

        input = input.mid(equalSignIndex + 1).trimmed();
        const int lastQuoteIndex = input.lastIndexOf(QLatin1Char('"'));
        if (lastQuoteIndex != -1) {
            propertyDefaultValue = input.left(lastQuoteIndex + 1);
            input = input.mid(lastQuoteIndex + 1);
        }
        const int whitespaceIndex = input.indexOf(whitespace);
        if (whitespaceIndex == -1) { // no flag given
            if (propertyDefaultValue.isEmpty())
                propertyDefaultValue = input;
            propertyModifier = ASTProperty::ReadPush;
        } else { // flag given
            if (propertyDefaultValue.isEmpty())
                propertyDefaultValue = input.left(whitespaceIndex).trimmed();

            const QString flag = input.mid(whitespaceIndex + 1).trimmed();
            if (!parseModifierFlag(flag, propertyModifier, persisted))
                return false;
        }
    } else { // there is no default value
        const int whitespaceIndex = input.indexOf(whitespace);
        if (whitespaceIndex == -1) { // no flag given
            propertyName = input;
            propertyModifier = ASTProperty::ReadPush;
        } else { // flag given
            propertyName = input.left(whitespaceIndex).trimmed();

            const QString flag = input.mid(whitespaceIndex + 1).trimmed();
            if (!parseModifierFlag(flag, propertyModifier, persisted))
                return false;
        }
    }

    astClass.properties << ASTProperty(propertyType, propertyName, propertyDefaultValue, propertyModifier, persisted);
    if (persisted)
        astClass.hasPersisted = true;
    return true;
}

bool RepParser::parseRoles(ASTModel &astModel, const QString &modelRoles)
{
    const QString input = modelRoles.trimmed();

    if (input.isEmpty())
        return true;

    const QStringList roleStrings = input.split(QChar(QLatin1Char(',')));
    for (auto role : roleStrings)
        astModel.roles << ASTModelRole(role.trimmed());
    return true;
}

AST RepParser::ast() const
{
    return m_ast;
}

void RepParser::TypeParser::parseArguments(const QString &arguments)
{
    const QString strippedArgs = stripArgs(arguments);
    int templateDepth = 0;
    bool inTemplate = false;
    bool inVariable = false;
    QString propertyType;
    QString variableName;
    ASTDeclaration::VariableTypes variableType = ASTDeclaration::None;
    int variableNameIndex = 0;
    for (int i = 0; i < strippedArgs.size(); ++i) {
        const QChar inputChar(strippedArgs.at(i));
        if (inputChar == QLatin1Char('<')) {
            propertyType += inputChar;
            inTemplate = true;
            ++templateDepth;
        } else if (inputChar == QLatin1Char('>')) {
            propertyType += inputChar;
            --templateDepth;
            if (templateDepth == 0)
                inTemplate = false;
        } else if (inputChar.isSpace()) {
            if (inTemplate)
                propertyType += inputChar;
            else if (!propertyType.isEmpty()) {
                if (propertyType == QLatin1String("const")) {
                    propertyType.clear();
                    variableType |= ASTDeclaration::Constant;
                } else {
                    inVariable = true;
                }
            }
        } else if (inputChar == QLatin1Char('&')) {
            variableType |= ASTDeclaration::Reference;
        } else if (inputChar == QLatin1Char(',')) {
            if (!inTemplate) {
                RepParser::TypeParser::generateFunctionParameter(variableName, propertyType, variableNameIndex, variableType);
                propertyType.clear();
                variableName.clear();
                variableType = ASTDeclaration::None;
                inVariable = false;
            } else {
                propertyType += inputChar;
            }
        } else {
            if (inVariable)
                variableName += inputChar;
            else
                propertyType += inputChar;
        }
    }
    if (!propertyType.isEmpty()) {
        RepParser::TypeParser::generateFunctionParameter(variableName, propertyType, variableNameIndex, variableType);
    }
}

void RepParser::TypeParser::generateFunctionParameter(QString variableName, const QString &propertyType, int &variableNameIndex, ASTDeclaration::VariableTypes variableType)
{
    if (!variableName.isEmpty())
        variableName = variableName.trimmed();
    else
        variableName = QString::fromLatin1("__repc_variable_%1").arg(++variableNameIndex);
    arguments.append(ASTDeclaration(propertyType, variableName, variableType));
}

void RepParser::TypeParser::appendParams(ASTFunction &slot)
{
    for (const ASTDeclaration &arg : qAsConst(arguments))
        slot.params << arg;
}

void RepParser::TypeParser::appendPods(POD &pods)
{
    for (const ASTDeclaration &arg : qAsConst(arguments)) {
        PODAttribute attr;
        attr.type = arg.type;
        attr.name = arg.name;
        pods.attributes.append(qMove(attr));
    }
}

bool RepParser::consumeRule(int ruleno)
{
    if (isDebug()) {
        qDebug() << "consumeRule:" << ruleno << spell[rule_info[rule_index[ruleno]]];
    }
    switch (ruleno) {
./

TopLevel: Types | Newlines Types | FileComments Types | Newlines FileComments Types;

FileComments: Comments;

Types: Type | Type Types;

Newlines: newline | newline Newlines;
Comments: Comment | Comment Comments;
Comment: comment | comment Newlines;
Type: PreprocessorDirective | PreprocessorDirective Newlines;
Type: Pod | Pod Newlines;
Type: Class;
Type: UseEnum | UseEnum Newlines;
Type: Comments | Comments Newlines;
Type: Enum;
/.
    case $rule_number:
    {
        m_ast.enums.append(m_astEnum);
    }
    break;
./

Comma: comma | comma Newlines;

PreprocessorDirective: preprocessor_directive;
/.
    case $rule_number:
    {
        m_ast.preprocessorDirectives.append(captured().value(QStringLiteral("preprocessor_directive")));
    }
    break;
./

Pod: pod;
/.
    case $rule_number:
    {
        POD pod;
        pod.name = captured().value(QStringLiteral("name")).trimmed();

        const QString argString = captured().value(QLatin1String("types")).trimmed();
        if (argString.isEmpty()) {
            qWarning() << "[repc] - Ignoring POD with no data members.  POD name: " << qPrintable(pod.name);
            return true;
        }

        RepParser::TypeParser parseType;
        parseType.parseArguments(argString);
        parseType.appendPods(pod);
        m_ast.pods.append(pod);
    }
    break;
./

Class: ClassStart Start ClassTypes Stop;
/.
    case $rule_number:
./
Class: ClassStart Start Comments Stop;
/.
    case $rule_number:
./
Class: ClassStart Start Stop;
/.
    case $rule_number:
    {
        m_ast.classes.append(m_astClass);
    }
    break;
./

ClassTypes: ClassType | ClassType ClassTypes;
ClassType: DecoratedProp | DecoratedSignal | DecoratedSlot | DecoratedModel | DecoratedClass | Comments;
ClassType: Enum;
/.
    case $rule_number:
    {
        m_astClass.enums.append(m_astEnum);
    }
    break;
./

DecoratedSlot: Slot | Comments Slot | Slot Newlines | Comments Slot Newlines;
DecoratedSignal: Signal | Comments Signal | Signal Newlines | Comments Signal Newlines;
DecoratedProp: Prop | Comments Prop | Prop Newlines | Comments Prop Newlines;
DecoratedModel: Model | Comments Model | Model Newlines | Comments Model Newlines;
DecoratedClass: ChildRep | Comments ChildRep | ChildRep Newlines | Comments ChildRep Newlines;
DecoratedEnumParam: EnumParam | Comments EnumParam | EnumParam Newlines | Comments EnumParam Newlines;

Start: start | Comments start | start Newlines | Comments start Newlines;
Stop: stop | stop Newlines;

Enum: EnumStart Start EnumParams Comments Stop;
Enum: EnumStart Start EnumParams Stop;

EnumStart: enum;
/.
    case $rule_number:
    {
        const QString name = captured().value(QLatin1String("name"));

        // new Class declaration
        m_astEnum = ASTEnum(name);
        m_astEnumValue = -1;
    }
    break;
./

EnumParams: DecoratedEnumParam | DecoratedEnumParam Comma EnumParams;

EnumParam: enum_param;
/.
    case $rule_number:
    {
        ASTEnumParam param;
        param.name = captured().value(QStringLiteral("name")).trimmed();
        QString value = captured().value(QStringLiteral("value"));
        value.remove(QLatin1Char('='));
        value = value.trimmed();
        if (value.isEmpty())
            param.value = ++m_astEnumValue;
        else if (value.startsWith(QLatin1String("0x"), Qt::CaseInsensitive))
            param.value = m_astEnumValue = value.toInt(0,16);
        else
            param.value = m_astEnumValue = value.toInt();
        if (param.value < 0) {
            m_astEnum.isSigned = true;
            if (m_astEnum.max < -param.value)
                m_astEnum.max = -param.value;
        } else if (m_astEnum.max < param.value)
            m_astEnum.max = param.value;
        m_astEnum.params << param;
    }
    break;
./

Prop: prop;
/.
    case $rule_number:
    {
        const QString args = captured().value(QLatin1String("args"));
        if (!parseProperty(m_astClass, args))
            return false;
    }
    break;
./

Signal: signal;
/.
    case $rule_number:
    {
        ASTFunction signal;
        signal.name = captured().value(QLatin1String("name")).trimmed();

        const QString argString = captured().value(QLatin1String("args")).trimmed();
        RepParser::TypeParser parseType;
        parseType.parseArguments(argString);
        parseType.appendParams(signal);
        m_astClass.signalsList << signal;
    }
    break;
./

Slot: slot;
/.
    case $rule_number:
    {
        QString returnTypeAndName = captured().value(QLatin1String("type")).trimmed();
        const QString argString = captured().value(QLatin1String("args")).trimmed();

        // compat code with old SLOT declaration: "SLOT(func(...))"
        const bool hasWhitespace = returnTypeAndName.indexOf(u' ') != -1;
        if (!hasWhitespace) {
            qWarning() << "[repc] - Adding 'void' for unspecified return type on" << qPrintable(returnTypeAndName);
            returnTypeAndName.prepend(QLatin1String("void "));
        }

        const int startOfFunctionName = returnTypeAndName.lastIndexOf(u' ') + 1;

        ASTFunction slot;
        slot.returnType = returnTypeAndName.mid(0, startOfFunctionName-1);
        slot.name = returnTypeAndName.mid(startOfFunctionName);

        RepParser::TypeParser parseType;
        parseType.parseArguments(argString);
        parseType.appendParams(slot);
        m_astClass.slotsList << slot;
    }
    break;
./

Model: model;
/.
    case $rule_number:
    {
        ASTModel model(m_astClass.properties.size());
        const QString name = captured().value(QLatin1String("name")).trimmed();
        const QString argString = captured().value(QLatin1String("args")).trimmed();

        if (!parseRoles(model, argString))
            return false;

        m_astClass.modelMetadata << model;
        m_astClass.properties << ASTProperty(QStringLiteral("QAbstractItemModel"), name, QStringLiteral("nullptr"), ASTProperty::SourceOnlySetter, false, true);
    }
    break;
./

ChildRep: childrep;
/.
case $rule_number:
{
    const QString name = captured().value(QLatin1String("name")).trimmed();
    const QString type = captured().value(QLatin1String("type")).trimmed();

    m_astClass.subClassPropertyIndices << m_astClass.properties.size();
    m_astClass.properties << ASTProperty(type, name, QStringLiteral("nullptr"), ASTProperty::SourceOnlySetter, false, true);
}
break;
./

ClassStart: class Newlines;
/.
    case $rule_number:
./
ClassStart: class;
/.
    case $rule_number:
    {
        const QString name = captured().value(QLatin1String("name"));

        // new Class declaration
        m_astClass = ASTClass(name);
    }
    break;
./

UseEnum: use_enum;
/.
    case $rule_number:
    {
        const QString name = captured().value(QLatin1String("name"));

        m_ast.enumUses.append(name);
    }
    break;
./

--Error conditions/messages
ClassType: ClassStart;
/.
    case $rule_number:
    {
        setErrorString(QStringLiteral("class: Cannot be nested"));
        return false;
    }
    break;
./
ClassType: Pod;
/.
    case $rule_number:
    {
        setErrorString(QStringLiteral("POD: Can only be used in global scope"));
        return false;
    }
    break;
./
ClassType: UseEnum;
/.
    case $rule_number:
    {
        setErrorString(QStringLiteral("USE_ENUM: Can only be used in global scope"));
        return false;
    }
    break;
./
Type: Signal;
/.
    case $rule_number:
    {
        setErrorString(QStringLiteral("SIGNAL: Can only be used in class scope"));
        return false;
    }
    break;
./
Type: Slot;
/.
    case $rule_number:
    {
        setErrorString(QStringLiteral("SLOT: Can only be used in class scope"));
        return false;
    }
    break;
./
Type: Prop;
/.
    case $rule_number:
    {
        setErrorString(QStringLiteral("PROP: Can only be used in class scope"));
        return false;
    }
    break;
./
Type: Model;
/.
    case $rule_number:
    {
        setErrorString(QStringLiteral("MODEL: Can only be used in class scope"));
        return false;
    }
    break;
./

Type: ChildRep;
/.
    case $rule_number:
    {
        setErrorString(QStringLiteral("CLASS: Can only be used in class scope"));
        return false;
    }
    break;
./

/.
    } // switch
    return true;
}
./
