/****************************************************************************
**
** Copyright (C) 2014 Klaralvdalens Datakonsult AB (KDAB).
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
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

#ifndef QT3DCORE_QNODEVISITOR_P_H
#define QT3DCORE_QNODEVISITOR_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of other Qt classes.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <Qt3DCore/qentity.h>
#include <Qt3DCore/qnode.h>
#include <Qt3DCore/qt3dcore_global.h>

QT_BEGIN_NAMESPACE

namespace Qt3DCore
{

class Q_3DCORESHARED_EXPORT QNodeVisitor
{
public:
    QNodeVisitor();
    virtual ~QNodeVisitor();

    template<typename NodeVisitorFunc>
    void traverse(QNode *rootNode_, NodeVisitorFunc fN)
    {
        startTraversing(rootNode_, createFunctor(fN));
    }

    template<typename Obj, typename NodeVisitorFunc>
    void traverse(QNode *rootNode_, Obj *instance, NodeVisitorFunc fN)
    {
        startTraversing(rootNode_, createFunctor(instance, fN));
    }

    template<typename NodeVisitorFunc, typename EntityVisitorFunc>
    void traverse(QNode *rootNode_, NodeVisitorFunc fN, EntityVisitorFunc fE)
    {
        startTraversing(rootNode_, createFunctor(fN), createFunctor(fE));
    }

    template<typename Obj, typename NodeVisitorFunc, typename EntityVisitorFunc>
    void traverse(QNode *rootNode_, Obj *instance, NodeVisitorFunc fN, EntityVisitorFunc fE)
    {
        startTraversing(rootNode_, createFunctor(instance, fN), createFunctor(instance, fE));
    }

    QNode *rootNode() const;
    QNode *currentNode() const;
    void setPath(QNodeVector path);
    QNodeVector path() const;
    void append(QNode *n);
    void pop_back();

private:
    Q_DISABLE_COPY(QNodeVisitor)
    QNodeVector m_path;

    template<typename NodeVisitorFunctor>
    void startTraversing(QNode *rootNode_, NodeVisitorFunctor fN)
    {
        setPath(QNodeVector() << rootNode_);
        if (rootNode_)
            visitNode(rootNode_, fN);
    }

    template<typename NodeVisitorFunctor, typename EntityVisitorFunctor>
    void startTraversing(QNode *rootNode_, NodeVisitorFunctor fN, EntityVisitorFunctor fE)
    {
        setPath(QNodeVector() << rootNode_);
        QEntity* rootEntity = qobject_cast<QEntity *>(rootNode_);

        if (rootEntity)
            visitEntity(rootEntity, fN, fE);
        else if (rootNode_)
            visitNode(rootNode_, fN, fE);
    }

    template<typename NodeVisitorFunctor>
    void visitNode(QNode *nd, NodeVisitorFunctor &fN)
    {
        fN(nd);
        traverseChildren(fN);
    }

    template<typename NodeVisitorFunctor, typename EntityVisitorFunctor>
    void visitNode(QNode *nd, NodeVisitorFunctor &fN, EntityVisitorFunctor &fE)
    {
        fN(nd);
        traverseChildren(fN, fE);
    }

    template<typename NodeVisitorFunctor, typename EntityVisitorFunctor>
    void visitEntity(QEntity *ent, NodeVisitorFunctor &fN, EntityVisitorFunctor &fE)
    {
        fE(ent);
        visitNode(ent, fN, fE);
    }

    template<typename NodeVisitorFunctor, typename EntityVisitorFunctor>
    void traverseChildren(NodeVisitorFunctor &fN, EntityVisitorFunctor &fE)
    {
        for (QObject *n : currentNode()->children()) {
            QNode *node = qobject_cast<QNode *>(n);
            if (node != nullptr)
                outerVisitNode(node, fN, fE);
        } // of children iteration
    }

    template<typename NodeVisitorFunctor>
    void traverseChildren(NodeVisitorFunctor &fN)
    {
        for (QObject *n : currentNode()->children()) {
            QNode *node = qobject_cast<QNode *>(n);
            if (node != nullptr)
                outerVisitNode(node, fN);
        } // of children iteration
    }

    template<typename NodeVisitorFunctor, typename EntityVisitorFunctor>
    void outerVisitNode(QNode *n, NodeVisitorFunctor &fN, EntityVisitorFunctor &fE)
    {
        append(n);
        QEntity* e = qobject_cast<QEntity *>(n);
        if (e) {
            visitEntity(e, fN, fE);
        } else {
            visitNode(n, fN, fE);
        }
        pop_back();
    }

    template<typename NodeVisitorFunctor>
    void outerVisitNode(QNode *n, NodeVisitorFunctor &fN)
    {
        append(n);
        visitNode(n, fN);
        pop_back();
    }

    template <typename ReturnType, typename NodeType>
    class FunctionFunctor {
    public:
        typedef ReturnType (*functionPtr)(NodeType);

        FunctionFunctor(functionPtr fPtr)
            : m_functionPointer(fPtr)
        {}

        void operator()(NodeType node)
        {
            (*m_functionPointer)(node);
        }

    private:
        functionPtr m_functionPointer;
    };

    template <typename C, typename ReturnType, typename NodeType>
    class MemberFunctionFunctor {
    public:
        typedef ReturnType (C::*functionPtr)(NodeType);

        MemberFunctionFunctor(C* instance, functionPtr fPtr)
            : m_instance(instance)
            , m_functionPointer(fPtr)
        {}

        void operator()(NodeType node)
        {
            (*m_instance.*m_functionPointer)(node);
        }

    private:
        C *m_instance;
        functionPtr m_functionPointer;
    };

    template <typename C, typename ReturnType, typename NodeType>
    class ConstMemberFunctionFunctor {
    public:
        typedef ReturnType (C::*functionPtr)(NodeType) const;

        ConstMemberFunctionFunctor(C* instance, functionPtr fPtr)
            : m_instance(instance)
            , m_functionPointer(fPtr)
        {}

        void operator()(NodeType node) const
        {
            (*m_instance.*m_functionPointer)(node);
        }

    private:
        C *m_instance;
        functionPtr m_functionPointer;
    };

    template <typename T>
    const T& createFunctor(const T& t)
    {
        return t;
    }

    template <typename ReturnType, typename NodeType>
    FunctionFunctor<ReturnType, NodeType> createFunctor(ReturnType (*fPtr)(NodeType))
    {
        return FunctionFunctor<ReturnType, NodeType>(fPtr);
    }

    template <typename C, typename ReturnType, typename NodeType>
    MemberFunctionFunctor<C, ReturnType, NodeType> createFunctor(C *instance, ReturnType (C::*fPtr)(NodeType))
    {
        return MemberFunctionFunctor<C, ReturnType, NodeType>(instance, fPtr);
    }

    template <typename C, typename ReturnType, typename NodeType>
    ConstMemberFunctionFunctor<C, ReturnType, NodeType> createFunctor(C *instance, ReturnType (C::*fPtr)(NodeType) const)
    {
        return ConstMemberFunctionFunctor<C, ReturnType, NodeType>(instance, fPtr);
    }
};

} // namespace Qt3DCore

QT_END_NAMESPACE

#endif // QT3DCORE_QNODEVISITOR_P_H
