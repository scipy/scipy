/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtXmlPatterns module of the Qt Toolkit.
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

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#ifndef Patternist_DocumentFN_H
#define Patternist_DocumentFN_H

#include <private/qfunctioncall_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{

    /**
     * @short Implements XSL-T's function <tt>fn:document()</tt>.
     *
     * @c fn:document() has no evaluation functions, because it rewrites
     * itself to a set of expressions that is the implementation.
     *
     * The two-argument version:
     *
     * <tt>document($uris as item()*, $baseURINode as node()) as node()*</tt>
     *
     * is rewritten into:
     *
     * <tt>for $uri in distinct-values($uris)
     * return doc(resolve-uri($uri, base-uri($baseURINode)))</tt>
     *
     * and the single version argument:
     *
     * <tt>document($uris as item()*) as node()*</tt>
     *
     * is rewritten into:
     *
     * <tt>for $uri in distinct-values($uris)
     * return doc($uri)</tt>
     *
     * The distinct-values() call ensures the node deduplication and sorting,
     * although it fails in the case that URIs resolve/directs in some way to
     * the same document. Some of those cases can be solved by wrapping the
     * whole expression with a node deduplication(conceptually the-for-loop/.).
     * One advantage with distinct-values() over generating traditional node
     * sorting/deduplication code is that the latter contains node sorting
     * which is uecessary and can be hard to analyze away. distinct-values()
     * doesn't have this problem due to its narrower task..
     *
     * This works without problems, assuming XTRE1160 is not raised and that
     * the recover action instead is ignore. In the case XTRE1160 is raised,
     * one must cater for this.
     *
     * In addition to this, both signatures has its first argument changed to
     * type <tt>xs:string*</tt>, in order to generate atomization code.
     *
     * One notable thing is that the expression for $baseURINode, is moved
     * inside a loop, and will be evaluated repeatedly, unless moved out as
     * part of optimization.
     *
     * @ingroup Patternist_functions
     * @author Frans Englich <frans.englich@nokia.com>
     * @since 4.5
     */
    class DocumentFN : public FunctionCall
    {
    public:
        virtual Expression::Ptr typeCheck(const StaticContext::Ptr &context,
                                          const SequenceType::Ptr &reqType);
    };
}

QT_END_NAMESPACE

#endif
