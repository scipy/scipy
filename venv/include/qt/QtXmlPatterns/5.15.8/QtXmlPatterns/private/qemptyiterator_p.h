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

#ifndef Patternist_EmptyIterator_H
#define Patternist_EmptyIterator_H

#include <private/qabstractxmlforwarditerator_p.h>
#include <private/qprimitives_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{

    /**
     * @short An QAbstractXmlForwardIterator which always is empty.
     *
     * EmptyIterator is an QAbstractXmlForwardIterator over the type @c T, which always is empty. Other
     * iterators can also be empty(or, at least behave as they are empty), but this
     * class is special designed for this purpose and is therefore fast.
     *
     * EmptyIterator's constructor is protected, instances is retrieved from CommonValues.
     *
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_iterators
     */
    template<typename T> class EmptyIterator : public QAbstractXmlForwardIterator<T>
    {
    public:
        /**
         * @returns always a default constructed value, T().
         */
        virtual T next()
        {
            return T();
        }

        /**
         * @returns always a default constructed value, T().
         */
        virtual T current() const
        {
            return T();
        }

        /**
         * @returns always 0.
         */
        virtual xsInteger position() const
        {
            return 0;
        }

        /**
         * @returns always @c this, the reverse of <tt>()</tt> is <tt>()</tt>.
         */
        virtual typename QAbstractXmlForwardIterator<T>::Ptr toReversed()
        {
            return typename QAbstractXmlForwardIterator<T>::Ptr(const_cast<EmptyIterator<T> *>(this));
        }

        /**
         * @returns always 0
         */
        virtual xsInteger count()
        {
            return 0;
        }

        /**
         * @returns @c this
         */
        virtual typename QAbstractXmlForwardIterator<T>::Ptr copy() const
        {
            return typename QAbstractXmlForwardIterator<T>::Ptr(const_cast<EmptyIterator *>(this));
        }

    protected:
        friend class CommonValues;
    };

    template<typename T>
    static inline
    typename QAbstractXmlForwardIterator<T>::Ptr
    makeEmptyIterator()
    {
        return typename QAbstractXmlForwardIterator<T>::Ptr(new EmptyIterator<T>());
    }

}

QT_END_NAMESPACE

#endif
