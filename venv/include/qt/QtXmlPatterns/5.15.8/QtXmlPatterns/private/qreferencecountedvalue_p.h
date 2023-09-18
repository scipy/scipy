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

#ifndef QPatternist_ReferenceCountedValue_p_h
#define QPatternist_ReferenceCountedValue_p_h

QT_BEGIN_NAMESPACE

namespace QPatternist
{
/*!
   \class ReferenceCountedValue
   \internal
   \since 4.4
   \brief A template class that reference counts a value.

   This class is useful when an instance needs to have ownership semantics
   as if it was value based. A typical examples is a QObject pointer, which
   doesn't have a single owner.

   This is achieved through storing a copy of the object as
   a member inside ReferenceCountedValue, which never is copied. It will
   stay in scope until the last reference to the ReferenceCountedValue instance
   is removed, and subsequently ReferenceCountedValue is deleted and hence also
   the contained value. One should use ReferenceCountedValue by passing around
   copies of Ptr, which is a typedef for the QExplicitlySharedDataPointer
   smart pointer.
*/
    template<typename T>
    class ReferenceCountedValue : public QSharedData
    {
    public:
        typedef QExplicitlySharedDataPointer<ReferenceCountedValue<T> > Ptr;

        inline ReferenceCountedValue(T *const v) : value(v)
        {
        }

        inline ~ReferenceCountedValue()
        {
            delete value;
        }

        T *const value;
    private:
        /*!
          Disabled, no implementation provided.
         */
        inline ReferenceCountedValue();
        Q_DISABLE_COPY(ReferenceCountedValue)
    };
}

QT_END_NAMESPACE

#endif
