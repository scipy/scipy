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

#ifndef QT3DCORE_GLOBAL_H
#define QT3DCORE_GLOBAL_H

#include <QtCore/qglobal.h>
#include <Qt3DCore/qt3dcore-config.h>

QT_BEGIN_NAMESPACE

#if defined(QT_SHARED) || !defined(QT_STATIC)
#  if defined(QT_BUILD_3DCORE_LIB)
#    define Q_3DCORESHARED_EXPORT Q_DECL_EXPORT
#  else
#    define Q_3DCORESHARED_EXPORT Q_DECL_IMPORT
#  endif
#else
#  define Q_3DCORESHARED_EXPORT
#endif

#ifdef BUILD_QT3D_MODULE
#define Q3D_DECL_DEPRECATED
#else
#define Q3D_DECL_DEPRECATED Q_DECL_DEPRECATED
#endif

#define QT3D_DECLARE_TYPEINFO(NS, Class, Flags) \
    } /* NS */ \
    Q_DECLARE_TYPEINFO(NS :: Class, Flags); \
    namespace NS { \
    /*end*/

#define QT3D_DECLARE_TYPEINFO_2(OuterNS, InnerNS, Class, Flags) \
    } /* InnerNS */ \
    QT3D_DECLARE_TYPEINFO(OuterNS, InnerNS :: Class, Flags) \
    namespace InnerNS { \
    /*end*/

#define QT3D_DECLARE_TYPEINFO_3(OuterNS, InnerNS, InnerMostNS, Class, Flags) \
    } /* InnerMostNS */ \
    QT3D_DECLARE_TYPEINFO_2(OuterNS, InnerNS, InnerMostNS :: Class, Flags) \
    namespace InnerMostNS { \
    /*end*/

#define QT3D_DECLARE_SHARED(NS, Class) \
    inline void swap(Class &lhs, Class &rhs) \
        Q_DECL_NOEXCEPT_EXPR(noexcept(lhs.swap(rhs))) \
    { lhs.swap(rhs); } \
    QT3D_DECLARE_TYPEINFO(NS, Class, Q_MOVABLE_TYPE) \
    /*end*/

#define QT3D_DECLARE_SHARED_2(OuterNS, InnerNS, Class) \
    inline void swap(Class &lhs, Class &rhs) \
        Q_DECL_NOEXCEPT_EXPR(noexcept(lhs.swap(rhs))) \
    { lhs.swap(rhs); } \
    QT3D_DECLARE_TYPEINFO_2(OuterNS, InnerNS, Class, Q_MOVABLE_TYPE) \
    /*end*/

QT_END_NAMESPACE

#endif // QT3DCORE_GLOBAL_H
