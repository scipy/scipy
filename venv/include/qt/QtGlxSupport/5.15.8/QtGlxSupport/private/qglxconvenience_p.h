/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the plugins of the Qt Toolkit.
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

#ifndef QGLXCONVENIENCE_H
#define QGLXCONVENIENCE_H

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

#include <QSurfaceFormat>
#include <QVector>

#include <X11/Xlib.h>
#include <GL/glx.h>

enum QGlxFlags
{
    QGLX_SUPPORTS_SRGB = 0x01
};

QVector<int> qglx_buildSpec(const QSurfaceFormat &format, int drawableBit = GLX_WINDOW_BIT, int flags = 0);
XVisualInfo *qglx_findVisualInfo(Display *display, int screen, QSurfaceFormat *format, int drawableBit = GLX_WINDOW_BIT, int flags = 0);
GLXFBConfig qglx_findConfig(Display *display, int screen, QSurfaceFormat format, bool highestPixelFormat = false, int drawableBit = GLX_WINDOW_BIT, int flags = 0);
void qglx_surfaceFormatFromGLXFBConfig(QSurfaceFormat *format, Display *display, GLXFBConfig config, int flags = 0);
void qglx_surfaceFormatFromVisualInfo(QSurfaceFormat *format, Display *display, XVisualInfo *visualInfo, int flags = 0);
bool qglx_reduceFormat(QSurfaceFormat *format);

#endif // QGLXCONVENIENCE_H
