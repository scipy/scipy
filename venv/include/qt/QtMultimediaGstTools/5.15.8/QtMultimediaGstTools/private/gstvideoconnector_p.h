/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Toolkit.
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

#ifndef QGSTVIDEOCONNECTOR_H
#define QGSTVIDEOCONNECTOR_H

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

#include <private/qgsttools_global_p.h>

#include <gst/gst.h>

G_BEGIN_DECLS

#define GST_TYPE_VIDEO_CONNECTOR \
  (gst_video_connector_get_type())
#define GST_VIDEO_CONNECTOR(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_VIDEO_CONNECTOR, GstVideoConnector))
#define GST_VIDEO_CONNECTOR_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_CAST ((klass), GST_TYPE_VIDEO_CONNECTOR, GstVideoConnectorClass))
#define GST_IS_VIDEO_CONNECTOR(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_VIDEO_CONNECTOR))
#define GST_IS_VIDEO_CONNECTOR_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_TYPE ((klass), GST_TYPE_VIDEO_CONNECTOR))

typedef struct _GstVideoConnector GstVideoConnector;
typedef struct _GstVideoConnectorClass GstVideoConnectorClass;

struct Q_GSTTOOLS_EXPORT _GstVideoConnector {
  GstElement element;

  GstPad *srcpad;
  GstPad *sinkpad;

  gboolean relinked;
  gboolean failedSignalEmited;
  GstSegment segment;
  GstBuffer *latest_buffer;
};

struct Q_GSTTOOLS_EXPORT _GstVideoConnectorClass {
  GstElementClass parent_class;

  /* action signal to resend new segment */
  void (*resend_new_segment) (GstElement * element, gboolean emitFailedSignal);
};

GType Q_GSTTOOLS_EXPORT gst_video_connector_get_type (void);

G_END_DECLS

#endif

