/****************************************************************************
**
** Copyright (C) 2018 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Virtual Keyboard module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QVIRTUALKEYBOARD_STATICPLUGIN_P_H
#define QVIRTUALKEYBOARD_STATICPLUGIN_P_H

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

#if defined(QT_STATICPLUGIN)
#include <QtPlugin>
#include <QPluginLoader>

QT_BEGIN_NAMESPACE

// This macro is similar to Q_IMPORT_PLUGIN, except it does not
// register duplicate entries as static plugins.
// The check is required since the application may already have
// initialized the plugin by its own dependencies.
#define Q_VKB_IMPORT_PLUGIN(PLUGIN) \
        extern const QT_PREPEND_NAMESPACE(QStaticPlugin) qt_static_plugin_##PLUGIN(); \
        class Static##PLUGIN##PluginInstance{ \
        public: \
                Static##PLUGIN##PluginInstance() { \
                    if (!QPluginLoader::staticInstances().contains(qt_static_plugin_##PLUGIN().instance())) \
                        qRegisterStaticPluginFunction(qt_static_plugin_##PLUGIN()); \
                } \
        }; \
       static Static##PLUGIN##PluginInstance static##PLUGIN##Instance;

QT_END_NAMESPACE

#endif

#endif // QVIRTUALKEYBOARD_STATICPLUGIN_P_H
