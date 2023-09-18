
#ifndef QVULKANFUNCTIONS_P_H
#define QVULKANFUNCTIONS_P_H

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

#include "qvulkanfunctions.h"

QT_BEGIN_NAMESPACE

class QVulkanInstance;

class QVulkanFunctionsPrivate
{
public:
    QVulkanFunctionsPrivate(QVulkanInstance *inst);

    PFN_vkVoidFunction m_funcs[14];
};

class QVulkanDeviceFunctionsPrivate
{
public:
    QVulkanDeviceFunctionsPrivate(QVulkanInstance *inst, VkDevice device);

    PFN_vkVoidFunction m_funcs[120];
};

QT_END_NAMESPACE

#endif // QVULKANFUNCTIONS_P_H
