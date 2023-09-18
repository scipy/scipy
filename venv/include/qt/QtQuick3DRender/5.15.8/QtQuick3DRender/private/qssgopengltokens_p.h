/****************************************************************************
**
** Copyright (C) 2008-2012 NVIDIA Corporation.
** Copyright (C) 2019 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of Qt Quick 3D.
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


#ifndef QSSGOPENGLTOKENS_H
#define QSSGOPENGLTOKENS_H

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

#include <QtQuick3DRender/private/qtquick3drenderglobal_p.h>
#include <qopengl.h>

#ifndef GL_TEXTURE_2D_MULTISAMPLE
#define GL_TEXTURE_2D_MULTISAMPLE 0x9100
#endif

#ifndef GL_IMAGE_2D
#define GL_IMAGE_2D 0x904D
#endif

#ifndef GL_MULTISAMPLE
#define GL_MULTISAMPLE 0x809D
#endif

#ifndef GL_COLOR_ATTACHMENT1
#define GL_COLOR_ATTACHMENT1 0x8CE1
#define GL_COLOR_ATTACHMENT2 0x8CE2
#define GL_COLOR_ATTACHMENT3 0x8CE3
#define GL_COLOR_ATTACHMENT4 0x8CE4
#define GL_COLOR_ATTACHMENT5 0x8CE5
#define GL_COLOR_ATTACHMENT6 0x8CE6
#define GL_COLOR_ATTACHMENT7 0x8CE7
#endif

#ifndef GL_RED
#define GL_RED 0x1903
#define GL_GREEN 0x1904
#define GL_BLUE 0x1905
#endif

#ifndef GL_PATCHES
#define GL_PATCHES 0x000E
#endif

#ifndef GL_READ_ONLY
#define GL_READ_ONLY 0x88B8
#define GL_WRITE_ONLY 0x88B9
#define GL_READ_WRITE 0x88BA
#endif

#ifndef GL_SHADER_STORAGE_BUFFER
#define GL_SHADER_STORAGE_BUFFER 0x90D2
#endif

#ifndef GL_ATOMIC_COUNTER_BUFFER
#define GL_ATOMIC_COUNTER_BUFFER 0x92C0
#endif

#ifndef GL_VERTEX_SHADER_BIT
#define GL_VERTEX_SHADER_BIT 0x00000001
#endif

#ifndef GL_FRAGMENT_SHADER_BIT
#define GL_FRAGMENT_SHADER_BIT 0x00000002
#endif

#ifndef GL_GEOMETRY_SHADER_BIT
#define GL_GEOMETRY_SHADER_BIT 0x00000004
#endif

#ifndef GL_TESS_CONTROL_SHADER_BIT
#define GL_TESS_CONTROL_SHADER_BIT 0x00000008
#endif

#ifndef GL_TESS_EVALUATION_SHADER_BIT
#define GL_TESS_EVALUATION_SHADER_BIT 0x00000010
#endif

#ifndef GL_UNIFORM_BUFFER
#define GL_UNIFORM_BUFFER                 0x8A11
#endif

#ifndef GL_ANY_SAMPLES_PASSED
#define GL_ANY_SAMPLES_PASSED             0x8C2F
#endif

#ifndef GL_QUERY_RESULT_AVAILABLE
#define GL_QUERY_RESULT_AVAILABLE         0x8867
#endif

#ifndef GL_QUERY_RESULT
#define GL_QUERY_RESULT                   0x8866
#endif

#ifndef GL_SYNC_GPU_COMMANDS_COMPLETE
#define GL_SYNC_GPU_COMMANDS_COMPLETE     0x9117
#endif

#ifndef GL_R8
#define GL_R8                             0x8229
#endif

#ifndef GL_RG8
#define GL_RG8                            0x822B
#endif

#ifndef GL_RG
#define GL_RG                             0x8227
#endif

#ifndef GL_RGBA8
#define GL_RGBA8                          0x8058
#endif

#ifndef GL_RGB8
#define GL_RGB8                           0x8051
#endif

#ifndef GL_R16F
#define GL_R16F                           0x822D
#endif

#ifndef GL_R32F
#define GL_R32F                           0x822E
#endif

#ifndef GL_RGBA16F
#define GL_RGBA16F                        0x881A
#endif

#ifndef GL_HALF_FLOAT
#define GL_HALF_FLOAT                     0x140B
#endif

#ifndef GL_RED_INTEGER
#define GL_RED_INTEGER                    0x8D94
#endif

#ifndef GL_R32UI
#define GL_R32UI                          0x8236
#endif

#ifndef GL_BGRA
#define GL_BGRA                           0x80E1
#endif

#ifndef GL_R16
#define GL_R16                            0x822A
#endif

#ifndef GL_RED
#define GL_RED                            0x1903
#endif

#ifndef GL_RGBA32F
#define GL_RGBA32F                        0x8814
#endif

#ifndef GL_DEPTH_COMPONENT16
#define GL_DEPTH_COMPONENT16              0x81A5
#endif

#ifndef GL_DEPTH_COMPONENT24
#define GL_DEPTH_COMPONENT24              0x81A6
#endif

#ifndef GL_DEPTH_COMPONENT32F
#define GL_DEPTH_COMPONENT32F             0x8CAC
#endif

#ifndef GL_STENCIL_INDEX
#define GL_STENCIL_INDEX                  0x1901
#endif

#ifndef GL_STENCIL_INDEX8
#define GL_STENCIL_INDEX8                 0x8D48
#endif

#ifndef GL_DEPTH24_STENCIL8
#define GL_DEPTH24_STENCIL8               0x88F0
#endif

#ifndef GL_DEPTH_STENCIL_ATTACHMENT
#define GL_DEPTH_STENCIL_ATTACHMENT       0x821A
#endif

#ifndef GL_DEPTH_STENCIL
#define GL_DEPTH_STENCIL                  0x84F9
#endif

#ifndef GL_PRIMITIVE_RESTART_FIXED_INDEX
#define GL_PRIMITIVE_RESTART_FIXED_INDEX  0x8D69
#endif

#ifndef GL_FRAMEBUFFER_SRGB
#define GL_FRAMEBUFFER_SRGB 0x8DB9
#endif

#ifndef GL_READ_FRAMEBUFFER
#define GL_READ_FRAMEBUFFER               0x8CA8
#endif

#ifndef GL_DRAW_FRAMEBUFFER
#define GL_DRAW_FRAMEBUFFER               0x8CA9
#endif

#ifndef GL_MAX_DRAW_BUFFERS
#define GL_MAX_DRAW_BUFFERS               0x8824
#endif

#ifndef GL_TEXTURE_COMPARE_MODE
#define GL_TEXTURE_COMPARE_MODE           0x884C
#endif

#ifndef GL_COMPARE_REF_TO_TEXTURE
#define GL_COMPARE_REF_TO_TEXTURE         0x884E
#endif

#ifndef GL_TEXTURE_COMPARE_FUNC
#define GL_TEXTURE_COMPARE_FUNC           0x884D
#endif

#ifndef GL_MAX_SAMPLES
#define GL_MAX_SAMPLES                    0x8D57
#endif

#ifndef GL_SHADER_STORAGE_BUFFER
#define GL_SHADER_STORAGE_BUFFER          0x90D2
#endif

#ifndef GL_READ_ONLY
#define GL_READ_ONLY                      0x88B8
#endif

#ifndef GL_WRITE_ONLY
#define GL_WRITE_ONLY                     0x88B9
#endif

#ifndef GL_READ_WRITE
#define GL_READ_WRITE                     0x88BA
#endif

#ifndef GL_COMPUTE_SHADER
#define GL_COMPUTE_SHADER                 0x91B9
#endif

#ifndef GL_VERTEX_PROGRAM_POINT_SIZE
#define GL_VERTEX_PROGRAM_POINT_SIZE      0x8642
#endif

#ifndef GL_POINT_SPRITE
#define GL_POINT_SPRITE                   0x8861
#endif

#ifndef GL_RG16F
#define GL_RG16F                          0x822F
#endif

#ifndef GL_RG32F
#define GL_RG32F                          0x8230
#endif

#ifndef GL_RGB32F
#define GL_RGB32F                         0x8815
#endif

#ifndef GL_R11F_G11F_B10F
#define GL_R11F_G11F_B10F                 0x8C3A
#endif

#ifndef GL_UNSIGNED_INT_10F_11F_11F_REV
#define GL_UNSIGNED_INT_10F_11F_11F_REV   0x8C3B
#endif

#ifndef GL_RGB9_E5
#define GL_RGB9_E5                        0x8C3D
#endif

#ifndef GL_UNSIGNED_INT_5_9_9_9_REV
#define GL_UNSIGNED_INT_5_9_9_9_REV       0x8C3E
#endif

#ifndef GL_SRGB8
#define GL_SRGB8                          0x8C41
#endif

#ifndef GL_SRGB8_ALPHA8
#define GL_SRGB8_ALPHA8                   0x8C43
#endif

#ifndef GL_UNSIGNED_INT_10_10_10_2
#define GL_UNSIGNED_INT_10_10_10_2        0x8DF6
#endif

#ifndef GL_UNSIGNED_INT_24_8
#define GL_UNSIGNED_INT_24_8              0x84FA
#endif

#ifndef GL_R32I
#define GL_R32I                           0x8235
#endif

#ifndef GL_MAP_INVALIDATE_BUFFER_BIT
#define GL_MAP_INVALIDATE_BUFFER_BIT      0x0008
#endif

#ifndef GL_MAP_INVALIDATE_RANGE_BIT
#define GL_MAP_INVALIDATE_RANGE_BIT       0x0004
#endif

#ifndef GL_UNSIGNED_INT_VEC2
#define GL_UNSIGNED_INT_VEC2              0x8DC6
#endif

#ifndef GL_UNSIGNED_INT_VEC3
#define GL_UNSIGNED_INT_VEC3              0x8DC7
#endif

#ifndef GL_UNSIGNED_INT_VEC4
#define GL_UNSIGNED_INT_VEC4              0x8DC8
#endif

#ifndef GL_SAMPLER_2D_SHADOW
#define GL_SAMPLER_2D_SHADOW              0x8B62
#endif

#ifndef GL_TEXTURE_BASE_LEVEL
#define GL_TEXTURE_BASE_LEVEL             0x813C
#endif

#ifndef GL_TEXTURE_MAX_LEVEL
#define GL_TEXTURE_MAX_LEVEL              0x813D
#endif

#ifndef GL_TEXTURE_SWIZZLE_R
#define GL_TEXTURE_SWIZZLE_R              0x8E42
#endif

#ifndef GL_TEXTURE_SWIZZLE_G
#define GL_TEXTURE_SWIZZLE_G              0x8E43
#endif

#ifndef GL_TEXTURE_SWIZZLE_B
#define GL_TEXTURE_SWIZZLE_B              0x8E44
#endif

#ifndef GL_TEXTURE_SWIZZLE_A
#define GL_TEXTURE_SWIZZLE_A              0x8E45
#endif

#ifndef GL_TEXTURE_MIN_LOD
#define GL_TEXTURE_MIN_LOD                0x813A
#endif

#ifndef GL_TEXTURE_MAX_LOD
#define GL_TEXTURE_MAX_LOD                0x813B
#endif

#ifndef GL_TEXTURE_WRAP_R
#define GL_TEXTURE_WRAP_R                 0x8072
#endif

#ifndef GL_INVALID_INDEX
#define GL_INVALID_INDEX                  0xFFFFFFFFu
#endif

#ifndef GL_UNIFORM_BLOCK_DATA_SIZE
#define GL_UNIFORM_BLOCK_DATA_SIZE        0x8A40
#endif

#ifndef GL_UNIFORM_BLOCK_NAME_LENGTH
#define GL_UNIFORM_BLOCK_NAME_LENGTH      0x8A41
#endif

#ifndef GL_UNIFORM_BLOCK_ACTIVE_UNIFORMS
#define GL_UNIFORM_BLOCK_ACTIVE_UNIFORMS  0x8A42
#endif

#ifndef GL_UNIFORM_BLOCK_ACTIVE_UNIFORM_INDICES
#define GL_UNIFORM_BLOCK_ACTIVE_UNIFORM_INDICES 0x8A43
#endif

#ifndef GL_UNIFORM_TYPE
#define GL_UNIFORM_TYPE                   0x8A37
#endif

#ifndef GL_UNIFORM_SIZE
#define GL_UNIFORM_SIZE                   0x8A38
#endif

#ifndef GL_UNIFORM_NAME_LENGTH
#define GL_UNIFORM_NAME_LENGTH            0x8A39
#endif

#ifndef GL_UNIFORM_BLOCK_INDEX
#define GL_UNIFORM_BLOCK_INDEX            0x8A3A
#endif

#ifndef GL_UNIFORM_OFFSET
#define GL_UNIFORM_OFFSET                 0x8A3B
#endif

#ifndef GL_ACTIVE_UNIFORM_BLOCKS
#define GL_ACTIVE_UNIFORM_BLOCKS          0x8A36
#endif

#ifndef GL_FRAMEBUFFER_ATTACHMENT_DEPTH_SIZE
#define GL_FRAMEBUFFER_ATTACHMENT_DEPTH_SIZE 0x8216
#endif

#ifndef GL_FRAMEBUFFER_ATTACHMENT_STENCIL_SIZE
#define GL_FRAMEBUFFER_ATTACHMENT_STENCIL_SIZE 0x8217
#endif

#ifndef GL_NUM_EXTENSIONS
#define GL_NUM_EXTENSIONS                 0x821D
#endif

#ifndef GL_TIMEOUT_IGNORED
#define GL_TIMEOUT_IGNORED                0xFFFFFFFFFFFFFFFFull
#endif

#ifndef GL_MAX_ARRAY_TEXTURE_LAYERS
#define GL_MAX_ARRAY_TEXTURE_LAYERS       0x88FF
#endif

#ifndef GL_MAX_UNIFORM_BUFFER_BINDINGS
#define GL_MAX_UNIFORM_BUFFER_BINDINGS    0x8A2F
#endif

#ifndef GL_MAX_UNIFORM_BLOCK_SIZE
#define GL_MAX_UNIFORM_BLOCK_SIZE         0x8A30
#endif

#ifndef GL_NVIDIA_PLATFORM_BINARY_NV
#define GL_NVIDIA_PLATFORM_BINARY_NV 0x890B
#endif

#ifndef GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS
#define GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS 0x8CD9
#endif


#ifndef GL_NV_blend_equation_advanced
#define GL_NV_blend_equation_advanced 1
#define GL_BLEND_OVERLAP_NV 0x9281
#define GL_BLEND_PREMULTIPLIED_SRC_NV 0x9280
#define GL_BLUE_NV 0x1905
#define GL_COLORBURN_NV 0x929A
#define GL_COLORDODGE_NV 0x9299
#define GL_CONJOINT_NV 0x9284
#define GL_CONTRAST_NV 0x92A1
#define GL_DARKEN_NV 0x9297
#define GL_DIFFERENCE_NV 0x929E
#define GL_DISJOINT_NV 0x9283
#define GL_DST_ATOP_NV 0x928F
#define GL_DST_IN_NV 0x928B
#define GL_DST_NV 0x9287
#define GL_DST_OUT_NV 0x928D
#define GL_DST_OVER_NV 0x9289
#define GL_EXCLUSION_NV 0x92A0
#define GL_GREEN_NV 0x1904
#define GL_HARDLIGHT_NV 0x929B
#define GL_HARDMIX_NV 0x92A9
#define GL_HSL_COLOR_NV 0x92AF
#define GL_HSL_HUE_NV 0x92AD
#define GL_HSL_LUMINOSITY_NV 0x92B0
#define GL_HSL_SATURATION_NV 0x92AE
#define GL_INVERT_OVG_NV 0x92B4
#define GL_INVERT_RGB_NV 0x92A3
#define GL_LIGHTEN_NV 0x9298
#define GL_LINEARBURN_NV 0x92A5
#define GL_LINEARDODGE_NV 0x92A4
#define GL_LINEARLIGHT_NV 0x92A7
#define GL_MINUS_CLAMPED_NV 0x92B3
#define GL_MINUS_NV 0x929F
#define GL_MULTIPLY_NV 0x9294
#define GL_OVERLAY_NV 0x9296
#define GL_PINLIGHT_NV 0x92A8
#define GL_PLUS_CLAMPED_ALPHA_NV 0x92B2
#define GL_PLUS_CLAMPED_NV 0x92B1
#define GL_PLUS_DARKER_NV 0x9292
#define GL_PLUS_NV 0x9291
#define GL_RED_NV 0x1903
#define GL_SCREEN_NV 0x9295
#define GL_SOFTLIGHT_NV 0x929C
#define GL_SRC_ATOP_NV 0x928E
#define GL_SRC_IN_NV 0x928A
#define GL_SRC_NV 0x9286
#define GL_SRC_OUT_NV 0x928C
#define GL_SRC_OVER_NV 0x9288
#define GL_UNCORRELATED_NV 0x9282
#define GL_VIVIDLIGHT_NV 0x92A6
#define GL_XOR_NV 0x1506
#endif /* GL_NV_blend_equation_advanced */

#ifndef GL_SHADER_STORAGE_BUFFER
#define GL_SHADER_STORAGE_BUFFER 0x90D2
#endif

#ifndef GL_ATOMIC_COUNTER_BUFFER
#define GL_ATOMIC_COUNTER_BUFFER 0x92C0
#endif

#ifndef GL_ALL_BARRIER_BITS
#define GL_VERTEX_ATTRIB_ARRAY_BARRIER_BIT 0x00000001
#define GL_ELEMENT_ARRAY_BARRIER_BIT 0x00000002
#define GL_UNIFORM_BARRIER_BIT 0x00000004
#define GL_TEXTURE_FETCH_BARRIER_BIT 0x00000008
#define GL_SHADER_IMAGE_ACCESS_BARRIER_BIT 0x00000020
#define GL_COMMAND_BARRIER_BIT 0x00000040
#define GL_PIXEL_BUFFER_BARRIER_BIT 0x00000080
#define GL_TEXTURE_UPDATE_BARRIER_BIT 0x00000100
#define GL_BUFFER_UPDATE_BARRIER_BIT 0x00000200
#define GL_FRAMEBUFFER_BARRIER_BIT 0x00000400
#define GL_TRANSFORM_FEEDBACK_BARRIER_BIT 0x00000800
#define GL_ATOMIC_COUNTER_BARRIER_BIT 0x00001000
#define GL_ALL_BARRIER_BITS 0xFFFFFFFF
#endif

#ifndef GL_SHADER_STORAGE_BARRIER_BIT
#define GL_SHADER_STORAGE_BARRIER_BIT 0x00002000
#endif

#ifndef GL_UNSIGNED_INT_ATOMIC_COUNTER
#define GL_UNSIGNED_INT_ATOMIC_COUNTER 0x92DB
#endif

#ifndef GL_UNSIGNED_INT_IMAGE_2D
#define GL_UNSIGNED_INT_IMAGE_2D 0x9063
#endif

#ifndef GL_SHADER_STORAGE_BLOCK
#define GL_SHADER_STORAGE_BLOCK 0x92E6
#endif

#ifndef GL_ACTIVE_RESOURCES
#define GL_ACTIVE_RESOURCES 0x92F5
#endif

#ifndef GL_BUFFER_BINDING
#define GL_ATOMIC_COUNTER_BUFFER_INDEX 0x9301
#define GL_BUFFER_BINDING 0x9302
#define GL_BUFFER_DATA_SIZE 0x9303
#define GL_NUM_ACTIVE_VARIABLES 0x9304
#define GL_ACTIVE_VARIABLES 0x9305
#endif

#ifndef GL_UNIFORM
#define GL_UNIFORM 0x92E1
#endif

#ifndef GL_COMPUTE_SHADER
#define GL_COMPUTE_SHADER 0x91B9
#endif

#ifndef GL_LUMINANCE16F_EXT
#define GL_LUMINANCE16F_EXT               0x881E
#endif

#ifndef GL_COMPRESSED_RGB_S3TC_DXT1_EXT
#define GL_COMPRESSED_RGB_S3TC_DXT1_EXT   0x83F0
#endif

#ifndef GL_COMPRESSED_RGBA_S3TC_DXT1_EXT
#define GL_COMPRESSED_RGBA_S3TC_DXT1_EXT  0x83F1
#endif

#ifndef GL_COMPRESSED_RGBA_S3TC_DXT3_EXT
#define GL_COMPRESSED_RGBA_S3TC_DXT3_EXT  0x83F2
#endif

#ifndef GL_COMPRESSED_RGBA_S3TC_DXT5_EXT
#define GL_COMPRESSED_RGBA_S3TC_DXT5_EXT  0x83F3
#endif

#ifndef GL_COMPRESSED_R11_EAC
#define GL_COMPRESSED_R11_EAC             0x9270
#endif

#ifndef GL_COMPRESSED_SIGNED_R11_EAC
#define GL_COMPRESSED_SIGNED_R11_EAC      0x9271
#endif

#ifndef GL_COMPRESSED_RG11_EAC
#define GL_COMPRESSED_RG11_EAC            0x9272
#endif

#ifndef GL_COMPRESSED_SIGNED_RG11_EAC
#define GL_COMPRESSED_SIGNED_RG11_EAC     0x9273
#endif

#ifndef GL_COMPRESSED_RGB8_ETC2
#define GL_COMPRESSED_RGB8_ETC2           0x9274
#endif

#ifndef GL_COMPRESSED_SRGB8_ETC2
#define GL_COMPRESSED_SRGB8_ETC2          0x9275
#endif

#ifndef GL_COMPRESSED_RGB8_PUNCHTHROUGH_ALPHA1_ETC2
#define GL_COMPRESSED_RGB8_PUNCHTHROUGH_ALPHA1_ETC2 0x9276
#endif

#ifndef GL_COMPRESSED_SRGB8_PUNCHTHROUGH_ALPHA1_ETC2
#define GL_COMPRESSED_SRGB8_PUNCHTHROUGH_ALPHA1_ETC2 0x9277
#endif

#ifndef GL_COMPRESSED_RGBA8_ETC2_EAC
#define GL_COMPRESSED_RGBA8_ETC2_EAC      0x9278
#endif

#ifndef GL_COMPRESSED_SRGB8_ALPHA8_ETC2_EAC
#define GL_COMPRESSED_SRGB8_ALPHA8_ETC2_EAC 0x9279
#endif

#ifndef GL_COMPRESSED_RGBA_ASTC_4x4_KHR
#define GL_COMPRESSED_RGBA_ASTC_4x4_KHR   0x93B0
#endif

#ifndef GL_COMPRESSED_RGBA_ASTC_5x4_KHR
#define GL_COMPRESSED_RGBA_ASTC_5x4_KHR   0x93B1
#endif

#ifndef GL_COMPRESSED_RGBA_ASTC_5x5_KHR
#define GL_COMPRESSED_RGBA_ASTC_5x5_KHR   0x93B2
#endif

#ifndef GL_COMPRESSED_RGBA_ASTC_6x5_KHR
#define GL_COMPRESSED_RGBA_ASTC_6x5_KHR   0x93B3
#endif

#ifndef GL_COMPRESSED_RGBA_ASTC_6x6_KHR
#define GL_COMPRESSED_RGBA_ASTC_6x6_KHR   0x93B4
#endif

#ifndef GL_COMPRESSED_RGBA_ASTC_8x5_KHR
#define GL_COMPRESSED_RGBA_ASTC_8x5_KHR   0x93B5
#endif

#ifndef GL_COMPRESSED_RGBA_ASTC_8x6_KHR
#define GL_COMPRESSED_RGBA_ASTC_8x6_KHR   0x93B6
#endif

#ifndef GL_COMPRESSED_RGBA_ASTC_8x8_KHR
#define GL_COMPRESSED_RGBA_ASTC_8x8_KHR   0x93B7
#endif

#ifndef GL_COMPRESSED_RGBA_ASTC_10x5_KHR
#define GL_COMPRESSED_RGBA_ASTC_10x5_KHR  0x93B8
#endif

#ifndef GL_COMPRESSED_RGBA_ASTC_10x6_KHR
#define GL_COMPRESSED_RGBA_ASTC_10x6_KHR  0x93B9
#endif

#ifndef GL_COMPRESSED_RGBA_ASTC_10x8_KHR
#define GL_COMPRESSED_RGBA_ASTC_10x8_KHR  0x93BA
#endif

#ifndef GL_COMPRESSED_RGBA_ASTC_10x10_KHR
#define GL_COMPRESSED_RGBA_ASTC_10x10_KHR 0x93BB
#endif

#ifndef GL_COMPRESSED_RGBA_ASTC_12x10_KHR
#define GL_COMPRESSED_RGBA_ASTC_12x10_KHR 0x93BC
#endif

#ifndef GL_COMPRESSED_RGBA_ASTC_12x12_KHR
#define GL_COMPRESSED_RGBA_ASTC_12x12_KHR 0x93BD
#endif

#ifndef GL_COMPRESSED_SRGB8_ALPHA8_ASTC_4x4_KHR
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_4x4_KHR 0x93D0
#endif

#ifndef GL_COMPRESSED_SRGB8_ALPHA8_ASTC_5x4_KHR
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_5x4_KHR 0x93D1
#endif

#ifndef GL_COMPRESSED_SRGB8_ALPHA8_ASTC_5x5_KHR
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_5x5_KHR 0x93D2
#endif

#ifndef GL_COMPRESSED_SRGB8_ALPHA8_ASTC_6x5_KHR
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_6x5_KHR 0x93D3
#endif

#ifndef GL_COMPRESSED_SRGB8_ALPHA8_ASTC_6x6_KHR
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_6x6_KHR 0x93D4
#endif

#ifndef GL_COMPRESSED_SRGB8_ALPHA8_ASTC_8x5_KHR
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_8x5_KHR 0x93D5
#endif

#ifndef GL_COMPRESSED_SRGB8_ALPHA8_ASTC_8x6_KHR
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_8x6_KHR 0x93D6
#endif

#ifndef GL_COMPRESSED_SRGB8_ALPHA8_ASTC_8x8_KHR
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_8x8_KHR 0x93D7
#endif

#ifndef GL_COMPRESSED_SRGB8_ALPHA8_ASTC_10x5_KHR
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_10x5_KHR 0x93D8
#endif

#ifndef GL_COMPRESSED_SRGB8_ALPHA8_ASTC_10x6_KHR
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_10x6_KHR 0x93D9
#endif

#ifndef GL_COMPRESSED_SRGB8_ALPHA8_ASTC_10x8_KHR
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_10x8_KHR 0x93DA
#endif

#ifndef GL_COMPRESSED_SRGB8_ALPHA8_ASTC_10x10_KHR
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_10x10_KHR 0x93DB
#endif

#ifndef GL_COMPRESSED_SRGB8_ALPHA8_ASTC_12x10_KHR
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_12x10_KHR 0x93DC
#endif

#ifndef GL_COMPRESSED_SRGB8_ALPHA8_ASTC_12x12_KHR
#define GL_COMPRESSED_SRGB8_ALPHA8_ASTC_12x12_KHR 0x93DD
#endif

#ifndef GL_RGB10_A2
#define GL_RGB10_A2                       0x8059
#endif

#ifndef GL_RGB16F
#define GL_RGB16F                         0x881B
#endif

#ifndef GL_RGBA32UI
#define GL_RGBA32UI                       0x8D70
#endif

#ifndef GL_RGB32UI
#define GL_RGB32UI                        0x8D71
#endif

#ifndef GL_RGBA16UI
#define GL_RGBA16UI                       0x8D76
#endif

#ifndef GL_RGB16UI
#define GL_RGB16UI                        0x8D77
#endif

#ifndef GL_RGBA8UI
#define GL_RGBA8UI                        0x8D7C
#endif

#ifndef GL_RGB8UI
#define GL_RGB8UI                         0x8D7D
#endif

#ifndef GL_RGBA32I
#define GL_RGBA32I                        0x8D82
#endif

#ifndef GL_RGB32I
#define GL_RGB32I                         0x8D83
#endif

#ifndef GL_RGBA16I
#define GL_RGBA16I                        0x8D88
#endif

#ifndef GL_RGB16I
#define GL_RGB16I                         0x8D89
#endif

#ifndef GL_RGBA8I
#define GL_RGBA8I                         0x8D8E
#endif

#ifndef GL_RGB8I
#define GL_RGB8I                          0x8D8F
#endif

#ifndef GL_RGB_INTEGER
#define GL_RGB_INTEGER                    0x8D98
#endif

#ifndef GL_RGBA_INTEGER
#define GL_RGBA_INTEGER                   0x8D99
#endif

#ifndef GL_TIMESTAMP
#define GL_TIMESTAMP                      0x8E28
#endif


QT_BEGIN_NAMESPACE

QT_END_NAMESPACE

#endif // QSSGOPENGLTOKENS_H
