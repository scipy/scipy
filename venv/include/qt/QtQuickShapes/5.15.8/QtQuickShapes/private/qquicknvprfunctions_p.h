/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtQuick module of the Qt Toolkit.
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

#ifndef QQUICKNVPRFUNCTIONS_P_H
#define QQUICKNVPRFUNCTIONS_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of a number of Qt sources files.  This header file may change from
// version to version without notice, or even be removed.
//
// We mean it.
//

#include <QtQuickShapes/private/qquickshapesglobal_p.h>
#include <qopengl.h>
#include <QtGui/qsurfaceformat.h>

#if QT_CONFIG(opengl)

QT_BEGIN_NAMESPACE

// note: fixed pipeline specific functions are removed - modern ES ext
// headers have all this, but not the fixed stuff

#ifndef GL_NV_path_rendering
#define GL_PATH_FORMAT_SVG_NV             0x9070
#define GL_PATH_FORMAT_PS_NV              0x9071
#define GL_STANDARD_FONT_NAME_NV          0x9072
#define GL_SYSTEM_FONT_NAME_NV            0x9073
#define GL_FILE_NAME_NV                   0x9074
#define GL_PATH_STROKE_WIDTH_NV           0x9075
#define GL_PATH_END_CAPS_NV               0x9076
#define GL_PATH_INITIAL_END_CAP_NV        0x9077
#define GL_PATH_TERMINAL_END_CAP_NV       0x9078
#define GL_PATH_JOIN_STYLE_NV             0x9079
#define GL_PATH_MITER_LIMIT_NV            0x907A
#define GL_PATH_DASH_CAPS_NV              0x907B
#define GL_PATH_INITIAL_DASH_CAP_NV       0x907C
#define GL_PATH_TERMINAL_DASH_CAP_NV      0x907D
#define GL_PATH_DASH_OFFSET_NV            0x907E
#define GL_PATH_CLIENT_LENGTH_NV          0x907F
#define GL_PATH_FILL_MODE_NV              0x9080
#define GL_PATH_FILL_MASK_NV              0x9081
#define GL_PATH_FILL_COVER_MODE_NV        0x9082
#define GL_PATH_STROKE_COVER_MODE_NV      0x9083
#define GL_PATH_STROKE_MASK_NV            0x9084
#define GL_COUNT_UP_NV                    0x9088
#define GL_COUNT_DOWN_NV                  0x9089
#define GL_PATH_OBJECT_BOUNDING_BOX_NV    0x908A
#define GL_CONVEX_HULL_NV                 0x908B
#define GL_BOUNDING_BOX_NV                0x908D
#define GL_TRANSLATE_X_NV                 0x908E
#define GL_TRANSLATE_Y_NV                 0x908F
#define GL_TRANSLATE_2D_NV                0x9090
#define GL_TRANSLATE_3D_NV                0x9091
#define GL_AFFINE_2D_NV                   0x9092
#define GL_AFFINE_3D_NV                   0x9094
#define GL_TRANSPOSE_AFFINE_2D_NV         0x9096
#define GL_TRANSPOSE_AFFINE_3D_NV         0x9098
#define GL_UTF8_NV                        0x909A
#define GL_UTF16_NV                       0x909B
#define GL_BOUNDING_BOX_OF_BOUNDING_BOXES_NV 0x909C
#define GL_PATH_COMMAND_COUNT_NV          0x909D
#define GL_PATH_COORD_COUNT_NV            0x909E
#define GL_PATH_DASH_ARRAY_COUNT_NV       0x909F
#define GL_PATH_COMPUTED_LENGTH_NV        0x90A0
#define GL_PATH_FILL_BOUNDING_BOX_NV      0x90A1
#define GL_PATH_STROKE_BOUNDING_BOX_NV    0x90A2
#define GL_SQUARE_NV                      0x90A3
#define GL_ROUND_NV                       0x90A4
#define GL_TRIANGULAR_NV                  0x90A5
#define GL_BEVEL_NV                       0x90A6
#define GL_MITER_REVERT_NV                0x90A7
#define GL_MITER_TRUNCATE_NV              0x90A8
#define GL_SKIP_MISSING_GLYPH_NV          0x90A9
#define GL_USE_MISSING_GLYPH_NV           0x90AA
#define GL_PATH_ERROR_POSITION_NV         0x90AB
#define GL_PATH_FOG_GEN_MODE_NV           0x90AC
#define GL_ACCUM_ADJACENT_PAIRS_NV        0x90AD
#define GL_ADJACENT_PAIRS_NV              0x90AE
#define GL_FIRST_TO_REST_NV               0x90AF
#define GL_PATH_GEN_MODE_NV               0x90B0
#define GL_PATH_GEN_COEFF_NV              0x90B1
#define GL_PATH_GEN_COLOR_FORMAT_NV       0x90B2
#define GL_PATH_GEN_COMPONENTS_NV         0x90B3
#define GL_PATH_STENCIL_FUNC_NV           0x90B7
#define GL_PATH_STENCIL_REF_NV            0x90B8
#define GL_PATH_STENCIL_VALUE_MASK_NV     0x90B9
#define GL_PATH_STENCIL_DEPTH_OFFSET_FACTOR_NV 0x90BD
#define GL_PATH_STENCIL_DEPTH_OFFSET_UNITS_NV 0x90BE
#define GL_PATH_COVER_DEPTH_FUNC_NV       0x90BF
#define GL_PATH_DASH_OFFSET_RESET_NV      0x90B4
#define GL_MOVE_TO_RESETS_NV              0x90B5
#define GL_MOVE_TO_CONTINUES_NV           0x90B6
#define GL_CLOSE_PATH_NV                  0x00
#define GL_MOVE_TO_NV                     0x02
#define GL_RELATIVE_MOVE_TO_NV            0x03
#define GL_LINE_TO_NV                     0x04
#define GL_RELATIVE_LINE_TO_NV            0x05
#define GL_HORIZONTAL_LINE_TO_NV          0x06
#define GL_RELATIVE_HORIZONTAL_LINE_TO_NV 0x07
#define GL_VERTICAL_LINE_TO_NV            0x08
#define GL_RELATIVE_VERTICAL_LINE_TO_NV   0x09
#define GL_QUADRATIC_CURVE_TO_NV          0x0A
#define GL_RELATIVE_QUADRATIC_CURVE_TO_NV 0x0B
#define GL_CUBIC_CURVE_TO_NV              0x0C
#define GL_RELATIVE_CUBIC_CURVE_TO_NV     0x0D
#define GL_SMOOTH_QUADRATIC_CURVE_TO_NV   0x0E
#define GL_RELATIVE_SMOOTH_QUADRATIC_CURVE_TO_NV 0x0F
#define GL_SMOOTH_CUBIC_CURVE_TO_NV       0x10
#define GL_RELATIVE_SMOOTH_CUBIC_CURVE_TO_NV 0x11
#define GL_SMALL_CCW_ARC_TO_NV            0x12
#define GL_RELATIVE_SMALL_CCW_ARC_TO_NV   0x13
#define GL_SMALL_CW_ARC_TO_NV             0x14
#define GL_RELATIVE_SMALL_CW_ARC_TO_NV    0x15
#define GL_LARGE_CCW_ARC_TO_NV            0x16
#define GL_RELATIVE_LARGE_CCW_ARC_TO_NV   0x17
#define GL_LARGE_CW_ARC_TO_NV             0x18
#define GL_RELATIVE_LARGE_CW_ARC_TO_NV    0x19
#define GL_RESTART_PATH_NV                0xF0
#define GL_DUP_FIRST_CUBIC_CURVE_TO_NV    0xF2
#define GL_DUP_LAST_CUBIC_CURVE_TO_NV     0xF4
#define GL_RECT_NV                        0xF6
#define GL_CIRCULAR_CCW_ARC_TO_NV         0xF8
#define GL_CIRCULAR_CW_ARC_TO_NV          0xFA
#define GL_CIRCULAR_TANGENT_ARC_TO_NV     0xFC
#define GL_ARC_TO_NV                      0xFE
#define GL_RELATIVE_ARC_TO_NV             0xFF
#define GL_BOLD_BIT_NV                    0x01
#define GL_ITALIC_BIT_NV                  0x02
#define GL_GLYPH_WIDTH_BIT_NV             0x01
#define GL_GLYPH_HEIGHT_BIT_NV            0x02
#define GL_GLYPH_HORIZONTAL_BEARING_X_BIT_NV 0x04
#define GL_GLYPH_HORIZONTAL_BEARING_Y_BIT_NV 0x08
#define GL_GLYPH_HORIZONTAL_BEARING_ADVANCE_BIT_NV 0x10
#define GL_GLYPH_VERTICAL_BEARING_X_BIT_NV 0x20
#define GL_GLYPH_VERTICAL_BEARING_Y_BIT_NV 0x40
#define GL_GLYPH_VERTICAL_BEARING_ADVANCE_BIT_NV 0x80
#define GL_GLYPH_HAS_KERNING_BIT_NV       0x100
#define GL_FONT_X_MIN_BOUNDS_BIT_NV       0x00010000
#define GL_FONT_Y_MIN_BOUNDS_BIT_NV       0x00020000
#define GL_FONT_X_MAX_BOUNDS_BIT_NV       0x00040000
#define GL_FONT_Y_MAX_BOUNDS_BIT_NV       0x00080000
#define GL_FONT_UNITS_PER_EM_BIT_NV       0x00100000
#define GL_FONT_ASCENDER_BIT_NV           0x00200000
#define GL_FONT_DESCENDER_BIT_NV          0x00400000
#define GL_FONT_HEIGHT_BIT_NV             0x00800000
#define GL_FONT_MAX_ADVANCE_WIDTH_BIT_NV  0x01000000
#define GL_FONT_MAX_ADVANCE_HEIGHT_BIT_NV 0x02000000
#define GL_FONT_UNDERLINE_POSITION_BIT_NV 0x04000000
#define GL_FONT_UNDERLINE_THICKNESS_BIT_NV 0x08000000
#define GL_FONT_HAS_KERNING_BIT_NV        0x10000000
#define GL_PRIMARY_COLOR_NV               0x852C
#define GL_SECONDARY_COLOR_NV             0x852D
#define GL_ROUNDED_RECT_NV                0xE8
#define GL_RELATIVE_ROUNDED_RECT_NV       0xE9
#define GL_ROUNDED_RECT2_NV               0xEA
#define GL_RELATIVE_ROUNDED_RECT2_NV      0xEB
#define GL_ROUNDED_RECT4_NV               0xEC
#define GL_RELATIVE_ROUNDED_RECT4_NV      0xED
#define GL_ROUNDED_RECT8_NV               0xEE
#define GL_RELATIVE_ROUNDED_RECT8_NV      0xEF
#define GL_RELATIVE_RECT_NV               0xF7
#define GL_FONT_GLYPHS_AVAILABLE_NV       0x9368
#define GL_FONT_TARGET_UNAVAILABLE_NV     0x9369
#define GL_FONT_UNAVAILABLE_NV            0x936A
#define GL_FONT_UNINTELLIGIBLE_NV         0x936B
#define GL_CONIC_CURVE_TO_NV              0x1A
#define GL_RELATIVE_CONIC_CURVE_TO_NV     0x1B
#define GL_FONT_NUM_GLYPH_INDICES_BIT_NV  0x20000000
#define GL_STANDARD_FONT_FORMAT_NV        0x936C
#define GL_2_BYTES_NV                     0x1407
#define GL_3_BYTES_NV                     0x1408
#define GL_4_BYTES_NV                     0x1409
#define GL_EYE_LINEAR_NV                  0x2400
#define GL_OBJECT_LINEAR_NV               0x2401
#define GL_CONSTANT_NV                    0x8576
#define GL_PATH_PROJECTION_NV             0x1701
#define GL_PATH_MODELVIEW_NV              0x1700
#define GL_PATH_MODELVIEW_STACK_DEPTH_NV  0x0BA3
#define GL_PATH_MODELVIEW_MATRIX_NV       0x0BA6
#define GL_PATH_MAX_MODELVIEW_STACK_DEPTH_NV 0x0D36
#define GL_PATH_TRANSPOSE_MODELVIEW_MATRIX_NV 0x84E3
#define GL_PATH_PROJECTION_STACK_DEPTH_NV 0x0BA4
#define GL_PATH_PROJECTION_MATRIX_NV      0x0BA7
#define GL_PATH_MAX_PROJECTION_STACK_DEPTH_NV 0x0D38
#define GL_PATH_TRANSPOSE_PROJECTION_MATRIX_NV 0x84E4
#define GL_FRAGMENT_INPUT_NV              0x936D

typedef GLuint (QOPENGLF_APIENTRYP PFNGLGENPATHSNVPROC) (GLsizei range);
typedef void (QOPENGLF_APIENTRYP PFNGLDELETEPATHSNVPROC) (GLuint path, GLsizei range);
typedef GLboolean (QOPENGLF_APIENTRYP PFNGLISPATHNVPROC) (GLuint path);
typedef void (QOPENGLF_APIENTRYP PFNGLPATHCOMMANDSNVPROC) (GLuint path, GLsizei numCommands, const GLubyte *commands, GLsizei numCoords, GLenum coordType, const void *coords);
typedef void (QOPENGLF_APIENTRYP PFNGLPATHCOORDSNVPROC) (GLuint path, GLsizei numCoords, GLenum coordType, const void *coords);
typedef void (QOPENGLF_APIENTRYP PFNGLPATHSUBCOMMANDSNVPROC) (GLuint path, GLsizei commandStart, GLsizei commandsToDelete, GLsizei numCommands, const GLubyte *commands, GLsizei numCoords, GLenum coordType, const void *coords);
typedef void (QOPENGLF_APIENTRYP PFNGLPATHSUBCOORDSNVPROC) (GLuint path, GLsizei coordStart, GLsizei numCoords, GLenum coordType, const void *coords);
typedef void (QOPENGLF_APIENTRYP PFNGLPATHSTRINGNVPROC) (GLuint path, GLenum format, GLsizei length, const void *pathString);
typedef void (QOPENGLF_APIENTRYP PFNGLPATHGLYPHSNVPROC) (GLuint firstPathName, GLenum fontTarget, const void *fontName, GLbitfield fontStyle, GLsizei numGlyphs, GLenum type, const void *charcodes, GLenum handleMissingGlyphs, GLuint pathParameterTemplate, GLfloat emScale);
typedef void (QOPENGLF_APIENTRYP PFNGLPATHGLYPHRANGENVPROC) (GLuint firstPathName, GLenum fontTarget, const void *fontName, GLbitfield fontStyle, GLuint firstGlyph, GLsizei numGlyphs, GLenum handleMissingGlyphs, GLuint pathParameterTemplate, GLfloat emScale);
typedef void (QOPENGLF_APIENTRYP PFNGLWEIGHTPATHSNVPROC) (GLuint resultPath, GLsizei numPaths, const GLuint *paths, const GLfloat *weights);
typedef void (QOPENGLF_APIENTRYP PFNGLCOPYPATHNVPROC) (GLuint resultPath, GLuint srcPath);
typedef void (QOPENGLF_APIENTRYP PFNGLINTERPOLATEPATHSNVPROC) (GLuint resultPath, GLuint pathA, GLuint pathB, GLfloat weight);
typedef void (QOPENGLF_APIENTRYP PFNGLTRANSFORMPATHNVPROC) (GLuint resultPath, GLuint srcPath, GLenum transformType, const GLfloat *transformValues);
typedef void (QOPENGLF_APIENTRYP PFNGLPATHPARAMETERIVNVPROC) (GLuint path, GLenum pname, const GLint *value);
typedef void (QOPENGLF_APIENTRYP PFNGLPATHPARAMETERINVPROC) (GLuint path, GLenum pname, GLint value);
typedef void (QOPENGLF_APIENTRYP PFNGLPATHPARAMETERFVNVPROC) (GLuint path, GLenum pname, const GLfloat *value);
typedef void (QOPENGLF_APIENTRYP PFNGLPATHPARAMETERFNVPROC) (GLuint path, GLenum pname, GLfloat value);
typedef void (QOPENGLF_APIENTRYP PFNGLPATHDASHARRAYNVPROC) (GLuint path, GLsizei dashCount, const GLfloat *dashArray);
typedef void (QOPENGLF_APIENTRYP PFNGLPATHSTENCILFUNCNVPROC) (GLenum func, GLint ref, GLuint mask);
typedef void (QOPENGLF_APIENTRYP PFNGLPATHSTENCILDEPTHOFFSETNVPROC) (GLfloat factor, GLfloat units);
typedef void (QOPENGLF_APIENTRYP PFNGLSTENCILFILLPATHNVPROC) (GLuint path, GLenum fillMode, GLuint mask);
typedef void (QOPENGLF_APIENTRYP PFNGLSTENCILSTROKEPATHNVPROC) (GLuint path, GLint reference, GLuint mask);
typedef void (QOPENGLF_APIENTRYP PFNGLSTENCILFILLPATHINSTANCEDNVPROC) (GLsizei numPaths, GLenum pathNameType, const void *paths, GLuint pathBase, GLenum fillMode, GLuint mask, GLenum transformType, const GLfloat *transformValues);
typedef void (QOPENGLF_APIENTRYP PFNGLSTENCILSTROKEPATHINSTANCEDNVPROC) (GLsizei numPaths, GLenum pathNameType, const void *paths, GLuint pathBase, GLint reference, GLuint mask, GLenum transformType, const GLfloat *transformValues);
typedef void (QOPENGLF_APIENTRYP PFNGLPATHCOVERDEPTHFUNCNVPROC) (GLenum func);
typedef void (QOPENGLF_APIENTRYP PFNGLCOVERFILLPATHNVPROC) (GLuint path, GLenum coverMode);
typedef void (QOPENGLF_APIENTRYP PFNGLCOVERSTROKEPATHNVPROC) (GLuint path, GLenum coverMode);
typedef void (QOPENGLF_APIENTRYP PFNGLCOVERFILLPATHINSTANCEDNVPROC) (GLsizei numPaths, GLenum pathNameType, const void *paths, GLuint pathBase, GLenum coverMode, GLenum transformType, const GLfloat *transformValues);
typedef void (QOPENGLF_APIENTRYP PFNGLCOVERSTROKEPATHINSTANCEDNVPROC) (GLsizei numPaths, GLenum pathNameType, const void *paths, GLuint pathBase, GLenum coverMode, GLenum transformType, const GLfloat *transformValues);
typedef void (QOPENGLF_APIENTRYP PFNGLGETPATHPARAMETERIVNVPROC) (GLuint path, GLenum pname, GLint *value);
typedef void (QOPENGLF_APIENTRYP PFNGLGETPATHPARAMETERFVNVPROC) (GLuint path, GLenum pname, GLfloat *value);
typedef void (QOPENGLF_APIENTRYP PFNGLGETPATHCOMMANDSNVPROC) (GLuint path, GLubyte *commands);
typedef void (QOPENGLF_APIENTRYP PFNGLGETPATHCOORDSNVPROC) (GLuint path, GLfloat *coords);
typedef void (QOPENGLF_APIENTRYP PFNGLGETPATHDASHARRAYNVPROC) (GLuint path, GLfloat *dashArray);
typedef void (QOPENGLF_APIENTRYP PFNGLGETPATHMETRICSNVPROC) (GLbitfield metricQueryMask, GLsizei numPaths, GLenum pathNameType, const void *paths, GLuint pathBase, GLsizei stride, GLfloat *metrics);
typedef void (QOPENGLF_APIENTRYP PFNGLGETPATHMETRICRANGENVPROC) (GLbitfield metricQueryMask, GLuint firstPathName, GLsizei numPaths, GLsizei stride, GLfloat *metrics);
typedef void (QOPENGLF_APIENTRYP PFNGLGETPATHSPACINGNVPROC) (GLenum pathListMode, GLsizei numPaths, GLenum pathNameType, const void *paths, GLuint pathBase, GLfloat advanceScale, GLfloat kerningScale, GLenum transformType, GLfloat *returnedSpacing);
typedef GLboolean (QOPENGLF_APIENTRYP PFNGLISPOINTINFILLPATHNVPROC) (GLuint path, GLuint mask, GLfloat x, GLfloat y);
typedef GLboolean (QOPENGLF_APIENTRYP PFNGLISPOINTINSTROKEPATHNVPROC) (GLuint path, GLfloat x, GLfloat y);
typedef GLfloat (QOPENGLF_APIENTRYP PFNGLGETPATHLENGTHNVPROC) (GLuint path, GLsizei startSegment, GLsizei numSegments);
typedef GLboolean (QOPENGLF_APIENTRYP PFNGLPOINTALONGPATHNVPROC) (GLuint path, GLsizei startSegment, GLsizei numSegments, GLfloat distance, GLfloat *x, GLfloat *y, GLfloat *tangentX, GLfloat *tangentY);
typedef void (QOPENGLF_APIENTRYP PFNGLMATRIXLOAD3X2FNVPROC) (GLenum matrixMode, const GLfloat *m);
typedef void (QOPENGLF_APIENTRYP PFNGLMATRIXLOAD3X3FNVPROC) (GLenum matrixMode, const GLfloat *m);
typedef void (QOPENGLF_APIENTRYP PFNGLMATRIXLOADTRANSPOSE3X3FNVPROC) (GLenum matrixMode, const GLfloat *m);
typedef void (QOPENGLF_APIENTRYP PFNGLMATRIXMULT3X2FNVPROC) (GLenum matrixMode, const GLfloat *m);
typedef void (QOPENGLF_APIENTRYP PFNGLMATRIXMULT3X3FNVPROC) (GLenum matrixMode, const GLfloat *m);
typedef void (QOPENGLF_APIENTRYP PFNGLMATRIXMULTTRANSPOSE3X3FNVPROC) (GLenum matrixMode, const GLfloat *m);
typedef void (QOPENGLF_APIENTRYP PFNGLSTENCILTHENCOVERFILLPATHNVPROC) (GLuint path, GLenum fillMode, GLuint mask, GLenum coverMode);
typedef void (QOPENGLF_APIENTRYP PFNGLSTENCILTHENCOVERSTROKEPATHNVPROC) (GLuint path, GLint reference, GLuint mask, GLenum coverMode);
typedef void (QOPENGLF_APIENTRYP PFNGLSTENCILTHENCOVERFILLPATHINSTANCEDNVPROC) (GLsizei numPaths, GLenum pathNameType, const void *paths, GLuint pathBase, GLenum fillMode, GLuint mask, GLenum coverMode, GLenum transformType, const GLfloat *transformValues);
typedef void (QOPENGLF_APIENTRYP PFNGLSTENCILTHENCOVERSTROKEPATHINSTANCEDNVPROC) (GLsizei numPaths, GLenum pathNameType, const void *paths, GLuint pathBase, GLint reference, GLuint mask, GLenum coverMode, GLenum transformType, const GLfloat *transformValues);
typedef GLenum (QOPENGLF_APIENTRYP PFNGLPATHGLYPHINDEXRANGENVPROC) (GLenum fontTarget, const void *fontName, GLbitfield fontStyle, GLuint pathParameterTemplate, GLfloat emScale, GLuint baseAndCount[2]);
typedef GLenum (QOPENGLF_APIENTRYP PFNGLPATHGLYPHINDEXARRAYNVPROC) (GLuint firstPathName, GLenum fontTarget, const void *fontName, GLbitfield fontStyle, GLuint firstGlyphIndex, GLsizei numGlyphs, GLuint pathParameterTemplate, GLfloat emScale);
typedef GLenum (QOPENGLF_APIENTRYP PFNGLPATHMEMORYGLYPHINDEXARRAYNVPROC) (GLuint firstPathName, GLenum fontTarget, GLsizeiptr fontSize, const void *fontData, GLsizei faceIndex, GLuint firstGlyphIndex, GLsizei numGlyphs, GLuint pathParameterTemplate, GLfloat emScale);
typedef void (QOPENGLF_APIENTRYP PFNGLPROGRAMPATHFRAGMENTINPUTGENNVPROC) (GLuint program, GLint location, GLenum genMode, GLint components, const GLfloat *coeffs);
typedef void (QOPENGLF_APIENTRYP PFNGLGETPROGRAMRESOURCEFVNVPROC) (GLuint program, GLenum programInterface, GLuint index, GLsizei propCount, const GLenum *props, GLsizei bufSize, GLsizei *length, GLfloat *params);
#endif

#ifndef GL_FLAT
#define GL_FLAT                           0x1D00
#endif

#ifndef GL_INVERT
#define GL_INVERT                         0x150A
#endif

// this one originates from fixed pipeline so may not be in GLES ext headers, but we need it still
#ifndef GL_OBJECT_LINEAR_NV
#define GL_OBJECT_LINEAR_NV               0x2401
#endif

#ifndef GL_EXT_direct_state_access
typedef void (QOPENGLF_APIENTRYP PFNGLMATRIXLOADFEXTPROC) (GLenum mode, const GLfloat *m);
typedef void (QOPENGLF_APIENTRYP PFNGLMATRIXLOADIDENTITYEXTPROC) (GLenum mode);
#endif

// When building on a system with GLES 2.0 or 3.0, we may still compile the NVPR
// code path even though it's never used. Keep it compiling by defining the
// necessary ES 3.1 separable program constants.
#ifndef GL_FRAGMENT_SHADER_BIT
#define GL_FRAGMENT_SHADER_BIT            0x00000002
#endif
#ifndef GL_UNIFORM
#define GL_UNIFORM                        0x92E1
#endif

class QQuickNvprFunctionsPrivate;

class QQuickNvprFunctions
{
public:
    QQuickNvprFunctions();
    ~QQuickNvprFunctions();

    static QSurfaceFormat format();
    static bool isSupported();

    bool create();

    bool createFragmentOnlyPipeline(const char *fragmentShaderSource, GLuint *pipeline, GLuint *program);

    PFNGLGENPATHSNVPROC genPaths = nullptr;
    PFNGLDELETEPATHSNVPROC deletePaths = nullptr;
    PFNGLISPATHNVPROC isPath = nullptr;
    PFNGLPATHCOMMANDSNVPROC pathCommands = nullptr;
    PFNGLPATHCOORDSNVPROC pathCoords = nullptr;
    PFNGLPATHSUBCOMMANDSNVPROC pathSubCommands = nullptr;
    PFNGLPATHSUBCOORDSNVPROC pathSubCoords = nullptr;
    PFNGLPATHSTRINGNVPROC pathString = nullptr;
    PFNGLPATHGLYPHSNVPROC pathGlyphs = nullptr;
    PFNGLPATHGLYPHRANGENVPROC pathGlyphRange = nullptr;
    PFNGLWEIGHTPATHSNVPROC weightPaths = nullptr;
    PFNGLCOPYPATHNVPROC copyPath = nullptr;
    PFNGLINTERPOLATEPATHSNVPROC interpolatePaths = nullptr;
    PFNGLTRANSFORMPATHNVPROC transformPath = nullptr;
    PFNGLPATHPARAMETERIVNVPROC pathParameteriv = nullptr;
    PFNGLPATHPARAMETERINVPROC pathParameteri = nullptr;
    PFNGLPATHPARAMETERFVNVPROC pathParameterfv = nullptr;
    PFNGLPATHPARAMETERFNVPROC pathParameterf = nullptr;
    PFNGLPATHDASHARRAYNVPROC pathDashArray = nullptr;
    PFNGLPATHSTENCILFUNCNVPROC pathStencilFunc = nullptr;
    PFNGLPATHSTENCILDEPTHOFFSETNVPROC pathStencilDepthOffset = nullptr;
    PFNGLSTENCILFILLPATHNVPROC stencilFillPath = nullptr;
    PFNGLSTENCILSTROKEPATHNVPROC stencilStrokePath = nullptr;
    PFNGLSTENCILFILLPATHINSTANCEDNVPROC stencilFillPathInstanced = nullptr;
    PFNGLSTENCILSTROKEPATHINSTANCEDNVPROC stencilStrokePathInstanced = nullptr;
    PFNGLPATHCOVERDEPTHFUNCNVPROC pathCoverDepthFunc = nullptr;
    PFNGLCOVERFILLPATHNVPROC coverFillPath = nullptr;
    PFNGLCOVERSTROKEPATHNVPROC coverStrokePath = nullptr;
    PFNGLCOVERFILLPATHINSTANCEDNVPROC coverFillPathInstanced = nullptr;
    PFNGLCOVERSTROKEPATHINSTANCEDNVPROC coverStrokePathInstanced = nullptr;
    PFNGLGETPATHPARAMETERIVNVPROC getPathParameteriv = nullptr;
    PFNGLGETPATHPARAMETERFVNVPROC getPathParameterfv = nullptr;
    PFNGLGETPATHCOMMANDSNVPROC getPathCommands = nullptr;
    PFNGLGETPATHCOORDSNVPROC getPathCoords = nullptr;
    PFNGLGETPATHDASHARRAYNVPROC getPathDashArray = nullptr;
    PFNGLGETPATHMETRICSNVPROC getPathMetrics = nullptr;
    PFNGLGETPATHMETRICRANGENVPROC getPathMetricRange = nullptr;
    PFNGLGETPATHSPACINGNVPROC getPathSpacing = nullptr;
    PFNGLISPOINTINFILLPATHNVPROC isPointInFillPath = nullptr;
    PFNGLISPOINTINSTROKEPATHNVPROC isPointInStrokePath = nullptr;
    PFNGLGETPATHLENGTHNVPROC getPathLength = nullptr;
    PFNGLPOINTALONGPATHNVPROC getPointAlongPath = nullptr;
    PFNGLMATRIXLOAD3X2FNVPROC matrixLoad3x2f = nullptr;
    PFNGLMATRIXLOAD3X3FNVPROC matrixLoad3x3f = nullptr;
    PFNGLMATRIXLOADTRANSPOSE3X3FNVPROC matrixLoadTranspose3x3f = nullptr;
    PFNGLMATRIXMULT3X2FNVPROC matrixMult3x2f = nullptr;
    PFNGLMATRIXMULT3X3FNVPROC matrixMult3x3f = nullptr;
    PFNGLMATRIXMULTTRANSPOSE3X3FNVPROC matrixMultTranspose3x3f = nullptr;
    PFNGLSTENCILTHENCOVERFILLPATHNVPROC stencilThenCoverFillPath = nullptr;
    PFNGLSTENCILTHENCOVERSTROKEPATHNVPROC stencilThenCoverStrokePath = nullptr;
    PFNGLSTENCILTHENCOVERFILLPATHINSTANCEDNVPROC stencilThenCoverFillPathInstanced = nullptr;
    PFNGLSTENCILTHENCOVERSTROKEPATHINSTANCEDNVPROC stencilThenCoverStrokePathInstanced = nullptr;
    PFNGLPATHGLYPHINDEXRANGENVPROC pathGlyphIndexRange = nullptr;
    PFNGLPATHGLYPHINDEXARRAYNVPROC pathGlyphIndexArray = nullptr;
    PFNGLPATHMEMORYGLYPHINDEXARRAYNVPROC pathMemoryGlyphIndexArray = nullptr;
    PFNGLPROGRAMPATHFRAGMENTINPUTGENNVPROC programPathFragmentInputGen = nullptr;
    PFNGLGETPROGRAMRESOURCEFVNVPROC getProgramResourcefv = nullptr;

    PFNGLMATRIXLOADFEXTPROC matrixLoadf = nullptr;
    PFNGLMATRIXLOADIDENTITYEXTPROC matrixLoadIdentity = nullptr;

private:
    QQuickNvprFunctionsPrivate *d;
};

QT_END_NAMESPACE

#endif // QT_CONFIG(opengl)

#endif // QQUICKNVPRFUNCTIONS_P_H
