/****************************************************************************
**
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

#ifndef QSSGMESHUTILITIES_P_H
#define QSSGMESHUTILITIES_P_H

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

#include <QtQuick3DAssetImport/private/qtquick3dassetimportglobal_p.h>

#include <QtQuick3DUtils/private/qssgbounds3_p.h>

#include <QtQuick3DRender/private/qssgrenderbasetypes_p.h>

#include <QtCore/QString>
#include <QtCore/QByteArray>
#include <QtCore/QIODevice>
#include <QtCore/QFile>

QT_BEGIN_NAMESPACE

namespace QSSGMeshUtilities {

struct MeshData
{
    // All enums must match the ones defined by QSSGRenderGeometry class in quick3d module
    enum PrimitiveType { // Must match also internal QSSGRenderDrawMode
        UnknownType = 0,
        Points,
        LineStrip,
        LineLoop,
        Lines,
        TriangleStrip,
        TriangleFan,
        Triangles, // Default primitive type
        Patches
    };

    struct Attribute {
        enum Semantic {
            UnknownSemantic = 0,
            IndexSemantic,
            PositionSemantic, // attr_pos
            NormalSemantic,   // attr_norm
            TexCoordSemantic, // attr_uv0
            TangentSemantic,  // attr_textan
            BinormalSemantic  // attr_binormal
        };
        enum ComponentType { // Must match also internal QSSGRenderComponentType
            DefaultType = 0,
            U8Type,
            I8Type,
            U16Type,
            I16Type,
            U32Type, // Default for IndexSemantic
            I32Type,
            U64Type,
            I64Type,
            F16Type,
            F32Type, // Default for other semantics
            F64Type
        };

        int typeSize() const
        {
            switch (componentType) {
            case U8Type:  return 1;
            case I8Type:  return 1;
            case U16Type: return 2;
            case I16Type: return 2;
            case U32Type: return 4;
            case I32Type: return 4;
            case U64Type: return 8;
            case I64Type: return 8;
            case F16Type: return 2;
            case F32Type: return 4;
            case F64Type: return 8;
            default:
                Q_ASSERT(false);
                return 0;
            }
        }

        int componentCount() const
        {
            switch (semantic) {
            case IndexSemantic:    return 1;
            case PositionSemantic: return 3;
            case NormalSemantic:   return 3;
            case TexCoordSemantic: return 2;
            case TangentSemantic:  return 3;
            case BinormalSemantic: return 3;
            default:
                Q_ASSERT(false);
                return 0;
            }
        }

        Semantic semantic = PositionSemantic;
        ComponentType componentType = F32Type;
        int offset = 0;
    };

    static const int MAX_ATTRIBUTES = 6;

    void clear()
    {
        m_vertexBuffer.clear();
        m_indexBuffer.clear();
        m_attributeCount = 0;
        m_primitiveType = Triangles;
    }

    QByteArray m_vertexBuffer;
    QByteArray m_indexBuffer;

    Attribute m_attributes[MAX_ATTRIBUTES];
    int m_attributeCount = 0;
    PrimitiveType m_primitiveType = Triangles;
    int m_stride = 0;
};

template<typename DataType>
struct OffsetDataRef
{
    quint32 m_offset;
    quint32 m_size;
    OffsetDataRef() : m_offset(0), m_size(0) {}
    DataType *begin(quint8 *inBase)
    {
        DataType *value = reinterpret_cast<DataType *>(inBase + m_offset);
        return value;
    }
    DataType *end(quint8 *inBase) { return begin(inBase) + m_size; }
    const DataType *begin(const quint8 *inBase) const { return reinterpret_cast<const DataType *>(inBase + m_offset); }
    const DataType *end(const quint8 *inBase) const { return begin(inBase) + m_size; }
    quint32 size() const { return m_size; }
    bool empty() const { return m_size == 0; }
    DataType &index(quint8 *inBase, quint32 idx)
    {
        Q_ASSERT(idx < m_size);
        return begin(inBase)[idx];
    }
    const DataType &index(const quint8 *inBase, quint32 idx) const
    {
        Q_ASSERT(idx < m_size);
        return begin(inBase)[idx];
    }
};

struct MeshVertexBufferEntry
{
    quint32 m_nameOffset;
    /** Datatype of the this entry points to in the buffer */
    QSSGRenderComponentType m_componentType;
    /** Number of components of each data member. 1,2,3, or 4.  Don't be stupid.*/
    quint32 m_numComponents;
    /** Offset from the beginning of the buffer of the first item */
    quint32 m_firstItemOffset;
    MeshVertexBufferEntry()
        : m_nameOffset(0), m_componentType(QSSGRenderComponentType::Float32), m_numComponents(3), m_firstItemOffset(0)
    {
    }
    QSSGRenderVertexBufferEntry toVertexBufferEntry(quint8 *inBaseAddress) const
    {
        const char *nameBuffer = "";
        if (m_nameOffset)
            nameBuffer = reinterpret_cast<const char *>(inBaseAddress + m_nameOffset);
        return QSSGRenderVertexBufferEntry(nameBuffer, m_componentType, m_numComponents, m_firstItemOffset);
    }
};

struct VertexBuffer
{
    OffsetDataRef<MeshVertexBufferEntry> m_entries;
    quint32 m_stride;
    OffsetDataRef<quint8> m_data;
    VertexBuffer(OffsetDataRef<MeshVertexBufferEntry> entries, quint32 stride, OffsetDataRef<quint8> data)
        : m_entries(entries), m_stride(stride), m_data(data)
    {
    }
    VertexBuffer() : m_stride(0) {}
};

struct IndexBuffer
{
    // Component types must be either UnsignedInt16 or UnsignedInt8 in order for the
    // graphics hardware to deal with the buffer correctly.
    QSSGRenderComponentType m_componentType;
    OffsetDataRef<quint8> m_data;
    // Either quint8 or quint16 component types are allowed by the underlying rendering
    // system, so you would be wise to stick with those.
    IndexBuffer(QSSGRenderComponentType compType, OffsetDataRef<quint8> data)
        : m_componentType(compType), m_data(data)
    {
    }
    IndexBuffer() : m_componentType(QSSGRenderComponentType::Unknown) {}
};

template<quint32 TNumBytes>
struct MeshPadding
{
    quint8 m_padding[TNumBytes];
    MeshPadding() { memZero(m_padding, TNumBytes); }
};

struct Vec3
{
    float x;
    float y;
    float z;
};

struct MeshSubset
{
    // std::numeric_limits<quint32>::max() means use all available items
    quint32 m_count;
    // Offset is in item size, not bytes.
    quint32 m_offset;
    // Bounds of this subset.  This is filled in by the builder
    // see AddMeshSubset
    QSSGBounds3 m_bounds;

    // Subsets have to be named else artists will be unable to use
    // a mesh with multiple subsets as they won't have any idea
    // while part of the model a given mesh actually maps to.
    OffsetDataRef<char16_t> m_name;

    MeshSubset(quint32 count, quint32 off, const QSSGBounds3 &bounds, OffsetDataRef<char16_t> inName)
        : m_count(count), m_offset(off), m_bounds(bounds), m_name(inName)
    {
    }
    MeshSubset() : m_count(quint32(-1)), m_offset(0), m_bounds() {}
    bool hasCount() const { return m_count != 4294967295U; } // AKA U_MAX 0xffffffff
};

struct Joint
{
    qint32 m_jointID;
    qint32 m_parentID;
    float m_invBindPose[16];
    float m_localToGlobalBoneSpace[16];

    Joint(qint32 jointID, qint32 parentID, const float *invBindPose, const float *localToGlobalBoneSpace)
        : m_jointID(jointID), m_parentID(parentID)
    {
        ::memcpy(m_invBindPose, invBindPose, sizeof(float) * 16);
        ::memcpy(m_localToGlobalBoneSpace, localToGlobalBoneSpace, sizeof(float) * 16);
    }
    Joint() : m_jointID(-1), m_parentID(-1)
    {
        ::memset(m_invBindPose, 0, sizeof(float) * 16);
        ::memset(m_localToGlobalBoneSpace, 0, sizeof(float) * 16);
    }
};

// Tells us what offset a mesh with this ID starts.
struct MeshMultiEntry
{
    quint64 m_meshOffset;
    quint32 m_meshId;
    quint32 m_padding;
    MeshMultiEntry() : m_meshOffset(0), m_meshId(0), m_padding(0) {}
    MeshMultiEntry(quint64 mo, quint32 meshId) : m_meshOffset(mo), m_meshId(meshId), m_padding(0) {}
};

// The multi headers are actually saved at the end of the file.
// Thus when you append to the file we overwrite the last header
// then write out a new header structure.
// The last 8 bytes of the file contain the multi header.
// The previous N*8 bytes contain the mesh entries.
struct MeshMultiHeader
{
    quint32 m_fileId;
    quint32 m_version;
    OffsetDataRef<MeshMultiEntry> m_entries;
    static quint32 getMultiStaticFileId() { return 555777497U; }
    static quint32 getMultiStaticVersion() { return 1; }

    MeshMultiHeader() : m_fileId(getMultiStaticFileId()), m_version(getMultiStaticVersion()) {}
};

struct Mesh;

// Result of a multi-load operation.  This returns both the mesh
// and the id of the mesh that was loaded.
struct MultiLoadResult
{
    Mesh *m_mesh;
    quint32 m_id;
    MultiLoadResult(Mesh *inMesh, quint32 inId) : m_mesh(inMesh), m_id(inId) {}
    MultiLoadResult() : m_mesh(nullptr), m_id(0) {}
    operator Mesh *() { return m_mesh; }
};

struct Q_QUICK3DASSETIMPORT_EXPORT Mesh
{
    static const char16_t *m_defaultName;

    VertexBuffer m_vertexBuffer;
    IndexBuffer m_indexBuffer;
    OffsetDataRef<MeshSubset> m_subsets;
    OffsetDataRef<Joint> m_joints;
    QSSGRenderDrawMode m_drawMode;
    QSSGRenderWinding m_winding;

    Mesh() : m_drawMode(QSSGRenderDrawMode::Triangles), m_winding(QSSGRenderWinding::CounterClockwise) {}
    Mesh(VertexBuffer vbuf,
         IndexBuffer ibuf,
         const OffsetDataRef<MeshSubset> &insts,
         const OffsetDataRef<Joint> &joints,
         QSSGRenderDrawMode drawMode = QSSGRenderDrawMode::Triangles,
         QSSGRenderWinding winding = QSSGRenderWinding::CounterClockwise)
        : m_vertexBuffer(vbuf), m_indexBuffer(ibuf), m_subsets(insts), m_joints(joints), m_drawMode(drawMode), m_winding(winding)
    {
    }

    quint8 *getBaseAddress() { return reinterpret_cast<quint8 *>(this); }
    const quint8 *getBaseAddress() const { return reinterpret_cast<const quint8 *>(this); }

    static const char *getPositionAttrName() { return "attr_pos"; }
    static const char *getNormalAttrName() { return "attr_norm"; }
    static const char *getUVAttrName() { return "attr_uv0"; }
    static const char *getUV2AttrName() { return "attr_uv1"; }
    static const char *getTexTanAttrName() { return "attr_textan"; }
    static const char *getTexBinormalAttrName() { return "attr_binormal"; }
    static const char *getWeightAttrName() { return "attr_weight"; }
    static const char *getBoneIndexAttrName() { return "attr_boneid"; }
    static const char *getColorAttrName() { return "attr_color"; }

    // Run through the vertex buffer items indicated by subset
    // Assume vbuf entry[posEntryIndex] is the position entry
    // This entry has to be QT3DSF32 and 3 components.
    // Using this entry and the (possibly empty) index buffer
    // along with the (possibly emtpy) logical vbuf data
    // return a bounds of the given vertex buffer.
    static QSSGBounds3 calculateSubsetBounds(const QSSGRenderVertexBufferEntry &inEntry,
                                               const QByteArray &inVertxData,
                                               quint32 inStride,
                                               const QByteArray &inIndexData,
                                               QSSGRenderComponentType inIndexCompType,
                                               quint32 inSubsetCount,
                                               quint32 inSubsetOffset);

    // Format is:
    // MeshDataHeader
    // mesh data.

    void save(QIODevice &outStream) const;

    // Save a mesh using fopen and fwrite
    bool save(const char *inFilePath) const;

    // read the header, then read the object.
    // Load a mesh using fopen and fread
    // Mesh needs to be freed by the caller using free
    static Mesh *load(QIODevice &inStream);
    static Mesh *load(const char *inFilePath);

    // Create a mesh given this header, and that data.  data.size() must match
    // header.SizeInBytes.  The mesh returned starts a data[0], so however data
    // was allocated is how the mesh should be deallocated.
    static Mesh *initialize(quint16 meshVersion, quint16 meshFlags, QSSGByteView data);

    // You can save multiple meshes in a file.  Each mesh returns an incrementing
    // integer for the multi file.  The original meshes aren't changed, and the file
    // is appended to.
    quint32 saveMulti(QIODevice &inStream, quint32 inId = 0) const;
    quint32 saveMulti(const char *inFilePath) const;

    // Load a single mesh using c file API and malloc/free.
    static MultiLoadResult loadMulti(QIODevice &inStream, quint32 inId);
    static MultiLoadResult loadMulti(const char *inFilePath, quint32 inId);

    // Returns true if this is a multimesh (several meshes in one file).
    static bool isMulti(QIODevice &inStream);

    // Load a multi header from a file using malloc.  Header needs to be freed using free.
    static MeshMultiHeader *loadMultiHeader(QIODevice &inStream);
    static MeshMultiHeader *loadMultiHeader(const char *inFilePath);

    // Get the highest mesh version from a file.
    static quint32 getHighestMultiVersion(QIODevice &inStream);
    static quint32 getHighestMultiVersion(const char *inFilePath);
};

struct MeshDataHeader
{
    static quint32 getFileId() { return quint32(-929005747); }
    static quint16 getCurrentFileVersion() { return 3; }
    quint32 m_fileId;
    quint16 m_fileVersion;
    quint16 m_headerFlags;
    quint32 m_sizeInBytes;
    MeshDataHeader(quint32 size = 0)
        : m_fileId(getFileId()), m_fileVersion(getCurrentFileVersion()), m_sizeInBytes(size)
    {
    }
};

struct MeshBuilderVBufEntry
{
    const char *m_name;
    QByteArray m_data;
    QSSGRenderComponentType m_componentType;
    quint32 m_numComponents;
    MeshBuilderVBufEntry() : m_name(nullptr), m_componentType(QSSGRenderComponentType::Unknown), m_numComponents(0)
    {
    }
    MeshBuilderVBufEntry(const char *name, const QByteArray &data, QSSGRenderComponentType componentType, quint32 numComponents)
        : m_name(name), m_data(data), m_componentType(componentType), m_numComponents(numComponents)
    {
    }
};

// Useful class to build up a mesh.  Necessary since meshes don't include that
// sort of utility.
class Q_QUICK3DASSETIMPORT_EXPORT QSSGMeshBuilder
{
public:
    QAtomicInt ref;
    virtual ~QSSGMeshBuilder();
    virtual void release() = 0;
    virtual void reset() = 0;
    // Set the draw parameters for any subsets.  Defaults to triangles and counter clockwise
    virtual void setDrawParameters(QSSGRenderDrawMode drawMode, QSSGRenderWinding winding) = 0;
    // Set the vertex buffer and have the mesh builder interleave the data for you
    virtual bool setVertexBuffer(const QVector<MeshBuilderVBufEntry> &entries) = 0;
    // Set the vertex buffer from interleaved data.
    virtual void setVertexBuffer(const QVector<QSSGRenderVertexBufferEntry> &entries, quint32 stride, QByteArray data) = 0;
    // The builder (and the majority of the rest of the product) only supports unsigned 16 bit
    // indexes
    virtual void setIndexBuffer(const QByteArray &data, QSSGRenderComponentType comp) = 0;
    // Assets if the supplied parameters are out of range.
    virtual void addJoint(qint32 jointID, qint32 parentID, const float *invBindPose, const float *localToGlobalBoneSpace) = 0;
    /**
     *  Add a subset, which equates roughly to a draw call.
     *  A logical vertex buffer allows you to have more that 64K vertexes but still
     *  use u16 index buffers.  In any case, if the mesh has an index buffer then this subset
     *  refers to that index buffer, else it is assumed to index into the vertex buffer.
     *  count and offset do exactly what they seem to do, while boundsPositionEntryIndex,
     *  if set to something other than std::numeric_limits<quint32>::max(),
     *  drives the calculation of the aa-bounds of the subset using mesh::CalculateSubsetBounds.
     */
    virtual void addMeshSubset(const char16_t *inSubsetName = Mesh::m_defaultName,
                               quint32 count = std::numeric_limits<quint32>::max(),
                               quint32 offset = 0,
                               quint32 boundsPositionEntryIndex = std::numeric_limits<quint32>::max()) = 0;

    virtual void addMeshSubset(const char16_t *inSubsetName, quint32 count, quint32 offset, const QSSGBounds3 &inBounds) = 0;

    // Call to optimize the index and vertex buffers.  This doesn't change the subset information,
    // each triangle is rendered precisely the same.
    // It just orders the vertex data so we iterate through it as linearly as possible.
    // This *only* works if the *entire* builder is using triangles as the draw mode.  This will be
    // a disaster if that condition is not met.
    virtual void optimizeMesh() = 0;

    /**
     * @brief This functions stitches together sub-meshes with the same material.
     *		 This re-writes the index buffer
     *
     * @return no return.
     */
    virtual void connectSubMeshes() = 0;

    // Return the current mesh.  This is only good for this function call, item may change or be
    // released
    // due to any further function calls.
    virtual Mesh &getMesh() = 0;

    virtual Mesh *buildMesh(const MeshData &data, QString &error, const QSSGBounds3 &inBounds) = 0;

    // Uses new/delete.
    static QSSGRef<QSSGMeshBuilder> createMeshBuilder();
};

} // end QSSGMeshUtilities namespace

QT_END_NAMESPACE

#endif // QSSGMESHUTILITIES_P_H
