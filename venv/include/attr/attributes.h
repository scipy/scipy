/*
 * Copyright (c) 2001-2002,2004 Silicon Graphics, Inc.
 * All Rights Reserved.
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef __ATTRIBUTES_H__
#define	__ATTRIBUTES_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

#ifndef ENOATTR
# define ENOATTR ENODATA
#endif

/*
 *	An almost-IRIX-compatible extended attributes API
 *	(the IRIX attribute "list" operation is missing, added ATTR_SECURE).
 */

/*
 * The maximum size (into the kernel or returned from the kernel) of an
 * attribute value or the buffer used for an attr_list() call.  Larger
 * sizes will result in an E2BIG return code.
 */
#define ATTR_MAX_VALUELEN	(64*1024)	/* max length of a value */


/*
 * Flags that can be used with any of the simple attribute calls.
 * All desired flags should be bit-wise OR'ed together.
 */
#define ATTR_DONTFOLLOW	0x0001	/* do not follow symlinks for a pathname */
#define ATTR_ROOT	0x0002	/* use root namespace attributes in op */
#define ATTR_TRUST	0x0004	/* tell server we can be trusted to properly
				   handle extended attributes */
#define ATTR_SECURE	0x0008	/* use security namespace attributes in op */

/*
 * Additional flags that can be used with the set() attribute call.
 * All desired flags (from both lists) should be bit-wise OR'ed together.
 */
#define ATTR_CREATE	0x0010	/* pure create: fail if attr already exists */
#define ATTR_REPLACE	0x0020	/* pure set: fail if attr does not exist */

/*
 * Define how lists of attribute names are returned to the user from
 * the attr_list() call.  A large, 32bit aligned, buffer is passed in
 * along with its size.  We put an array of offsets at the top that each
 * reference an attrlist_ent_t and pack the attrlist_ent_t's at the bottom.
 */
typedef struct attrlist {
	int32_t		al_count;	/* number of entries in attrlist */
	int32_t		al_more;	/* T/F: more attrs (do call again) */
	int32_t		al_offset[1];	/* byte offsets of attrs [var-sized] */
} attrlist_t;

/*
 * Show the interesting info about one attribute.  This is what the
 * al_offset[i] entry points to.
 */
typedef struct attrlist_ent {	/* data from attr_list() */
	uint32_t	a_valuelen;	/* number bytes in value of attr */
	char		a_name[1];	/* attr name (NULL terminated) */
} attrlist_ent_t;

/*
 * Given a pointer to the (char*) buffer containing the attr_list() result,
 * and an index, return a pointer to the indicated attribute in the buffer.
 */
#define	ATTR_ENTRY(buffer, index)		\
	((attrlist_ent_t *)			\
	 &((char *)buffer)[ ((attrlist_t *)(buffer))->al_offset[index] ])

/*
 * Implement a "cursor" for use in successive attr_list() calls.
 * It provides a way to find the last attribute that was returned in the
 * last attr_list() call so that we can get the next one without missing
 * any.  This should be zeroed before use and whenever it is desired to
 * start over from the beginning of the attribute list.  The only valid
 * operation on a cursor is to zero it.
 */
typedef struct attrlist_cursor {
	uint32_t	opaque[4];	/* an opaque cookie */
} attrlist_cursor_t;

/*
 * Multi-attribute operation vector.
 */
typedef struct attr_multiop {
	int32_t	am_opcode;	/* operation to perform (ATTR_OP_GET, etc.) */
	int32_t	am_error;	/* [out arg] result of this sub-op (an errno) */
	char	*am_attrname;	/* attribute name to work with */
	char	*am_attrvalue;	/* [in/out arg] attribute value (raw bytes) */
	int32_t	am_length;	/* [in/out arg] length of value */
	int32_t	am_flags;	/* flags (bit-wise OR of #defines above) */
} attr_multiop_t;
#define	ATTR_MAX_MULTIOPS	128	/* max number ops in an oplist array */

/*
 * Valid values of am_opcode.
 */
#define ATTR_OP_GET	1	/* return the indicated attr's value */
#define ATTR_OP_SET	2	/* set/create the indicated attr/value pair */
#define ATTR_OP_REMOVE	3	/* remove the indicated attr */

/*
 * Get the value of an attribute.
 * Valuelength must be set to the maximum size of the value buffer, it will
 * be set to the actual number of bytes used in the value buffer upon return.
 * The return value is -1 on error (w/errno set appropriately), 0 on success.
 */
extern int attr_get (const char *__path, const char *__attrname,
			char *__attrvalue, int *__valuelength, int __flags)
	__attribute__ ((deprecated ("Use getxattr or lgetxattr instead")));
extern int attr_getf (int __fd, const char *__attrname, char *__attrvalue,
			int *__valuelength, int __flags)
	__attribute__ ((deprecated ("Use fgetxattr instead")));

/*
 * Set the value of an attribute, creating the attribute if necessary.
 * The return value is -1 on error (w/errno set appropriately), 0 on success.
 */
extern int attr_set (const char *__path, const char *__attrname,
			const char *__attrvalue, const int __valuelength,
			int __flags)
	__attribute__ ((deprecated ("Use setxattr or lsetxattr instead")));
extern int attr_setf (int __fd, const char *__attrname,
			const char *__attrvalue, const int __valuelength,
			int __flags)
	__attribute__ ((deprecated ("Use fsetxattr instead")));

/*
 * Remove an attribute.
 * The return value is -1 on error (w/errno set appropriately), 0 on success.
 */
extern int attr_remove (const char *__path, const char *__attrname,
			int __flags)
	__attribute__ ((deprecated ("Use removexattr or lremovexattr instead")));
extern int attr_removef (int __fd, const char *__attrname, int __flags)
	__attribute__ ((deprecated ("Use fremovexattr instead")));

/*
 * List the names and sizes of the values of all the attributes of an object.
 * "Cursor" must be allocated and zeroed before the first call, it is used
 * to maintain context between system calls if all the attribute names won't
 * fit into the buffer on the first system call.
 * The return value is -1 on error (w/errno set appropriately), 0 on success.
 */
extern int attr_list(const char *__path, char *__buffer, const int __buffersize,
		int __flags, attrlist_cursor_t *__cursor)
	__attribute__ ((deprecated ("Use listxattr or llistxattr instead")));
extern int attr_listf(int __fd, char *__buffer, const int __buffersize,
		int __flags, attrlist_cursor_t *__cursor)
	__attribute__ ((deprecated ("Use flistxattr instead")));

/*
 * Operate on multiple attributes of the same object simultaneously.
 *
 * This call will save on system call overhead when many attributes are
 * going to be operated on.
 *
 * The return value is -1 on error (w/errno set appropriately), 0 on success.
 * Note that this call will not return -1 as a result of failure of any
 * of the sub-operations, their return value is stored in each element
 * of the operation array.  This call will return -1 for a failure of the
 * call as a whole, eg: if the pathname doesn't exist, or the fd is bad.
 *
 * The semantics and allowable values for the fields in a attr_multiop_t
 * are the same as the semantics and allowable values for the arguments to
 * the corresponding "simple" attribute interface.  For example: the args
 * to a ATTR_OP_GET are the same as the args to an attr_get() call.
 */
extern int attr_multi (const char *__path, attr_multiop_t *__oplist,
			int __count, int __flags)
	__attribute__ ((deprecated ("Use getxattr, setxattr, listxattr, removexattr instead")));
extern int attr_multif (int __fd, attr_multiop_t *__oplist,
			int __count, int __flags)
	__attribute__ ((deprecated ("Use getxattr, setxattr, listxattr, removexattr instead")));

#ifdef __cplusplus
}
#endif

#endif	/* __ATTRIBUTES_H__ */
