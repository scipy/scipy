/* Public API for GNU gettext PO files - contained in libgettextpo.
   Copyright (C) 2003-2008, 2010, 2012-2016, 2019-2022 Free Software Foundation, Inc.
   Written by Bruno Haible <bruno@clisp.org>, 2003.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.  */

#ifndef _GETTEXT_PO_H
#define _GETTEXT_PO_H 1

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif


/* =========================== Meta Information ============================ */

/* Version number: (major<<16) + (minor<<8) + subminor */
#define LIBGETTEXTPO_VERSION 0x001501
extern int libgettextpo_version;

/* ================================= Types ================================= */

/* A po_file_t represents the contents of a PO file.  */
typedef struct po_file *po_file_t;

/* A po_message_iterator_t represents an iterator through a domain of a
   PO file.  */
typedef struct po_message_iterator *po_message_iterator_t;

/* A po_message_t represents a message in a PO file.  */
typedef struct po_message *po_message_t;

/* A po_filepos_t represents a string's position within a source file.  */
typedef struct po_filepos *po_filepos_t;

/* A po_error_handler handles error situations.  */
struct po_error_handler
{
  /* Signal an error.  The error message is built from FORMAT and the following
     arguments.  ERRNUM, if nonzero, is an errno value.
     Must increment the error_message_count variable declared in error.h.
     Must not return if STATUS is nonzero.  */
  void (*error) (int status, int errnum,
                 const char *format, ...)
#if ((__GNUC__ == 3 && __GNUC_MINOR__ >= 1) || __GNUC__ > 3) && !__STRICT_ANSI__
  __attribute__ ((__format__ (__printf__, 3, 4)))
#endif
  ;

  /* Signal an error.  The error message is built from FORMAT and the following
     arguments.  The error location is at FILENAME line LINENO. ERRNUM, if
     nonzero, is an errno value.
     Must increment the error_message_count variable declared in error.h.
     Must not return if STATUS is nonzero.  */
  void (*error_at_line) (int status, int errnum,
                         const char *filename, unsigned int lineno,
                         const char *format, ...)
#if ((__GNUC__ == 3 && __GNUC_MINOR__ >= 1) || __GNUC__ > 3) && !__STRICT_ANSI__
  __attribute__ ((__format__ (__printf__, 5, 6)))
#endif
  ;

  /* Signal a multiline warning.  The PREFIX applies to all lines of the
     MESSAGE.  Free the PREFIX and MESSAGE when done.  */
  void (*multiline_warning) (char *prefix, char *message);

  /* Signal a multiline error.  The PREFIX applies to all lines of the
     MESSAGE.  Free the PREFIX and MESSAGE when done.
     Must increment the error_message_count variable declared in error.h if
     PREFIX is non-NULL.  */
  void (*multiline_error) (char *prefix, char *message);
};
typedef const struct po_error_handler *po_error_handler_t;

/* A po_xerror_handler handles warnings, error and fatal error situations.  */
#define PO_SEVERITY_WARNING     0 /* just a warning, tell the user */
#define PO_SEVERITY_ERROR       1 /* an error, the operation cannot complete */
#define PO_SEVERITY_FATAL_ERROR 2 /* an error, the operation must be aborted */
struct po_xerror_handler
{
  /* Signal a problem of the given severity.
     MESSAGE and/or FILENAME + LINENO indicate where the problem occurred.
     If FILENAME is NULL, FILENAME and LINENO and COLUMN should be ignored.
     If LINENO is (size_t)(-1), LINENO and COLUMN should be ignored.
     If COLUMN is (size_t)(-1), it should be ignored.
     MESSAGE_TEXT is the problem description (if MULTILINE_P is true,
     multiple lines of text, each terminated with a newline, otherwise
     usually a single line).
     Must not return if SEVERITY is PO_SEVERITY_FATAL_ERROR.  */
  void (*xerror) (int severity,
                  po_message_t message,
                  const char *filename, size_t lineno, size_t column,
                  int multiline_p, const char *message_text);
  /* Signal a problem that refers to two messages.
     Similar to two calls to xerror.
     If possible, a "..." can be appended to MESSAGE_TEXT1 and prepended to
     MESSAGE_TEXT2.  */
  void (*xerror2) (int severity,
                   po_message_t message1,
                   const char *filename1, size_t lineno1, size_t column1,
                   int multiline_p1, const char *message_text1,
                   po_message_t message2,
                   const char *filename2, size_t lineno2, size_t column2,
                   int multiline_p2, const char *message_text2);
};
typedef const struct po_xerror_handler *po_xerror_handler_t;

/* Memory allocation:
   The memory allocations performed by these functions use xmalloc(),
   therefore will cause a program exit if memory is exhausted.
   The memory allocated by po_file_read, and implicitly returned through
   the po_message_* functions, lasts until freed with po_file_free.  */


/* ============================= po_file_t API ============================= */

/* Create an empty PO file representation in memory.  */
extern po_file_t po_file_create (void);

/* Read a PO file into memory.
   Return its contents.  Upon failure, call function from handler.  */
#define po_file_read po_file_read_v3
extern po_file_t po_file_read (const char *filename,
                               po_xerror_handler_t handler);

/* Write an in-memory PO file to a file.
   Upon failure, call function from handler.  */
#define po_file_write po_file_write_v2
extern po_file_t po_file_write (po_file_t file, const char *filename,
                                po_xerror_handler_t handler);

/* Free a PO file from memory.  */
extern void po_file_free (po_file_t file);

/* Return the names of the domains covered by a PO file in memory.  */
extern const char * const * po_file_domains (po_file_t file);


/* =========================== Header entry API ============================ */

/* Return the header entry of a domain of a PO file in memory.
   The domain NULL denotes the default domain.
   Return NULL if there is no header entry.  */
extern const char * po_file_domain_header (po_file_t file, const char *domain);

/* Return the value of a field in a header entry.
   The return value is either a freshly allocated string, to be freed by the
   caller, or NULL.  */
extern char * po_header_field (const char *header, const char *field);

/* Return the header entry with a given field set to a given value.  The field
   is added if necessary.
   The return value is a freshly allocated string.  */
extern char * po_header_set_field (const char *header, const char *field, const char *value);


/* ======================= po_message_iterator_t API ======================= */

/* Create an iterator for traversing a domain of a PO file in memory.
   The domain NULL denotes the default domain.  */
extern po_message_iterator_t po_message_iterator (po_file_t file, const char *domain);

/* Free an iterator.  */
extern void po_message_iterator_free (po_message_iterator_t iterator);

/* Return the next message, and advance the iterator.
   Return NULL at the end of the message list.  */
extern po_message_t po_next_message (po_message_iterator_t iterator);

/* Insert a message in a PO file in memory, in the domain and at the position
   indicated by the iterator.  The iterator thereby advances past the freshly
   inserted message.  */
extern void po_message_insert (po_message_iterator_t iterator, po_message_t message);


/* =========================== po_message_t API ============================ */

/* Return a freshly constructed message.
   To finish initializing the message, you must set the msgid and msgstr.  */
extern po_message_t po_message_create (void);

/* Return the context of a message, or NULL for a message not restricted to a
   context.  */
extern const char * po_message_msgctxt (po_message_t message);

/* Change the context of a message. NULL means a message not restricted to a
   context.  */
extern void po_message_set_msgctxt (po_message_t message, const char *msgctxt);

/* Return the msgid (untranslated English string) of a message.  */
extern const char * po_message_msgid (po_message_t message);

/* Change the msgid (untranslated English string) of a message.  */
extern void po_message_set_msgid (po_message_t message, const char *msgid);

/* Return the msgid_plural (untranslated English plural string) of a message,
   or NULL for a message without plural.  */
extern const char * po_message_msgid_plural (po_message_t message);

/* Change the msgid_plural (untranslated English plural string) of a message.
   NULL means a message without plural.  */
extern void po_message_set_msgid_plural (po_message_t message, const char *msgid_plural);

/* Return the msgstr (translation) of a message.
   Return the empty string for an untranslated message.  */
extern const char * po_message_msgstr (po_message_t message);

/* Change the msgstr (translation) of a message.
   Use an empty string to denote an untranslated message.  */
extern void po_message_set_msgstr (po_message_t message, const char *msgstr);

/* Return the msgstr[index] for a message with plural handling, or
   NULL when the index is out of range or for a message without plural.  */
extern const char * po_message_msgstr_plural (po_message_t message, int index);

/* Change the msgstr[index] for a message with plural handling.
   Use a NULL value at the end to reduce the number of plural forms.  */
extern void po_message_set_msgstr_plural (po_message_t message, int index, const char *msgstr);

/* Return the comments for a message.  */
extern const char * po_message_comments (po_message_t message);

/* Change the comments for a message.
   comments should be a multiline string, ending in a newline, or empty.  */
extern void po_message_set_comments (po_message_t message, const char *comments);

/* Return the extracted comments for a message.  */
extern const char * po_message_extracted_comments (po_message_t message);

/* Change the extracted comments for a message.
   comments should be a multiline string, ending in a newline, or empty.  */
extern void po_message_set_extracted_comments (po_message_t message, const char *comments);

/* Return the i-th file position for a message, or NULL if i is out of
   range.  */
extern po_filepos_t po_message_filepos (po_message_t message, int i);

/* Remove the i-th file position from a message.
   The indices of all following file positions for the message are decremented
   by one.  */
extern void po_message_remove_filepos (po_message_t message, int i);

/* Add a file position to a message, if it is not already present for the
   message.
   file is the file name.
   start_line is the line number where the string starts, or (size_t)(-1) if no
   line number is available.  */
extern void po_message_add_filepos (po_message_t message, const char *file, size_t start_line);

/* Return the previous context of a message, or NULL for none.  */
extern const char * po_message_prev_msgctxt (po_message_t message);

/* Change the previous context of a message.  NULL is allowed.  */
extern void po_message_set_prev_msgctxt (po_message_t message, const char *prev_msgctxt);

/* Return the previous msgid (untranslated English string) of a message, or
   NULL for none.  */
extern const char * po_message_prev_msgid (po_message_t message);

/* Change the previous msgid (untranslated English string) of a message.
   NULL is allowed.  */
extern void po_message_set_prev_msgid (po_message_t message, const char *prev_msgid);

/* Return the previous msgid_plural (untranslated English plural string) of a
   message, or NULL for none.  */
extern const char * po_message_prev_msgid_plural (po_message_t message);

/* Change the previous msgid_plural (untranslated English plural string) of a
   message.  NULL is allowed.  */
extern void po_message_set_prev_msgid_plural (po_message_t message, const char *prev_msgid_plural);

/* Return true if the message is marked obsolete.  */
extern int po_message_is_obsolete (po_message_t message);

/* Change the obsolete mark of a message.  */
extern void po_message_set_obsolete (po_message_t message, int obsolete);

/* Return true if the message is marked fuzzy.  */
extern int po_message_is_fuzzy (po_message_t message);

/* Change the fuzzy mark of a message.  */
extern void po_message_set_fuzzy (po_message_t message, int fuzzy);

/* Return true if the message is marked as being a format string of the given
   type (e.g. "c-format").  */
extern int po_message_is_format (po_message_t message, const char *format_type);

/* Change the format string mark for a given type of a message.  */
extern void po_message_set_format (po_message_t message, const char *format_type, /*bool*/int value);

/* If a numeric range of a message is set, return true and store the minimum
   and maximum value in *MINP and *MAXP.  */
extern int po_message_is_range (po_message_t message, int *minp, int *maxp);

/* Change the numeric range of a message.  MIN and MAX must be non-negative,
   with MIN < MAX.  Use MIN = MAX = -1 to remove the numeric range of a
   message.  */
extern void po_message_set_range (po_message_t message, int min, int max);


/* =========================== po_filepos_t API ============================ */

/* Return the file name.  */
extern const char * po_filepos_file (po_filepos_t filepos);

/* Return the line number where the string starts, or (size_t)(-1) if no line
   number is available.  */
extern size_t po_filepos_start_line (po_filepos_t filepos);


/* ============================ Format type API ============================= */

/* Return a NULL terminated array of the supported format types.  */
extern const char * const * po_format_list (void);

/* Return the pretty name associated with a format type.
   For example, for "csharp-format", return "C#".
   Return NULL if the argument is not a supported format type.  */
extern const char * po_format_pretty_name (const char *format_type);


/* ============================= Checking API ============================== */

/* Test whether an entire file PO file is valid, like msgfmt does it.
   If it is invalid, pass the reasons to the handler.  */
extern void po_file_check_all (po_file_t file, po_xerror_handler_t handler);

/* Test a single message, to be inserted in a PO file in memory, like msgfmt
   does it.  If it is invalid, pass the reasons to the handler.  The iterator
   is not modified by this call; it only specifies the file and the domain.  */
extern void po_message_check_all (po_message_t message, po_message_iterator_t iterator, po_xerror_handler_t handler);

/* Test whether the message translation is a valid format string if the message
   is marked as being a format string.  If it is invalid, pass the reasons to
   the handler.  */
#define po_message_check_format po_message_check_format_v2
extern void po_message_check_format (po_message_t message, po_xerror_handler_t handler);


#ifdef __cplusplus
}
#endif

#endif /* _GETTEXT_PO_H */
