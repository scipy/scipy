/*
  Copyright (C) 2009  Andreas Gruenbacher <agruen@suse.de>

  This program is free software: you can redistribute it and/or modify it
  under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __ERROR_CONTEXT_T
#define __ERROR_CONTEXT_T

#ifdef __cplusplus
extern "C" {
#endif

struct error_context {
	/* Process an error message */
	void (*error) (struct error_context *, const char *, ...);

	/* Quote a file name for including in an error message */
	const char *(*quote) (struct error_context *, const char *);

	/* Free a quoted name */
	void (*quote_free) (struct error_context *, const char *);
};

#ifdef ERROR_CONTEXT_MACROS
# define error(ctx, args...) do { \
	if ((ctx) && (ctx)->error) \
		(ctx)->error((ctx), args); \
	} while(0)
# define quote(ctx, name) \
	( ((ctx) && (ctx)->quote) ? (ctx)->quote((ctx), (name)) : (name) )
# define quote_free(ctx, name) do { \
	if ((ctx) && (ctx)->quote_free) \
		(ctx)->quote_free((ctx), (name)); \
	} while(0)
#endif

#ifdef __cplusplus
}
#endif

#endif  /* __ERROR_CONTEXT_T */
