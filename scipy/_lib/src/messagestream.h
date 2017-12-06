#ifndef MESSAGESTREAM_H_
#define MESSAGESTREAM_H_

#include <stdio.h>

#include "messagestream_config.h"

#if HAVE_OPEN_MEMSTREAM
FILE *messagestream_open_memstream(char **ptr, size_t *sizeloc)
{
    return open_memstream(ptr, sizeloc);
}
#else
FILE *messagestream_open_memstream(char **ptr, size_t *sizeloc)
{
    return NULL;
}
#endif

#endif /* MESSAGESTREAM_H_ */
