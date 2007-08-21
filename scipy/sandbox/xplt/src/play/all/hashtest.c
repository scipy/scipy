#include "phash.h"
#include "pstdlib.h"
#include "pstdio.h"
#include <string.h>

static int mxsyms = 50000;
static void test_delete(void *value);
static char nm1[] = "  name 1  ";
static char nm2[] = "  name 2  ";

int
main(int argc, char *argv[])
{
  p_hashkey *id, collist[32];
  p_hashtab *tab;
  int i, ncoll;
  int nsyms = 0;
  char **syms = p_malloc(sizeof(char *)*mxsyms);

  if (argc>1) {
    char line[120];
    FILE *f = fopen(argv[1], "r");
    if (!f) {
      printf("hashtest: unable to open %s\n", argv[1]);
      exit(1);
    }
    for (nsyms=0 ; nsyms<mxsyms ; nsyms++) {
      if (!fgets(line, 120, f)) break;
      i = strlen(line);
      if (i<1) continue;
      if (line[i-1]=='\n') {
        i--;
        line[i] = '\0';
        if (i<1) continue;
      }
      syms[nsyms] = p_strcpy(line);
    }
    fclose(f);
  } else {
    syms[0] = p_strcpy("first");
    syms[1] = p_strcpy("second");
    syms[2] = p_strcpy("third");
    syms[3] = p_strcpy("fourth");
    syms[4] = p_strcpy("fifth");
    nsyms = 5;
  }

  tab = p_halloc(4);
  for (i=0 ; i<nsyms ; i++) p_hinsert(tab, P_PHASH(syms[i]), syms[i]);
  for (i=0 ; i<nsyms ; i++)
    if (p_hfind(tab, P_PHASH(syms[i])) != syms[i]) {
      puts("hashtab: p_hfind failed first test");
      exit(1);
    }

  for (i=0 ; i<nsyms/2 ; i++) p_hinsert(tab, P_PHASH(syms[i]), (void *)0);
  for (i=0 ; i<nsyms/2 ; i++)
    if (p_hfind(tab, P_PHASH(syms[i]))) {
      puts("hashtab: p_hfind failed first removal test");
      exit(1);
    }
  for ( ; i<nsyms ; i++)
    if (p_hfind(tab, P_PHASH(syms[i])) != syms[i]) {
      puts("hashtab: p_hfind failed second removal test");
      exit(1);
    }

  p_hfree(tab, &test_delete);
  for (i=0 ; i<nsyms ; i++) {
    if (((syms[i][0] & 0x80)==0) == (i<nsyms/2)) syms[i][0] &= 0x7f;
    else {
      puts("hashtab: p_hfree deletion callback failed");
      exit(1);
    }
  }

  id = p_malloc(sizeof(p_hashkey)*(nsyms+2));
  for (i=0 ; i<nsyms ; i++) id[i] = p_idmake(syms[i],
                                             (i&1)? strlen(syms[i]) : 0);
  for (i=0 ; i<nsyms ; i++)
    if (p_id(syms[i], 0) != id[i]) {
      puts("hashtab: p_id failed first test");
      exit(1);
    }
  id[nsyms] = p_idstatic(nm1);
  id[nsyms+1] = p_idstatic(nm2);
  for (i=0 ; i<nsyms ; i++)
    if (strcmp(p_idname(id[i]), syms[i])) {
      puts("hashtab: p_idname failed first test");
      exit(1);
    }
  if (p_id("#%  #$&\t^",0)) {
    puts("hashtab: p_id failed second test");
    exit(1);
  }
  if (p_id(nm1,0)!=id[nsyms] ||
      p_id(nm2,0)!=id[nsyms+1]) {
    puts("hashtab: p_id failed third test");
    exit(1);
  }
  if (strcmp(p_idname(id[nsyms]),nm1) ||
      strcmp(p_idname(id[nsyms+1]),nm2)) {
    puts("hashtab: p_idname failed third test");
    exit(1);
  }
  for (i=0 ; i<nsyms ; i+=2)
    if (p_idmake(syms[i],0)!=id[i]) {
      puts("hashtab: p_idmake failed repeat test");
      exit(1);
    }
  for (i=0 ; i<nsyms ; i+=3) {
    p_idfree(id[i]);
    if (i&1) id[i] = 0;
  }
  for (i=0 ; i<nsyms ; i++)
    if (p_id(syms[i], 0) != id[i]) {
      puts("hashtab: p_id failed fourth test");
      exit(1);
    }
  for (i=0 ; i<nsyms ; i++) if (id[i]) p_idfree(id[i]);
  ncoll = 0;
  for (i=0 ; i<nsyms ; i++)
    if ((p_id(syms[i], 0)!=0) == (i&1 || !((i/2)%3))) {
      if (ncoll<32) collist[ncoll] = id[i];
      ncoll++;
    }
  {
    extern int p_id_collisions;
    printf("hashtab: p_id had %d collisions out of %d symbols",
           p_id_collisions, nsyms);
  }

  for (i=0 ; i<nsyms ; i++) p_setctx(syms[i], syms[i]+1);
  for (i=0 ; i<nsyms ; i++)
    if (p_getctx(syms[i]) != syms[i]+1) {
      puts("hashtab: p_getctx failed first test");
      exit(1);
    }
  for (i=0 ; i<nsyms ; i+=2) p_setctx(syms[i], (void *)0);
  for (i=0 ; i<nsyms ; i++)
    if (p_getctx(syms[i]) != ((i&1)? syms[i]+1 : 0)) {
      puts("hashtab: p_getctx failed second test");
      exit(1);
    }

  for (i=0 ; i<nsyms ; i++) p_free(syms[i]);
  p_free(syms);
  exit(0);
  return 0;
}

static void
test_delete(void *value)
{
  char *sym = value;
  sym[0] |= 0x80;
}
