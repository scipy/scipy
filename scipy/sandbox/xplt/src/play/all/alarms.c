/*
 * alarms.c -- $Id$
 * alarm event functions, implemented using play interface
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "play.h"
#include "pstdlib.h"

typedef struct p_alarm p_alarm;
struct p_alarm {
  p_alarm *next;
  double time;
  void (*on_alarm)(void *c);
  void *context;
};

static p_alarm *alarm_next = 0;
static p_alarm *alarm_free = 0;
static double alarm_query(void);
static int idle_eligible = 1;

static int p_dflt_idle(void);
static int (*p_app_idle)(void)= &p_dflt_idle;

void
p_idler(int (*on_idle)(void))
{
  p_app_idle = on_idle;
}

static int
p_dflt_idle(void)
{
  return 0;
}

void
p_on_idle(int reset)
{
  if (!reset) {
    if (alarm_next && !alarm_query()) {
      /* alarm has rung - unlink it and call its on_alarm */
      p_alarm *next = alarm_next;
      alarm_next = next->next;
      next->next = alarm_free;
      alarm_free = next;
      next->on_alarm(next->context);
      idle_eligible = 1;
    } else {
      idle_eligible = p_app_idle();
    }
  } else {
    idle_eligible = 1;
  }
}

double
p_timeout(void)
{
  int eligible = idle_eligible;
  idle_eligible = 1;
  return eligible? 0.0 : (alarm_next? alarm_query() : -1.0);
}

void
p_set_alarm(double secs, void (*on_alarm)(void *c), void *context)
{
  p_alarm *me;
  p_alarm *next = alarm_next;
  p_alarm **prev = &alarm_next;
  double time;
  if (!alarm_free) {
    int n = 8;
    alarm_free = p_malloc(sizeof(p_alarm)*n);
    alarm_free[--n].next = 0;
    while (n--) alarm_free[n].next = &alarm_free[n+1];
  }
  me = alarm_free;
  me->time = time = p_wall_secs() + secs;
  me->on_alarm = on_alarm;
  me->context = context;
  /* insert me into alarm_next list, kept in order of time */
  while (next && next->time<=time) {
    prev = &next->next;
    next = next->next;
  }
  alarm_free = alarm_free->next;
  me->next = next;
  *prev = me;
}

void
p_clr_alarm(void (*on_alarm)(void *c), void *context)
{
  p_alarm *next, **prev = &alarm_next;
  for (next=alarm_next ; next ; next=*prev) {
    if ((!on_alarm || on_alarm==next->on_alarm) &&
        (!context || context==next->context)) {
      *prev = next->next;
      next->next = alarm_free;
      alarm_free = next;
    } else {
      prev = &next->next;
    }
  }
}

static double
alarm_query(void)
{
  if (alarm_next->time != -1.e35) {
    double time = p_wall_secs();
    p_alarm *next = alarm_next;
    /* if no alarms need to ring yet, return earliest */
    if (next->time > time)
      return next->time - time;
    do {
      next->time = -1.e35;   /* mark all alarms that need to ring */
      next = next->next;
    } while (next && next->time<=time);
  }
  return 0.0;
}
