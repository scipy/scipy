#!/usr/bin/env python

import sys
import argparse
import dbus
import time

def get_cmdline(pid):
  cmdline = ''
  if pid > 0:
    try:
      procpath = '/proc/' + str(pid) + '/cmdline'
      with open(procpath, 'r') as f:
        cmdline = " ".join(f.readline().split('\0'))
    except:
      pass
  return cmdline

# Parsing parameters

parser = argparse.ArgumentParser(description='Testing D-Bus match rules')
parser.add_argument('--session', help='session bus', action="store_true")
parser.add_argument('--system', help='system bus', action="store_true")
parser.add_argument('--all', help='print all match rules', action="store_true")
args = parser.parse_args()

if args.system and args.session:
  parser.print_help()
  sys.exit(1)

# Fetch data from the bus driver

if args.system:
  bus = dbus.SystemBus()
else:
  bus = dbus.SessionBus()

remote_object = bus.get_object("org.freedesktop.DBus",
                               "/org/freedesktop/DBus")
bus_iface = dbus.Interface(remote_object, "org.freedesktop.DBus")
stats_iface = dbus.Interface(remote_object, "org.freedesktop.DBus.Debug.Stats")

try:
  match_rules = stats_iface.GetAllMatchRules()
except:
  print("GetConnectionMatchRules failed: did you enable the Stats interface?")
  sys.exit(1)

names = bus_iface.ListNames()
unique_names = [ a for a in names if a.startswith(":") ]
pids = dict((name, bus_iface.GetConnectionUnixProcessID(name)) for name in unique_names)
cmds = dict((name, get_cmdline(pids[name])) for name in unique_names)
well_known_names = [ a for a in names if a not in unique_names ]
owners = dict((wkn, bus_iface.GetNameOwner(wkn)) for wkn in well_known_names)

rules = dict((k_rules,
              dict({
                'wkn': [k for k, v in owners.items() if v == k_rules],
                'pid': pids[k_rules],
                'cmd': cmds[k_rules] or "",
                'rules': v_rules,
                'warnings': dict({
                  'not_signal': [a for a in v_rules if "type='signal'" not in a],
                  'no_sender': [a for a in v_rules if "sender=" not in a],
                  'local': [a for a in v_rules if "org.freedesktop.DBus.Local" in a],
                  'NameOwnerChanged_arg0': [a for a in v_rules if "member='NameOwnerChanged'" in a and "arg0" not in a]
                })
              })
             ) for k_rules, v_rules in match_rules.items())

warnings = dict({
             'not_signal': 'Match rule without selecting signals',
             'no_sender': 'Match rule without a sender criteria',
             'local': 'Match rule on the org.freedesktop.DBus.Local interface',
             'NameOwnerChanged_arg0': 'Match rule on NameOwnerChanged without a arg0* criteria'
           })

# Print the match rules

# print all match rules without analysing them
if args.all:
  for name in rules:
    print("Connection %s with pid %d '%s' (%s): %d match rules, %d warnings"
          % (name, rules[name]['pid'], rules[name]['cmd'],
             ' '.join(rules[name]['wkn']), len(rules[name]['rules']),
             len(sum(rules[name]['warnings'].values(), []))))
    for rule in rules[name]['rules']:
      print("\t%s" % (rule))
    print("")
  sys.exit(0)

# analyse match rules and print only the suspicious ones
for conn,data in rules.items():
  warnings_count = len(sum(data['warnings'].values(), []))
  if warnings_count == 0:
    continue

  print("Connection %s with pid %d '%s' (%s): %d match rules, %d warnings"
        % (conn, data['pid'], data['cmd'], ' '.join(data['wkn']),
           len(data['rules']), warnings_count))

  for warn_code,rule_list in [(warn_code,rule_list) \
                              for warn_code, rule_list \
                              in data['warnings'].items() \
                              if len(rule_list) > 0]:
    print("   - %s:" % (warnings[warn_code]))
    for rule in rule_list:
      print("         - %s" % (rule))

  print("")

sys.exit(0)
