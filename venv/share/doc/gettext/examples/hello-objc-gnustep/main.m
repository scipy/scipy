/* Example for use of GNU gettext.
   This file is in the public domain.

   Source code of the main program.  */

#include <AppKit/AppKit.h>
#include "AppController.h"

#define APP_NAME @"Hello"

/* Create the application's menu.  */
static void
createMenu ()
{
  NSMenu *menu;
  NSMenu *info;
  NSMenu *edit;
  NSMenu *services;
  NSMenu *windows;

  SEL action = @selector(method:);

  menu = [[NSMenu alloc] initWithTitle: APP_NAME];
  [menu addItemWithTitle: @"Info"
        action: action
        keyEquivalent: @""];
  [menu addItemWithTitle: @"Edit"
        action: action
        keyEquivalent: @""];
  [menu addItemWithTitle: @"Hello..."
        action: @selector(showHelloWindow:)
        keyEquivalent: @""];
  [menu addItemWithTitle: @"Windows"
        action: action
        keyEquivalent: @""];
  [menu addItemWithTitle: @"Services"
        action: action
        keyEquivalent: @""];
  [menu addItemWithTitle: @"Hide"
        action: @selector(hide:)
        keyEquivalent: @"h"];
  [menu addItemWithTitle: @"Quit"
        action: @selector(terminate:)
        keyEquivalent: @"q"];

  info = AUTORELEASE ([[NSMenu alloc] init]);
  [info addItemWithTitle: @"Info Panel..."
        action: @selector(showInfoPanel:)
        keyEquivalent: @""];
  [info addItemWithTitle: @"Preferences"
        action: @selector(showPrefPanel:)
        keyEquivalent: @""];
  [info addItemWithTitle: @"Help"
        action: action
        keyEquivalent: @"?"];
  [menu setSubmenu: info forItem: [menu itemWithTitle: @"Info"]];

  edit = AUTORELEASE ([[NSMenu alloc] init]);
  [edit addItemWithTitle: @"Cut"
        action: @selector(cut:)
        keyEquivalent: @"x"];
  [edit addItemWithTitle: @"Copy"
        action: @selector(copy:)
        keyEquivalent: @"c"];
  [edit addItemWithTitle: @"Paste"
        action: @selector(paste:)
        keyEquivalent: @"v"];
  [edit addItemWithTitle: @"Delete"
        action: @selector(delete:)
        keyEquivalent: @""];
  [edit addItemWithTitle: @"Select All"
        action: @selector(selectAll:)
        keyEquivalent: @"a"];
  [menu setSubmenu: edit forItem: [menu itemWithTitle: @"Edit"]];

  windows = AUTORELEASE ([[NSMenu alloc] init]);
  [windows addItemWithTitle: @"Arrange"
           action: @selector(arrangeInFront:)
           keyEquivalent: @""];
  [windows addItemWithTitle: @"Miniaturize"
           action: @selector(performMiniaturize:)
           keyEquivalent: @"m"];
  [windows addItemWithTitle: @"Close"
           action: @selector(performClose:)
           keyEquivalent: @"w"];
  [menu setSubmenu: windows forItem: [menu itemWithTitle: @"Windows"]];

  services = AUTORELEASE ([[NSMenu alloc] init]);
  [menu setSubmenu: services forItem: [menu itemWithTitle: @"Services"]];

  [[NSApplication sharedApplication] setMainMenu: menu];
  [[NSApplication sharedApplication] setServicesMenu: services];
  [[NSApplication sharedApplication] setWindowsMenu: windows];
}

/* Initialise and go!  */
int
main(int argc, const char *argv[])
{
  NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
  AppController *controller;

  [NSApplication sharedApplication];

  createMenu ();

  controller = [[AppController alloc] init];
  [NSApp setDelegate:controller];

  RELEASE (pool);

  return NSApplicationMain (argc, argv);
}
