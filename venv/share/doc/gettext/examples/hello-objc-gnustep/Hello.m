/* Example for use of GNU gettext.
   This file is in the public domain.

   Source code of the Hello class.  */

#include <unistd.h>
#include "Hello.h"
#include <GNUstepGUI/GSHbox.h>
#include <GNUstepGUI/GSVbox.h>

@implementation Hello

- (id)init
{
  if ((self = [super init]))
    [self createUI];
  return self;
}

- (void)dealloc
{
  RELEASE (window);
  [super dealloc];
}

- (void)makeKeyAndOrderFront
{
  if (![window isVisible])
    [window center];
  [window makeKeyAndOrderFront:self];
}

- (void)done
{
  [window close];
}

@end

@implementation Hello (UIBuilder)

- (void)createUI
{
  GSVbox *cview;
  GSHbox *buttonbar;
  int i;

  label1 = [NSTextField new];
  [label1 setStringValue: _(@"Hello, world!")];
  [label1 setAlignment: NSLeftTextAlignment];
  [label1 setBordered: NO];
  [label1 setEditable: NO];
  [label1 setBezeled: NO];
  [label1 setDrawsBackground: NO];
  [label1 sizeToFit];

  label2 = [NSTextField new];
  [label2 setStringValue: [NSString stringWithFormat: _(@"This program is running as process number %d."), [[NSProcessInfo processInfo] processIdentifier]]];
  [label2 setAlignment: NSLeftTextAlignment];
  [label2 setBordered: NO];
  [label2 setEditable: NO];
  [label2 setBezeled: NO];
  [label2 setDrawsBackground: NO];
  [label2 sizeToFit];

  okButton = [NSButton new];
  [okButton setTitle: @"OK"];
  [okButton setTarget: self];
  [okButton setAction: @selector(done)];
  [okButton setFrameSize: NSMakeSize(60,22)];
  [okButton setAutoresizingMask: 7];

  buttonbar = [GSHbox new];
  [buttonbar setAutoresizingMask: NSViewMinXMargin];
  [buttonbar addView: okButton];
  AUTORELEASE (okButton);

  cview = [GSVbox new];
  // GSVbox is flawed: We have to add the controls bottom-up, and mark the
  // last one (= the topmost one) as non-resizable in Y direction, so that the
  // Y space becomes equally distributed _between_ (not above) the subviews.
  [cview addView: buttonbar];
  AUTORELEASE (buttonbar);
  [cview addView: label2];
  AUTORELEASE (label2);
  [cview addView: label1 enablingYResizing: NO];
  AUTORELEASE (label1);

  window = [[NSWindow alloc] initWithContentRect: NSMakeRect(0,0, [cview frame].size.width, [cview frame].size.height)
                             styleMask: (NSTitledWindowMask | NSClosableWindowMask | NSMiniaturizableWindowMask)
                             backing: NSBackingStoreBuffered
                             defer: NO];
  [window setDelegate: self];
  [window setTitle: @"Hello example"];
  [window setReleasedWhenClosed: NO];
  [window setContentView: cview];
}

@end
