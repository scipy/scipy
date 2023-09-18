/* Example for use of GNU gettext.
   This file is in the public domain.

   Source code of the AppController class.  */

#include "AppController.h"
#include "Hello.h"

@implementation AppController

static NSDictionary *infoDict = nil;

+ (void)initialize
{
  NSMutableDictionary *defaults = [NSMutableDictionary dictionary];

  [[NSUserDefaults standardUserDefaults] registerDefaults: defaults];
  [[NSUserDefaults standardUserDefaults] synchronize];
}

- (id)init
{
  self = [super init];
  return self;
}

- (void)dealloc
{
  if (hello)
    RELEASE (hello);

  [super dealloc];
}

- (void)awakeFromNib
{
}

- (void)applicationDidFinishLaunching:(NSNotification *)notif
{
}

- (BOOL)applicationShouldTerminate:(id)sender
{
  return YES;
}

- (void)applicationWillTerminate:(NSNotification *)notification
{
}

- (BOOL)application:(NSApplication *)application openFile:(NSString *)fileName
{
}

- (void)showPrefPanel:(id)sender
{
}

- (void)showInfoPanel:(id)sender
{
  if (!infoDict)
    {
      NSString *fp;
      NSBundle *bundle = [NSBundle mainBundle];

      fp = [bundle pathForResource: @"Info-project" ofType: @"plist"];
      infoDict = [[NSDictionary dictionaryWithContentsOfFile: fp] retain];
    }
  [[NSApplication sharedApplication] orderFrontStandardInfoPanelWithOptions: infoDict];
}

- (void)showHelloWindow:(id)sender
{
  if (!hello)
    hello = [[Hello alloc] init];

  [hello makeKeyAndOrderFront];
}

@end
