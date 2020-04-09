#include <cstdio>
#include <cstring>

#include "HighsIO.h"
#include "catch.hpp"

char printedmsg[100000];
void* receiveddata = NULL;

// callback that saves message away for comparison
static void myprintmsgcb(int level, const char* msg, void* msgcb_data) {
  strcpy(printedmsg, msg);
  receiveddata = msgcb_data;
}

static void mylogmsgcb(HighsMessageType type, const char* msg,
                       void* msgcb_data) {
  strcpy(printedmsg, msg);
  receiveddata = msgcb_data;
}

TEST_CASE("msgcb", "[highs_io]") {
  int dummydata = 42;

  HighsSetMessageCallback(myprintmsgcb, mylogmsgcb, (void*)&dummydata);

  int message_level = ML_MINIMAL;
  HighsPrintMessage(stdout, message_level, 4, "Hi %s!", "HiGHS");
  REQUIRE(strcmp(printedmsg, "Hi HiGHS!") == 0);
  REQUIRE(receiveddata == &dummydata);

  /* printed at level 4 when level is 3 should not print */
  *printedmsg = '\0';
  message_level = 3;
  HighsPrintMessage(stdout, message_level, 4, "Hi %s!", "HiGHS");
  REQUIRE(*printedmsg == '\0');

  {
    char longmsg[sizeof(printedmsg)];
    memset(longmsg, 'H', sizeof(longmsg));
    longmsg[sizeof(longmsg) - 1] = '\0';
    HighsPrintMessage(stdout, message_level, 2, longmsg);
    REQUIRE(strncmp(printedmsg, "HHHH", 4) == 0);
    REQUIRE(strlen(printedmsg) <= sizeof(printedmsg));
  }

  HighsLogMessage(stdout, HighsMessageType::INFO, "Hello %s!", "HiGHS");
  REQUIRE(strlen(printedmsg) > 8);
  REQUIRE(strcmp(printedmsg + 8, " [INFO   ] Hello HiGHS!\n") ==
          0);  // begin of printedmsg is a timestamp, which we skip over
  REQUIRE(receiveddata == &dummydata);

  {
    char longmsg[sizeof(printedmsg)];
    memset(longmsg, 'H', sizeof(longmsg));
    longmsg[sizeof(longmsg) - 1] = '\0';
    HighsLogMessage(stdout, HighsMessageType::WARNING, longmsg);
    REQUIRE(strstr(printedmsg, "HHHH") != NULL);
    REQUIRE(strlen(printedmsg) <= sizeof(printedmsg));
  }
}
