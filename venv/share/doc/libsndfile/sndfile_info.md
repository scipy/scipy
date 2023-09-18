---
layout: page
title: sndfile-info
---

Here is an example of the output from the **sndfile-info** program
distributed with libsndfile.

This file was opened and parsed correctly but had been truncated so that
the values in the **FORM** and **SSND** chunks were incorrect.

```
        erikd@hendrix > examples/sndfile-info truncated.aiff
        truncated.aiff
        size : 200000
        FORM : 307474 (should be 199992)
         AIFF
         COMM : 18
          Sample Rate : 16000
          Samples     : 76857
          Channels    : 2
          Sample Size : 16
         SSND : 307436 (should be 199946)
          Offset     : 0
          Block Size : 0

        --------------------------------
        Sample Rate : 16000
        Frames      : 76857
        Channels    : 2
        Bit Width   : 16
        Format      : 0x00020001
        Sections    : 1
        Seekable    : TRUE
        Signal Max  : 32766
```
