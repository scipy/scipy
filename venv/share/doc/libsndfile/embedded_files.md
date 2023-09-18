---
layout: page
---

# Embedded Sound Files

By using the open SNDFILE with a file descriptor function:

    SNDFILE*  sf_open_fd (int fd, int mode, SF_INFO *sfinfo, int close_desc) ;

it is possible to open sound files embedded within larger files. There are
however a couple of caveats:

* Read/Write mode (SFM_RDWR) is not supported.
* Writing of embedded files is only supported at the end of the file.
* Reading of embedded files is only supported at file offsets greater than zero.
* Not all file formats are supported (currently only WAV, AIFF and AU).

The test program **multi_file_test.c** in the **tests/** directory of the source
code tarball shows how this functionality is used to read and write embedded
files.
