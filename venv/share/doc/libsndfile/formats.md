---
layout: page
title: libsndfile: Supported formats
---

# libsndfile : Supported formats

The following table lists the file formats and encodings that libsndfile can
read and write. The file formats are arranged across the top and encodings along
the left edge.

{:.formats-table}
|                    | Microsoft WAV | SGI / Apple AIFF / AIFC |Sun / DEC / NeXT AU / SND | Headerless RAW | Paris Audio File PAF | Commodore Amiga IFF / SVX | Sphere Nist WAV | IRCAM SF | Creative VOC | Soundforge W64 | GNU Octave 2.0 MAT4 | GNU Octave 2.1 MAT5 | Portable Voice Format PVF | Fasttracker 2 XI | HMM Tool Kit HTK | Apple CAF | Sound Designer II SD2 | Free Lossless Audio Codec FLAC |
|:-------------------|:-------------:|:-----------------------:|:------------------------:|:--------------:|:--------------------:|:-------------------------:|:---------------:|:--------:|:------------:|:--------------:|:-------------------:|:-------------------:|:-------------------------:|:----------------:|:----------------:|:---------:|:---------------------:|:------------------------------:|
| Unsigned 8 bit PCM | R/W           | R/W                     |                          | R/W            |                      |                           |                 |          | R/W          | R/W            |                     | R/W                 |                           |                  |                  |           |                       |                                |
| Signed 8 bit PCM   |               | R/W                     | R/W                      | R/W            | R/W                  | R/W                       | R/W             |          |              |                |                     |                     | R/W                       |                  |                  | R/W       | R/W                   | R/W                            |
| Signed 16 bit PCM  | R/W           | R/W                     | R/W                      | R/W            | R/W                  | R/W                       | R/W             | R/W      | R/W          | R/W            | R/W                 | R/W                 | R/W                       |                  | R/W              | R/W       | R/W                   | R/W                            |
| Signed 24 bit PCM  | R/W           | R/W                     | R/W                      | R/W            | R/W                  |                           | R/W             |          |              | R/W            |                     |                     |                           |                  |                  | R/W       | R/W                   | R/W                            |
| Signed 32 bit PCM  | R/W           | R/W                     | R/W                      | R/W            |                      |                           | R/W             | R/W      |              | R/W            | R/W                 | R/W                 | R/W                       |                  |                  | R/W       |                       |                                |
| 32 bit float       | R/W           | R/W                     | R/W                      | R/W            |                      |                           |                 | R/W      |              | R/W            | R/W                 | R/W                 |                           |                  |                  | R/W       |                       |                                |
| 64 bit double      | R/W           | R/W                     | R/W                      | R/W            |                      |                           |                 |          |              | R/W            | R/W                 | R/W                 |                           |                  |                  | R/W       |                       |                                |
| u-law encoding     | R/W           | R/W                     | R/W                      | R/W            |                      |                           | R/W             | R/W      | R/W          | R/W            |                     |                     |                           |                  |                  | R/W       |                       |                                |
| A-law encoding     | R/W           | R/W                     | R/W                      | R/W            |                      |                           | R/W             | R/W      | R/W          | R/W            |                     |                     |                           |                  |                  | R/W       |                       |                                |
| IMA ADPCM          | R/W           |                         |                          |                |                      |                           |                 |          |              | R/W            |                     |                     |                           |                  |                  |           |                       |                                |
| MS ADPCM           | R/W           |                         |                          |                |                      |                           |                 |          |              | R/W            |                     |                     |                           |                  |                  |           |                       |                                |
| GSM 6.10           | R/W           | R/W                     |                          | R/W            |                      |                           |                 |          |              | R/W            |                     |                     |                           |                  |                  |           |                       |                                |
| G721 ADPCM 32kbps  | R/W           |                         | R/W                      |                |                      |                           |                 |          |              |                |                     |                     |                           |                  |                  |           |                       |                                |
| G723 ADPCM 24kbps  |               |                         | R/W                      |                |                      |                           |                 |          |              |                |                     |                     |                           |                  |                  |           |                       |                                |
| G723 ADPCM 40kbps  |               |                         | R/W                      |                |                      |                           |                 |          |              |                |                     |                     |                           |                  |                  |           |                       |                                |
| 12 bit DWVW        |               | R/W                     |                          | R/W            |                      |                           |                 |          |              |                |                     |                     |                           |                  |                  |           |                       |                                |
| 16 bit DWVW        |               | R/W                     |                          | R/W            |                      |                           |                 |          |              |                |                     |                     |                           |                  |                  |           |                       |                                |
| 24 bit DWVW        |               | R/W                     |                          | R/W            |                      |                           |                 |          |              |                |                     |                     |                           |                  |                  |           |                       |                                |
| Ok Dialogic ADPCM  |               |                         |                          | R/W            |                      |                           |                 |          |              |                |                     |                     |                           |                  |                  |           |                       |                                |
| 8 bit DPCM         |               |                         |                          |                |                      |                           |                 |          |              |                |                     |                     |                           | R/W              |                  |           |                       |                                |
| 16 bit DPCM        |               |                         |                          |                |                      |                           |                 |          |              |                |                     |                     |                           | R/W              |                  |           |                       |                                |

From version 1.0.18, libsndfile also reads and writes
[FLAC](https://xiph.org/flac/) and [Ogg/Vorbis](https://xiph.org/vorbis/).

From version 1.0.29, libsndfile can read and write
[Ogg/Opus](https://opus-codec.org/).

From version 1.1.0, libsndfile can read and write MP3.

Some of the file formats I am also interested in adding are:

- Kurzweil K2000 sampler files.
- Ogg Speex.

Other file formats may also be added on request.

If you are interested in how to add a new format to a libsndfile, you may find
this [FAQ](new_file_type_howto.md) helpful.
