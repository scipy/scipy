/***
  This file is part of PulseAudio.

  Copyright 2012 Alexander Kurtz <kurtz.alex@googlemail.com>

  PulseAudio is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published
  by the Free Software Foundation; either version 2.1 of the License,
  or (at your option) any later version.

  PulseAudio is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with PulseAudio; if not, see <http://www.gnu.org/licenses/>.
***/

namespace PulseAudio {
        [Compact]
        [CCode (cheader_filename="pulse/simple.h", cname="pa_simple", cprefix="pa_simple_")]
        class Simple {
                public Simple(string? server = null, string? name = null, Stream.Direction dir = Stream.Direction.PLAYBACK,
                              string? dev = null, string stream_name = "",
                              SampleSpec ss = SampleSpec(){ format = SampleFormat.S16NE, rate = 44100, channels = 2 },
                              ChannelMap? map = null, Stream.BufferAttr? attr = null, out int error = null);
                public int write(void* data, size_t bytes, out int error = null);
                public int drain(out int error = null);
                public int read(void* data, size_t bytes, out int error = null);
                public usec get_latency(out int error = null);
                public int flush(out int error = null);
        }
}
