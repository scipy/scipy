/***
  This file is part of PulseAudio.

  Copyright 2009 Lennart Poettering

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

using GLib;
using Posix;

[CCode (cheader_filename="pulse/pulseaudio.h")]
namespace PulseAudio {

        [CCode (cname="pa_get_library_version")]
        public unowned string get_library_version();

        [CCode (cname="PA_API_VERSION")]
        public const int API_VERSION;

        [CCode (cname="PA_PROTOCOL_VERSION")]
        public const int PROTOCOL_VERSION;

        [CCode (cname="PA_MAJOR")]
        public const int MAJOR;

        [CCode (cname="PA_MINOR")]
        public const int MINOR;

        [CCode (cname="PA_MICRO")]
        public const int MICRO;

        [CCode (cname="PA_CHECK_VERSION")]
        public bool CHECK_VERSION(int major, int minor, int micro);

        [CCode (cname="PA_INVALID_INDEX")]
        public const uint32 INVALID_INDEX;

        [CCode (cname="pa_free_cb_t", has_target=false)]
        public delegate void FreeCb(void *p);

        [CCode (cname="pa_sample_format_t", cprefix="PA_SAMPLE_", has_type_id=false)]
        public enum SampleFormat {
                U8,
                ALAW,
                ULAW,
                S16LE,
                S16BE,
                FLOAT32LE,
                FLOAT32BE,
                S32LE,
                S32BE,
                S24LE,
                S24BE,
                S24_32LE,
                S24_32BE,
                MAX,
                S16NE,
                S16RE,
                FLOAT32NE,
                FLOAT32RE,
                S32NE,
                S32RE,
                S24NE,
                S24RE,
                S24_32NE,
                S24_32RE;

                [CCode (cname="pa_sample_size_of_format")]
                public size_t size();

                [CCode (cname="pa_sample_format_to_string")]
                public unowned string? to_string();

                [CCode (cname="pa_sample_format_is_le")]
                public int is_le();

                [CCode (cname="pa_sample_format_is_be")]
                public int is_be();

                [CCode (cname="pa_sample_format_is_ne")]
                public int is_ne();

                [CCode (cname="pa_sample_format_is_re")]
                public int is_re();

                [CCode (cname="pa_parse_sample_format")]
                public static SampleFormat parse(string b);
        }

        [CCode (cname="pa_usec_t")]
        public struct usec : uint64 {
        }

        [CCode (cname="pa_sample_spec", has_type_id=false)]
        public struct SampleSpec {
                public SampleFormat format;
                public uint32 rate;
                public uint8 channels;

                [CCode (cname="PA_SAMPLE_SPEC_SNPRINT_MAX")]
                public const size_t SNPRINT_MAX;

                [CCode (cname="pa_bytes_per_second")]
                public size_t bytes_per_second();

                [CCode (cname="pa_frame_size")]
                public size_t frame_size();

                [CCode (cname="pa_sample_size")]
                public size_t sample_size();

                [CCode (cname="pa_bytes_to_usec", instance_pos=1.1)]
                public usec bytes_to_usec(size_t size);

                [CCode (cname="pa_usec_to_bytes", instance_pos=1.1)]
                public size_t usec_to_bytes(usec u);

                [CCode (cname="pa_sample_spec_init")]
                public unowned SampleSpec? init();

                [CCode (cname="pa_sample_spec_valid")]
                public bool valid();

                [CCode (cname="pa_sample_spec_equal")]
                public bool equal(SampleSpec other);

                [CCode (cname="pa_sample_spec_snprint", instance_pos=3.1)]
                public unowned string snprint(char[] buf);

                public string sprint() {
                        var buffer = new char[SNPRINT_MAX];
                        this.snprint(buffer);
                        return (string) buffer;
                }

                public string to_string() {
                        return sprint();
                }

                [CCode (cname="pa_sample_spec_init")]
                public SampleSpec();
        }

        // [CCode (cname="PA_BYTES_SNPRINT_MAX")]
        [CCode (cname="PA_SAMPLE_SPEC_SNPRINT_MAX")]
        public const size_t BYTES_SNPRINT_MAX;

        [CCode (cname="pa_bytes_snprint")]
        public unowned string bytes_snprint(char[] buf, uint bytes);

        public string bytes_sprint(uint bytes) {
                var buffer = new char[BYTES_SNPRINT_MAX];
                bytes_snprint(buffer, bytes);
                return (string) buffer;
        }

        [CCode (cname="pa_volume_t", has_type_id=false)]
        public struct Volume : uint32 {

                [CCode (cname="PA_SW_VOLUME_SNPRINT_DB_MAX")]
                public const size_t SW_SNPRINT_DB_MAX;

                [CCode (cname="PA_VOLUME_SNPRINT_MAX")]
                public const size_t SNPRINT_MAX;

                [CCode (cname="PA_VOLUME_MAX")]
                public const Volume MAX;

                [CCode (cname="PA_VOLUME_NORM")]
                public const Volume NORM;

                [CCode (cname="PA_VOLUME_MUTED")]
                public const Volume MUTED;

                // [CCode (cname="PA_VOLUME_INVALID")]
                [CCode (cname="PA_VOLUME_MAX")]
                public const Volume INVALID;

                [CCode (cname="pa_volume_snprint", instance_pos = 3.1)]
                public unowned string snprint(char[] s);

                public string sprint() {
                        var buffer = new char[SNPRINT_MAX];
                        this.snprint(buffer);
                        return (string) buffer;
                }

                public string to_string() {
                        return sprint();
                }

                [CCode (cname="pa_sw_volume_snprint_dB", instance_pos = 3.1)]
                public unowned string sw_snprint_dB(char[] s);

                public string sw_sprint_dB() {
                        var buffer = new char[SW_SNPRINT_DB_MAX];
                        this.sw_snprint_dB(buffer);
                        return (string) buffer;
                }

                [CCode (cname="pa_sw_volume_multiply")]
                public Volume sw_multiply(Volume other);

                [CCode (cname="pa_sw_volume_divide")]
                public Volume sw_divide(Volume other);

                [CCode (cname="pa_sw_volume_from_dB")]
                public static Volume sw_from_dB(double f);

                [CCode (cname="pa_sw_volume_to_dB")]
                public double sw_to_dB();

                [CCode (cname="pa_sw_volume_from_linear")]
                public static Volume sw_from_linear(double f);

                [CCode (cname="pa_sw_volume_to_linear")]
                public double sw_to_linear();
        }

        [CCode (cname="PA_DECIBEL_MININFTY")]
        public const double DECIBEL_MININFTY;

        [CCode (cname="PA_CHANNELS_MAX")]
        public const int CHANNELS_MAX;

        [CCode (cname="PA_RATE_MAX")]
        public const int RATE_MAX;

        [CCode (cname="pa_cvolume", has_type_id=false)]
        public struct CVolume {
                public uint8 channels;
                public Volume values[CHANNELS_MAX];

                [CCode (cname="PA_SW_CVOLUME_SNPRINT_DB_MAX")]
                public const size_t SW_SNPRINT_DB_MAX;

                [CCode (cname="PA_CVOLUME_SNPRINT_MAX")]
                public const size_t SNPRINT_MAX;

                [CCode (cname="pa_cvolume_equal")]
                public bool equal(CVolume other);

                [CCode (cname="pa_cvolume_init")]
                public unowned CVolume? init();

                [CCode (cname="pa_cvolume_reset")]
                public unowned CVolume? reset(uint8 channels);

                [CCode (cname="pa_cvolume_mute")]
                public unowned CVolume? mute(uint8 channels);

                [CCode (cname="pa_cvolume_snprint", instance_pos = 3.1)]
                public unowned string snprint(char[] s);

                public string sprint() {
                        var buffer = new char[SNPRINT_MAX];
                        this.snprint(buffer);
                        return (string) buffer;
                }

                public string to_string() {
                        return sprint();
                }

                [CCode (cname="pa_sw_cvolume_snprint_dB", instance_pos = 3.1)]
                public unowned string sw_snprint_dB(char [] s);

                public string sw_sprint_dB() {
                        var buffer = new char[SW_SNPRINT_DB_MAX];
                        this.sw_snprint_dB(buffer);
                        return (string) buffer;
                }

                [CCode (cname="pa_cvolume_init")]
                public CVolume();

                [CCode (cname="pa_cvolume_avg")]
                public Volume avg();

                [CCode (cname="pa_cvolume_max")]
                public Volume max();

                [CCode (cname="pa_cvolume_min")]
                public Volume min();

                [CCode (cname="pa_cvolume_avg_mask")]
                public Volume avg_mask(ChannelMap map, ChannelPositionMask mask);

                [CCode (cname="pa_cvolume_max_mask")]
                public Volume max_mask(ChannelMap map, ChannelPositionMask mask);

                [CCode (cname="pa_cvolume_min_mask")]
                public Volume min_mask(ChannelMap map, ChannelPositionMask mask);

                [CCode (cname="pa_cvolume_valid")]
                public bool valid();

                [CCode (cname="pa_cvolume_channels_equal_to")]
                public bool channels_equal_to(Volume other);

                [CCode (cname="pa_cvolume_is_muted")]
                public bool is_muted();

                [CCode (cname="pa_cvolume_is_norm")]
                public bool is_norm();

                [CCode (cname="pa_sw_cvolume_multiply")]
                public unowned CVolume? multiply(CVolume other);

                [CCode (cname="pa_sw_cvolume_divide")]
                public unowned CVolume? divide(CVolume other);

                [CCode (cname="pa_sw_cvolume_multiply_scalar")]
                public unowned CVolume? multiply_scalar(Volume other);

                [CCode (cname="pa_sw_cvolume_divide_scalar")]
                public unowned CVolume? divide_scalar(Volume other);

                [CCode (cname="pa_cvolume_remap")]
                public unowned CVolume? remap(ChannelMap from, ChannelMap to);

                [CCode (cname="pa_cvolume_compatible")]
                public bool compatible(SampleSpec ss);

                [CCode (cname="pa_cvolume_compatible_with_channel_map")]
                public bool compatible_with_channel_map(ChannelMap cm);

                [CCode (cname="pa_cvolume_set")]
                public unowned CVolume? set(uint8 channels, Volume v);

                [CCode (cname="pa_cvolume_get_balance")]
                public float get_balance(ChannelMap map);

                [CCode (cname="pa_cvolume_set_balance")]
                public unowned CVolume? set_balance(ChannelMap map, float b);

                [CCode (cname="pa_cvolume_get_fade")]
                public float get_fade(ChannelMap map);

                [CCode (cname="pa_cvolume_set_fade")]
                public unowned CVolume? set_fade(ChannelMap map, float f);

                [CCode (cname="pa_cvolume_scale")]
                public unowned CVolume? scale(Volume max);

                [CCode (cname="pa_cvolume_scale_mask")]
                public unowned CVolume? scale_mask(Volume max, ChannelMap map, ChannelPositionMask mask);

                [CCode (cname="pa_cvolume_set_position")]
                public unowned CVolume? set_position(ChannelMap map, ChannelPosition p, Volume v);

                [CCode (cname="pa_cvolume_get_position")]
                public Volume get_position(ChannelMap map, ChannelPosition p);

                [CCode (cname="pa_cvolume_merge")]
                public unowned CVolume? merge(CVolume other);

                [CCode (cname="pa_cvolume_inc")]
                public unowned CVolume? inc(Volume plus = 1);

                [CCode (cname="pa_cvolume_dec")]
                public unowned CVolume? dec(Volume minus = 1);
        }

        [CCode (cname="pa_channel_map", has_type_id=false)]
        public struct ChannelMap {
                public uint8 channels;
                // TODO: Replace array length with CHANNELS_MAX once vala's bug #647788 is fixed
                public ChannelPosition map[32];

                [CCode (cname="PA_CHANNEL_MAP_SNPRINT_MAX")]
                public const size_t SNPRINT_MAX;

                [CCode (cname="pa_channel_map_init")]
                public ChannelMap();

                [CCode (cname="pa_channel_map_init")]
                public unowned ChannelMap? init();

                [CCode (cname="pa_channel_map_init_mono")]
                public unowned ChannelMap? init_mono();

                [CCode (cname="pa_channel_map_init_stereo")]
                public unowned ChannelMap? init_stereo();

                [CCode (cname="pa_channel_map_init_auto")]
                public unowned ChannelMap? init_auto(uint8 channels, ChannelMapDef def = ChannelMapDef.DEFAULT);

                [CCode (cname="pa_channel_map_init_extend")]
                public unowned ChannelMap? init_extend(uint8 channels, ChannelMapDef def = ChannelMapDef.DEFAULT);

                [CCode (cname="pa_channel_map_snprint", instance_pos = 3.1)]
                public unowned string snprint(char[] s);

                public string sprint() {
                        var buffer = new char[SNPRINT_MAX];
                        this.snprint(buffer);
                        return (string) buffer;
                }

                public string to_string() {
                        return sprint();
                }

                [CCode (cname="pa_channel_map_parse")]
                public unowned ChannelMap? parse(string s);

                [CCode (cname="pa_channel_map_equal")]
                public bool equal(ChannelMap other);

                [CCode (cname="pa_channel_map_superset")]
                public bool superset(ChannelMap other);

                [CCode (cname="pa_channel_map_valid")]
                public bool valid();

                [CCode (cname="pa_channel_map_compatible")]
                public bool compatible(SampleSpec ss);

                [CCode (cname="pa_channel_map_can_balance")]
                public bool can_balance();

                [CCode (cname="pa_channel_map_can_fade")]
                public bool can_fade();

                [CCode (cname="pa_channel_map_to_name")]
                public unowned string? to_name();

                [CCode (cname="pa_channel_map_to_pretty_name")]
                public unowned string? to_pretty_name();

                [CCode (cname="pa_channel_map_has_position")]
                public bool has_position(ChannelPosition p);

                [CCode (cname="pa_channel_map_mask")]
                public ChannelPositionMask mask();
        }

        [CCode (cname="pa_channel_position_mask_t", has_type_id=false)]
        public struct ChannelPositionMask : uint64 {
        }

        [CCode (cname="pa_channel_position_t", cprefix="PA_CHANNEL_POSITION_", has_type_id=false)]
        public enum ChannelPosition {
                INVALID,
                MONO,
                FRONT_LEFT,
                FRONT_RIGHT,
                FRONT_CENTER,
                REAR_CENTER,
                REAR_LEFT,
                REAR_RIGHT,
                LFE,
                FRONT_LEFT_OF_CENTER,
                FRONT_RIGHT_OF_CENTER,
                SIDE_LEFT,
                SIDE_RIGHT,
                TOP_CENTER,
                AUX0,
                AUX1,
                AUX2,
                AUX3,
                AUX4,
                AUX5,
                AUX6,
                AUX7,
                AUX8,
                AUX9,
                AUX10,
                AUX11,
                AUX12,
                AUX13,
                AUX14,
                AUX15,
                AUX16,
                AUX17,
                AUX18,
                AUX19,
                AUX20,
                AUX21,
                AUX22,
                AUX23,
                AUX24,
                AUX25,
                AUX26,
                AUX27,
                AUX28,
                AUX29,
                AUX30,
                AUX31,
                MAX;

                [CCode (cname="PA_CHANNEL_POSITION_MASK")]
                public ChannelPositionMask mask();

                [CCode (cname="pa_channel_position_to_string")]
                public unowned string? to_string();

                [CCode (cname="pa_channel_position_to_pretty_string")]
                public unowned string? to_pretty_string();

                [CCode (cname="pa_channel_position_from_string")]
                public static ChannelPosition from_string(string s);
        }

        [CCode (cname="pa_channel_map_def_t", cprefix="PA_CHANNEL_MAP_", has_type_id=false)]
        public enum ChannelMapDef {
                AIFF,
                WAVEEX,
                AUX,
                DEFAULT,

                [CCode (cname="PA_CHANNEL_MAP_DEF_MAX")]
                MAX
        }

        [Compact]
        [CCode (cname="pa_proplist", cprefix="pa_proplist_", free_function="pa_proplist_free", has_type_id=false)]
        public class Proplist {

                [CCode (cname="PA_PROP_MEDIA_NAME")]
                public const string PROP_MEDIA_NAME;
                [CCode (cname="PA_PROP_MEDIA_TITLE")]
                public const string PROP_MEDIA_TITLE;
                [CCode (cname="PA_PROP_MEDIA_ARTIST")]
                public const string PROP_MEDIA_ARTIST;
                [CCode (cname="PA_PROP_MEDIA_COPYRIGHT")]
                public const string PROP_MEDIA_COPYRIGHT;
                [CCode (cname="PA_PROP_MEDIA_SOFTWARE")]
                public const string PROP_MEDIA_SOFTWARE;
                [CCode (cname="PA_PROP_MEDIA_LANGUAGE")]
                public const string PROP_MEDIA_LANGUAGE;
                [CCode (cname="PA_PROP_MEDIA_FILENAME")]
                public const string PROP_MEDIA_FILENAME;
                [CCode (cname="PA_PROP_MEDIA_ICON_NAME")]
                public const string PROP_MEDIA_ICON_NAME;
                [CCode (cname="PA_PROP_MEDIA_ROLE")]
                public const string PROP_MEDIA_ROLE;
                [CCode (cname="PA_PROP_EVENT_ID")]
                public const string PROP_EVENT_ID;
                [CCode (cname="PA_PROP_EVENT_DESCRIPTION")]
                public const string PROP_EVENT_DESCRIPTION;
                [CCode (cname="PA_PROP_EVENT_MOUSE_X")]
                public const string PROP_EVENT_MOUSE_X;
                [CCode (cname="PA_PROP_EVENT_MOUSE_Y")]
                public const string PROP_EVENT_MOUSE_Y;
                [CCode (cname="PA_PROP_EVENT_MOUSE_HPOS")]
                public const string PROP_EVENT_MOUSE_HPOS;
                [CCode (cname="PA_PROP_EVENT_MOUSE_VPOS")]
                public const string PROP_EVENT_MOUSE_VPOS;
                [CCode (cname="PA_PROP_EVENT_MOUSE_BUTTON")]
                public const string PROP_EVENT_MOUSE_BUTTON;
                [CCode (cname="PA_PROP_WINDOW_NAME")]
                public const string PROP_WINDOW_NAME;
                [CCode (cname="PA_PROP_WINDOW_ID")]
                public const string PROP_WINDOW_ID;
                [CCode (cname="PA_PROP_WINDOW_ICON_NAME")]
                public const string PROP_WINDOW_ICON_NAME;
                [CCode (cname="PA_PROP_WINDOW_X11_DISPLAY")]
                public const string PROP_WINDOW_X11_DISPLAY;
                [CCode (cname="PA_PROP_WINDOW_X11_SCREEN")]
                public const string PROP_WINDOW_X11_SCREEN;
                [CCode (cname="PA_PROP_WINDOW_X11_MONITOR")]
                public const string PROP_WINDOW_X11_MONITOR;
                [CCode (cname="PA_PROP_WINDOW_X11_XID")]
                public const string PROP_WINDOW_X11_XID;
                [CCode (cname="PA_PROP_APPLICATION_NAME")]
                public const string PROP_APPLICATION_NAME;
                [CCode (cname="PA_PROP_APPLICATION_ID")]
                public const string PROP_APPLICATION_ID;
                [CCode (cname="PA_PROP_APPLICATION_VERSION")]
                public const string PROP_APPLICATION_VERSION;
                [CCode (cname="PA_PROP_APPLICATION_ICON_NAME")]
                public const string PROP_APPLICATION_ICON_NAME;
                [CCode (cname="PA_PROP_APPLICATION_LANGUAGE")]
                public const string PROP_APPLICATION_LANGUAGE;
                [CCode (cname="PA_PROP_APPLICATION_PROCESS_ID")]
                public const string PROP_APPLICATION_PROCESS_ID;
                [CCode (cname="PA_PROP_APPLICATION_PROCESS_BINARY")]
                public const string PROP_APPLICATION_PROCESS_BINARY;
                [CCode (cname="PA_PROP_APPLICATION_PROCESS_USER")]
                public const string PROP_APPLICATION_PROCESS_USER;
                [CCode (cname="PA_PROP_APPLICATION_PROCESS_HOST")]
                public const string PROP_APPLICATION_PROCESS_HOST;
                [CCode (cname="PA_PROP_APPLICATION_PROCESS_MACHINE_ID")]
                public const string PROP_APPLICATION_PROCESS_MACHINE_ID;
                [CCode (cname="PA_PROP_APPLICATION_PROCESS_SESSION_ID")]
                public const string PROP_APPLICATION_PROCESS_SESSION_ID;
                [CCode (cname="PA_PROP_DEVICE_STRING")]
                public const string PROP_DEVICE_STRING;
                [CCode (cname="PA_PROP_DEVICE_API")]
                public const string PROP_DEVICE_API;
                [CCode (cname="PA_PROP_DEVICE_DESCRIPTION")]
                public const string PROP_DEVICE_DESCRIPTION;
                [CCode (cname="PA_PROP_DEVICE_BUS_PATH")]
                public const string PROP_DEVICE_BUS_PATH;
                [CCode (cname="PA_PROP_DEVICE_SERIAL")]
                public const string PROP_DEVICE_SERIAL;
                [CCode (cname="PA_PROP_DEVICE_VENDOR_ID")]
                public const string PROP_DEVICE_VENDOR_ID;
                [CCode (cname="PA_PROP_DEVICE_VENDOR_NAME")]
                public const string PROP_DEVICE_VENDOR_NAME;
                [CCode (cname="PA_PROP_DEVICE_PRODUCT_ID")]
                public const string PROP_DEVICE_PRODUCT_ID;
                [CCode (cname="PA_PROP_DEVICE_PRODUCT_NAME")]
                public const string PROP_DEVICE_PRODUCT_NAME;
                [CCode (cname="PA_PROP_DEVICE_CLASS")]
                public const string PROP_DEVICE_CLASS;
                [CCode (cname="PA_PROP_DEVICE_FORM_FACTOR")]
                public const string PROP_DEVICE_FORM_FACTOR;
                [CCode (cname="PA_PROP_DEVICE_BUS")]
                public const string PROP_DEVICE_BUS;
                [CCode (cname="PA_PROP_DEVICE_ICON_NAME")]
                public const string PROP_DEVICE_ICON_NAME;
                [CCode (cname="PA_PROP_DEVICE_ACCESS_MODE")]
                public const string PROP_DEVICE_ACCESS_MODE;
                [CCode (cname="PA_PROP_DEVICE_MASTER_DEVICE")]
                public const string PROP_DEVICE_MASTER_DEVICE;
                [CCode (cname="PA_PROP_DEVICE_BUFFERING_BUFFER_SIZE")]
                public const string PROP_DEVICE_BUFFERING_BUFFER_SIZE;
                [CCode (cname="PA_PROP_DEVICE_BUFFERING_FRAGMENT_SIZE")]
                public const string PROP_DEVICE_BUFFERING_FRAGMENT_SIZE;
                [CCode (cname="PA_PROP_DEVICE_PROFILE_NAME")]
                public const string PROP_DEVICE_PROFILE_NAME;
                [CCode (cname="PA_PROP_DEVICE_INTENDED_ROLES")]
                public const string PROP_DEVICE_INTENDED_ROLES;
                [CCode (cname="PA_PROP_DEVICE_PROFILE_DESCRIPTION")]
                public const string PROP_DEVICE_PROFILE_DESCRIPTION;
                [CCode (cname="PA_PROP_MODULE_AUTHOR")]
                public const string PROP_MODULE_AUTHOR;
                [CCode (cname="PA_PROP_MODULE_DESCRIPTION")]
                public const string PROP_MODULE_DESCRIPTION;
                [CCode (cname="PA_PROP_MODULE_USAGE")]
                public const string PROP_MODULE_USAGE;
                [CCode (cname="PA_PROP_MODULE_VERSION")]
                public const string PROP_MODULE_VERSION;

                [CCode (cname="pa_proplist_new")]
                public Proplist();

                public int sets(string key, string value);
                public int setp(string pair);

                [PrintfFormat]
                public int setf(string key, string format, ...);

                public int set(string key, void* data, size_t size);

                public unowned string? gets(string key);

                public int get(string key, out void* data, out size_t size);

                public void update(UpdateMode mode, Proplist other);

                public void unset(string key);

                public void unset_many([CCode (array_length = false)] string[] key);

                public unowned string? iterate(ref void* state);

                public string to_string();

                public string to_string_sep(string sep);

                public static Proplist? from_string(string s);

                public int contains(string key);

                public void clear();

                public Proplist copy();

                public uint size();

                [CCode (cname="pa_proplist_isempty")]
                public bool is_empty();
        }

        [CCode (cname="pa_update_mode_t", cprefix="PA_UPDATE_", has_type_id=false)]
        public enum UpdateMode {
                SET,
                MERGE,
                REPLACE
        }

        [CCode (cname="PA_OK")]
        public const int OK;

        [CCode (cname="int", cprefix="PA_ERR_", has_type_id=false)]
        public enum Error {
                ACCESS,
                COMMAND,
                INVALID,
                EXIST,
                NOENTITY,
                CONNECTIONREFUSED,
                PROTOCOL,
                TIMEOUT,
                AUTHKEY,
                INTERNAL,
                CONNECTIONTERMINATED,
                KILLED,
                INVALIDSERVER,
                MODINITFAILED,
                BADSTATE,
                NODATA,
                VERSION,
                TOOLARGE,
                NOTSUPPORTED,
                UNKNOWN,
                NOEXTENSION,
                OBSOLETE,
                NOTIMPLEMENTED,
                FORKED,
                IO,
                MAX
        }

        [CCode (cname="pa_strerror")]
        public unowned string? strerror(Error e);

        [CCode (has_target = false, has_typedef = false)]
        public delegate void VoidFunc();

        [CCode (cname="pa_spawn_api", has_type_id=false)]
        public struct SpawnApi {
                VoidFunc? prefork;
                VoidFunc? postfork;
                VoidFunc? atfork;
        }

        [CCode (cname="pa_io_event_flags_t", cprefix="PA_IO_EVENT_", has_type_id=false)]
        public enum IoEventFlags {
                NULL,
                INPUT,
                OUTPUT,
                HANGUP,
                ERROR
        }

        [Compact]
        [CCode (cname="pa_io_event", unref_function="", ref_function="", has_type_id=false)]
        public class IoEvent {
        }

        [Compact]
        [CCode (cname="pa_time_event", unref_function="", ref_function="", has_type_id=false)]
        public class TimeEvent {
        }

        [Compact]
        [CCode (cname="pa_defer_event", unref_function="", ref_function="", has_type_id=false)]
        public class DeferEvent {
        }

        [Compact]
        [CCode (cname="pa_signal_event", cprefix="pa_signal_", free_function="pa_signal_free", has_type_id=false)]
        public class SignalEvent {

                [CCode (cname="pa_signal_new")]
                public SignalEvent(int sig, MainLoopApi.SignalEventCb cb);

                public void set_destroy([CCode (delegate_target = false)] MainLoopApi.SignalEventDestroyCb cb);
        }

        [Compact]
        [CCode (cname="pa_mainloop_api", unref_function="", ref_function="", has_type_id=false)]
        public class MainLoopApi {
                public void* userdata;

                /* Callbacks for the consumer to implement*/
                [CCode (cname = "pa_io_event_cb_t")]
                public delegate void IoEventCb(MainLoopApi a, IoEvent e, int fd, IoEventFlags flags);
                [CCode (cname = "pa_io_event_destroy_cb_t")]
                public delegate void IoEventDestroyCb(MainLoopApi a, IoEvent e);

                [CCode (cname = "pa_time_event_cb_t")]
                public delegate void TimeEventCb(MainLoopApi a, TimeEvent e, ref timeval t);
                [CCode (cname = "pa_time_event_destroy_cb_t")]
                public delegate void TimeEventDestroyCb(MainLoopApi a, TimeEvent e);

                [CCode (cname = "pa_defer_event_cb_t")]
                public delegate void DeferEventCb(MainLoopApi a, DeferEvent e);
                [CCode (cname = "pa_defer_event_destroy_cb_t")]
                public delegate void DeferEventDestroyCb(MainLoopApi a, DeferEvent e);

                [CCode (cname = "pa_signal_cb_t")]
                public delegate void SignalEventCb(MainLoopApi a, SignalEvent e, int sig);
                [CCode (cname = "pa_signal_destroy_cb_t")]
                public delegate void SignalEventDestroyCb(MainLoopApi a, SignalEvent e);

                /* Callbacks for the provider to implement */
                [CCode (has_target = false, has_typedef = false)]
                public delegate IoEvent IoNewCb(MainLoopApi a, int fd, IoEventFlags flags, IoEventCb cb);
                [CCode (has_target = false, has_typedef = false)]
                public delegate void IoEnableCb(IoEvent e, IoEventFlags flags);
                [CCode (has_target = false, has_typedef = false)]
                public delegate void IoFreeCb(IoEvent e);
                [CCode (has_target = false, has_typedef = false)]
                public delegate void IoSetDestroyCb(IoEvent e, [CCode (delegate_target = false)] IoEventDestroyCb? cb);

                [CCode (has_target = false, has_typedef = false)]
                public delegate TimeEvent TimeNewCb(MainLoopApi a, timeval? t, TimeEventCb cb);
                [CCode (has_target = false, has_typedef = false)]
                public delegate void TimeRestartCb(TimeEvent e, timeval? t);
                [CCode (has_target = false, has_typedef = false)]
                public delegate void TimeFreeCb(TimeEvent e);
                [CCode (has_target = false, has_typedef = false)]
                public delegate void TimeSetDestroyCb(TimeEvent e, [CCode (delegate_target = false)] TimeEventDestroyCb? cb);

                [CCode (has_target = false, has_typedef = false)]
                public delegate DeferEvent DeferNewCb(MainLoopApi a, DeferEventCb cb);
                [CCode (has_target = false, has_typedef = false)]
                public delegate void DeferEnableCb(DeferEvent e, bool b);
                [CCode (has_target = false, has_typedef = false)]
                public delegate void DeferFreeCb(DeferEvent e);
                [CCode (has_target = false, has_typedef = false)]
                public delegate void DeferSetDestroyCb(DeferEvent e, [CCode (delegate_target = false)] DeferEventDestroyCb? cb);

                [CCode (has_target = false, has_typedef = false)]
                public delegate void QuitCb(MainLoopApi a, int retval);

                [CCode (has_typedef = false)]
                public delegate void OnceCb(MainLoopApi a);

                public IoNewCb io_new;
                public IoEnableCb io_enable;
                public IoFreeCb io_free;
                public IoSetDestroyCb io_set_destroy;

                public TimeNewCb time_new;
                public TimeRestartCb time_restart;
                public TimeFreeCb time_free;
                public TimeSetDestroyCb time_set_destroy;

                public DeferNewCb defer_new;
                public DeferEnableCb defer_enable;
                public DeferFreeCb defer_free;
                public DeferSetDestroyCb defer_set_destroy;

                public QuitCb quit;

                [CCode (cname="pa_mainloop_api_once")]
                public void once(OnceCb cb);
        }

        [CCode (cname="pa_signal_init")]
        public void signal_init(MainLoopApi api);

        [CCode (cname="pa_signal_done")]
        public void signal_done();

        [CCode (cname="pa_poll_func")]
        public delegate int PollFunc(pollfd[] ufds, int timeout);

        [Compact]
        [CCode (cname="pa_mainloop", cprefix="pa_mainloop_", free_function="pa_mainloop_free", has_type_id=false)]
        public class MainLoop {

                [CCode (cname="pa_mainloop_new")]
                public MainLoop();

                public int prepare(int timeout = -1);
                public int poll();
                public int dispatch();
                public int get_retval();
                public int iterate(bool block = true, out int retval = null);
                public int run(out int retval = null);
                public unowned MainLoopApi get_api();
                public void quit(int retval);
                public void wakeup();
                public void set_poll_func(PollFunc poll_func);
        }

        [Compact]
        [CCode (cname="pa_threaded_mainloop", cprefix="pa_threaded_mainloop_", free_function="pa_threaded_mainloop_free")]
        public class ThreadedMainLoop {

                [CCode (cname="pa_threaded_mainloop_new")]
                public ThreadedMainLoop();

                public int start();
                public void stop();
                public void lock();
                public void unlock();
                public void wait();
                public void signal(bool WaitForAccept = false);
                public void accept();
                public int get_retval();
                public unowned MainLoopApi get_api();
                public bool in_thread();
        }

        [Compact]
        [CCode (cname="pa_operation", cprefix="pa_operation_", unref_function="pa_operation_unref", ref_function="pa_operation_ref", has_type_id=false)]
        public class Operation {

                [CCode (cname="pa_operation_state_t", cprefix="PA_OPERATION_", has_type_id=false)]
                public enum State {
                        RUNNING,
                        DONE,
                        CANCELED
                }

                public void cancel();
                public State get_state();
        }

        [Compact]
        [CCode (cname="pa_context", cprefix="pa_context_", unref_function="pa_context_unref", ref_function="pa_context_ref", has_type_id=false)]
        public class Context {

                [CCode (cname="pa_context_flags_t", cprefix="PA_CONTEXT_", has_type_id=false)]
                public enum Flags {
                        NOAUTOSPAWN,
                        NOFAIL
                }

                [CCode (cname="pa_context_state_t", cprefix="PA_CONTEXT_", has_type_id=false)]
                public enum State {
                        UNCONNECTED,
                        CONNECTING,
                        AUTHORIZING,
                        SETTING_NAME,
                        READY,
                        FAILED,
                        TERMINATED;

                        [CCode (cname="PA_CONTEXT_IS_GOOD", has_type_id=false)]
                        public bool IS_GOOD();
                }

                [CCode (cname="pa_subscription_mask_t", cprefix="PA_SUBSCRIPTION_MASK_", has_type_id=false)]
                public enum SubscriptionMask {
                        NULL,
                        SINK,
                        SOURCE,
                        SINK_INPUT,
                        SOURCE_OUTPUT,
                        MODULE,
                        CLIENT,
                        SAMPLE_CACHE,
                        SERVER,
                        CARD,
                        ALL;

                        [CCode (cname="pa_subscription_match_flags")]
                        public bool match_flags(SubscriptionEventType t);
                }

                [CCode (cname="pa_subscription_event_type_t", cprefix="PA_SUBSCRIPTION_EVENT_", has_type_id=false)]
                public enum SubscriptionEventType {
                        SINK,
                        SOURCE,
                        SINK_INPUT,
                        SOURCE_OUTPUT,
                        MODULE,
                        CLIENT,
                        SAMPLE_CACHE,
                        SERVER,
                        CARD,
                        FACILITY_MASK,
                        NEW,
                        CHANGE,
                        REMOVE,
                        TYPE_MASK
                }

                [CCode (cname = "pa_context_notify_cb_t")]
                public delegate void NotifyCb(Context c);
                [CCode (cname = "pa_context_success_cb_t")]
                public delegate void SuccessCb(Context c, int success);
                [CCode (cname = "pa_context_event_cb_t")]
                public delegate void EventCb(Context c, string name, Proplist? proplist);
                [CCode (cname = "pa_context_subscribe_cb_t")]
                public delegate void SubscribeCb(Context c, SubscriptionEventType t, uint32 idx);
                [CCode (cname = "pa_sink_info_cb_t")]
                public delegate void SinkInfoCb(Context c, SinkInfo? i, int eol);
                [CCode (cname = "pa_source_info_cb_t")]
                public delegate void SourceInfoCb(Context c, SourceInfo? i, int eol);
                [CCode (cname = "pa_card_info_cb_t")]
                public delegate void CardInfoCb(Context c, CardInfo? i, int eol);
                [CCode (cname = "pa_sink_input_info_cb_t")]
                public delegate void SinkInputInfoCb(Context c, SinkInputInfo? i, int eol);
                [CCode (cname = "pa_source_output_info_cb_t")]
                public delegate void SourceOutputInfoCb(Context c, SourceOutputInfo? i, int eol);
                [CCode (cname = "pa_server_info_cb_t")]
                public delegate void ServerInfoCb(Context c, ServerInfo? i);
                [CCode (cname = "pa_stat_info_cb_t")]
                public delegate void StatInfoCb(Context c, ServerInfo? i);
                [CCode (cname = "pa_module_info_cb_t")]
                public delegate void ModuleInfoCb(Context c, ModuleInfo? i, int eol);
                [CCode (cname = "pa_client_info_cb_t")]
                public delegate void ClientInfoCb(Context c, ClientInfo? i, int eol);
                [CCode (cname = "pa_sample_info_cb_t")]
                public delegate void SampleInfoCb(Context c, SampleInfo? i, int eol);
                [CCode (cname = "pa_context_index_cb_t")]
                public delegate void IndexCb(Context c, uint32 idx);

                [CCode (cname="pa_context_new_with_proplist")]
                public Context(MainLoopApi api, string? name, Proplist? proplist = null);

                public void set_state_callback(NotifyCb? cb = null);
                public void set_event_callback(EventCb? cb = null);
                public void set_subscribe_callback(SubscribeCb? cb = null);

                public Error errno();

                public int is_pending();
                public State get_state();
                public int is_local();
                public unowned string? get_server();
                public uint32 get_protocol_version();
                public uint32 get_server_protocol_version();
                public uint32 get_index();

                public int connect(string? server = null, Flags flags = 0, SpawnApi? api = null);
                public void disconnect();

                public Operation? drain(NotifyCb? cb = null);
                public Operation? exit_daemon(SuccessCb? cb = null);
                public Operation? set_default_sink(string name, SuccessCb? cb = null);
                public Operation? set_default_source(string name, SuccessCb? cb = null);
                public Operation? set_name(string name, SuccessCb? cb = null);

                public Operation? proplist_remove([CCode (array_length = false)] string[] keys, SuccessCb? cb = null);
                public Operation? proplist_update(UpdateMode mode, Proplist pl, SuccessCb? cb = null);

                public Operation? subscribe(SubscriptionMask mask, SuccessCb? cb = null);

                public Operation? get_sink_info_by_name(string name, SinkInfoCb cb);
                public Operation? get_sink_info_by_index(uint32 idx, SinkInfoCb cb);
                public Operation? get_sink_info_list(SinkInfoCb cb);

                public Operation? set_sink_volume_by_name(string name, CVolume volume, SuccessCb? cb = null);
                public Operation? set_sink_volume_by_index(uint32 idx, CVolume volume, SuccessCb? cb = null);
                public Operation? set_sink_mute_by_name(string name, bool mute, SuccessCb? cb = null);
                public Operation? set_sink_mute_by_index(uint32 idx, bool mute, SuccessCb? cb = null);

                public Operation? suspend_sink_by_name(string name, bool suspend, SuccessCb? cb = null);
                public Operation? suspend_sink_by_index(uint32 idx, bool suspend, SuccessCb? cb = null);

                public Operation? set_sink_port_by_name(string name, string port, SuccessCb? cb = null);
                public Operation? set_sink_port_by_index(uint32 idx, string port, SuccessCb? cb = null);

                public Operation? get_source_info_by_name(string name, SourceInfoCb cb);
                public Operation? get_source_info_by_index(uint32 idx, SourceInfoCb cb);
                public Operation? get_source_info_list(SourceInfoCb cb);

                public Operation? set_source_volume_by_name(string name, CVolume volume, SuccessCb? cb = null);
                public Operation? set_source_volume_by_index(uint32 idx, CVolume volume, SuccessCb? cb = null);
                public Operation? set_source_mute_by_name(string name, bool mute, SuccessCb? cb = null);
                public Operation? set_source_mute_by_index(uint32 idx, bool mute, SuccessCb? cb = null);

                public Operation? suspend_source_by_name(string name, bool suspend, SuccessCb? cb = null);
                public Operation? suspend_source_by_index(uint32 idx, bool suspend, SuccessCb? cb = null);

                public Operation? set_source_port_by_name(string name, string port, SuccessCb? cb = null);
                public Operation? set_source_port_by_index(uint32 idx, string port, SuccessCb? cb = null);

                public Operation? get_server_info(ServerInfoCb cb);

                public Operation? get_module_info(uint32 idx, ModuleInfoCb cb);
                public Operation? get_module_info_list(ModuleInfoCb cb);

                public Operation? load_module(string name, string? argument, IndexCb? cb = null);
                public Operation? unload_module(uint32 idx, SuccessCb? cb = null);

                public Operation? get_client_info(uint32 idx, ClientInfoCb cb);
                public Operation? get_client_info_list(ClientInfoCb cb);

                public Operation? kill_client(uint32 idx, SuccessCb? cb = null);

                public Operation? get_card_info_by_name(string name, CardInfoCb cb);
                public Operation? get_card_info_by_index(uint32 idx, CardInfoCb cb);
                public Operation? get_card_info_list(CardInfoCb cb);

                public Operation? set_card_profile_by_index(uint32 idx, string profile, SuccessCb? cb = null);
                public Operation? set_card_profile_by_name(string name, string profile, SuccessCb? cb = null);

                public Operation? get_sink_input_info(uint32 idx, SinkInputInfoCb cb);
                public Operation? get_sink_input_info_list(SinkInputInfoCb cb);

                public Operation? move_sink_input_by_index(uint32 idx, uint32 sink_idx, SuccessCb? cb = null);
                public Operation? move_sink_input_by_name(uint32 idx, string sink_name, SuccessCb? cb = null);

                public Operation? set_sink_input_volume(uint32 idx, CVolume volume, SuccessCb? cb = null);
                public Operation? set_sink_input_mute(uint32 idx, bool mute, SuccessCb? cb = null);

                public Operation? set_source_output_volume(uint32 idx, CVolume volume, SuccessCb? cb = null);
                public Operation? set_source_output_mute(uint32 idx, bool mute, SuccessCb? cb = null);

                public Operation? kill_sink_input(uint32 idx, SuccessCb? cb = null);

                public Operation? get_source_output_info(uint32 idx, SourceOutputInfoCb cb);
                public Operation? get_source_output_info_list(SourceOutputInfoCb cb);

                public Operation? move_source_output_by_index(uint32 idx, uint32 source_idx, SuccessCb? cb = null);
                public Operation? move_source_output_by_name(uint32 idx, string source_name, SuccessCb? cb = null);

                public Operation? kill_source_output(uint32 idx, SuccessCb? cb = null);

                public Operation? stat(StatInfoCb cb);

                public Operation? get_sample_info_by_name(string name, SampleInfoCb cb);
                public Operation? get_sample_info_by_index(uint32 idx, SampleInfoCb cb);
                public Operation? get_sample_info_list(SampleInfoCb cb);

                public Operation? remove_sample(string name, SuccessCb? cb = null);

                public Operation? play_sample(string name, string? device = null, Volume volume = Volume.INVALID, SuccessCb? cb = null);
                public Operation? play_sample_with_proplist(string name, string? device = null, Volume volume = Volume.INVALID, Proplist? p = null, IndexCb? cb = null);
        }

        [Compact]
        [CCode (cname="pa_stream", cprefix="pa_stream_", unref_function="pa_stream_unref", ref_function="pa_stream_ref", has_type_id=false)]
        public class Stream {

                [CCode (cname="pa_stream_flags_t", cprefix="PA_STREAM_", has_type_id=false)]
                public enum Flags {
                        START_CORKED,
                        INTERPOLATE_TIMING,
                        NOT_MONOTONIC,
                        AUTO_TIMING_UPDATE,
                        NO_REMAP_CHANNELS,
                        NO_REMIX_CHANNELS,
                        FIX_FORMAT,
                        FIX_RATE,
                        FIX_CHANNELS,
                        DONT_MOVE,
                        VARIABLE_RATE,
                        PEAK_DETECT,
                        START_MUTED,
                        ADJUST_LATENCY,
                        EARLY_REQUESTS,
                        DONT_INHIBIT_AUTO_SUSPEND,
                        START_UNMUTED,
                        FAIL_ON_SUSPEND
                }

                [CCode (cname="pa_stream_state_t", cprefix="PA_STREAM_", has_type_id=false)]
                public enum State {
                        UNCONNECTED,
                        CREATING,
                        READY,
                        FAILED,
                        TERMINATED;

                        [CCode (cname="PA_STREAM_IS_GOOD")]
                        public bool IS_GOOD();
                }

                [CCode (cname="pa_stream_direction_t", cprefix="PA_STREAM_", has_type_id=false)]
                public enum Direction {
                        NODIRECTION,
                        PLAYBACK,
                        RECORD,
                        UPLOAD
                }

                [CCode (cname="pa_seek_mode_t", cprefix="PA_SEEK_", has_type_id=false)]
                public enum SeekMode {
                        RELATIVE,
                        ABSOLUTE,
                        RELATIVE_ON_READ,
                        RELATIVE_END
                }

                [CCode (cname="pa_buffer_attr", has_type_id=false)]
                public struct BufferAttr {
                        uint32 maxlength;
                        uint32 tlength;
                        uint32 prebuf;
                        uint32 minreq;
                        uint32 fragsize;
                }

                [CCode (cname="pa_timing_info", has_type_id=false)]
                public struct TimingInfo {
                        timeval timestamp;
                        int synchronized_clocks;
                        usec sink_usec;
                        usec source_usec;
                        usec transport_usec;
                        int playing;
                        int write_index_corrupt;
                        int64 write_index;
                        int read_index_corrupt;
                        int64 read_index;
                        usec configured_sink_usec;
                        usec configured_source_usec;
                        int64 since_underrun;
                }

                [CCode (cname="PA_STREAM_EVENT_REQUEST_CORK")]
                public const string EVENT_REQUEST_CORK;

                [CCode (cname="PA_STREAM_EVENT_REQUEST_UNCORK")]
                public const string EVENT_REQUEST_UNCORK;

                [CCode (cname = "pa_stream_success_cb_t")]
                public delegate void SuccessCb(Stream s, int success);
                [CCode (cname = "pa_stream_request_cb_t")]
                public delegate void RequestCb(Stream s, size_t nbytes);
                [CCode (cname = "pa_stream_notify_cb_t")]
                public delegate void NotifyCb(Stream s);
                [CCode (cname = "pa_stream_event_cb_t")]
                public delegate void EventCb(Stream s, string name, Proplist proplist);

                [CCode (cname="pa_stream_new_with_proplist")]
                public Stream(Context c, string name, SampleSpec ss, ChannelMap? map = null, Proplist? proplist = null);

                public State get_state();
                public Context get_context();
                public uint32 get_index();
                public uint32 get_device_index();
                public unowned string? get_device_name();
                public int is_suspended();
                public int is_corked();

                public int connect_playback(string? dev = null, BufferAttr? a = null, Flags flags = 0, CVolume? volume = null, Stream? sync_stream = null);
                public int connect_record(string? dev = null, BufferAttr? a = null, Flags flags = 0);
                public int connect_upload(size_t length);
                public int disconnect();
                public int finish_upload();

                public int begin_write(out void* data, out size_t nbytes);
                public int cancel_write();
                public int write(void *data, size_t bytes, FreeCb? free_cb = null, int64 offset = 0, SeekMode mode = SeekMode.RELATIVE);
                public int peek(out void *data, out size_t nbytes);
                public int drop();
                public size_t writable_size();
                public size_t readable_size();

                public void set_state_callback(NotifyCb? cb = null);
                public void set_write_callback(RequestCb? cb = null);
                public void set_read_callback(RequestCb? cb = null);
                public void set_overflow_callback(NotifyCb? cb = null);
                public void set_underflow_callback(NotifyCb? cb = null);
                public void set_started_callback(NotifyCb? cb = null);
                public void set_latency_update_callback(NotifyCb? cb = null);
                public void set_moved_callback(NotifyCb? cb = null);
                public void set_suspended_callback(NotifyCb? cb = null);
                public void set_event_callback(EventCb? cb = null);
                public void set_buffer_attr_callback(NotifyCb? cb = null);

                public Operation? drain(SuccessCb? cb = null);
                public Operation? update_timing_info(SuccessCb? cb = null);

                public Operation? cork(bool b, SuccessCb? cb = null);
                public Operation? flush(SuccessCb? cb = null);
                public Operation? prebuf(SuccessCb? cb = null);
                public Operation? trigger(SuccessCb? cb = null);

                public Operation? set_name(string name, SuccessCb? cb = null);
                public Operation? set_buffer_attr(BufferAttr attr, SuccessCb? cb = null);
                public Operation? update_sample_rate(uint32 rate, SuccessCb? cb = null);

                public Operation? proplist_remove([CCode (array_length = false)] string[] keys, SuccessCb? cb = null);
                public Operation? proplist_update(UpdateMode mode, Proplist pl, SuccessCb? cb = null);

                public unowned TimingInfo? get_timing_info();
                public int get_time(out usec u);
                public int get_latency(out usec u, out bool negative = null);

                public unowned SampleSpec? get_sample_spec();
                public unowned ChannelMap? get_channel_map();
                public unowned BufferAttr? get_buffer_attr();

                public int set_monitor_stream(uint32 sink_input);
                public uint32 get_monitor_stream();
        }
        
        [CCode (cname="pa_direction_t", cprefix="PA_DIRECTION_", has_type_id=false)]
        [Flags]
        public enum Direction {
                OUTPUT = 0x0001U,
                INPUT = 0x0002U
        }

        [CCode (cname="pa_port_available_t", cprefix="PA_PORT_AVAILABLE_", has_type_id=false)]
        public enum PortAvailable {
                UNKNOWN = 0,
                NO = 1,
                YES = 2
        }

        [CCode (cname="pa_sink_port_info", has_type_id=false)]
        public struct SinkPortInfo {
                public string name;
                public string description;
                public uint32 priority;
        }

        [CCode (cname="pa_sink_info", has_type_id=false)]
        public struct SinkInfo {
                public unowned string name;
                public uint32 index;
                public unowned string description;
                public SampleSpec sample_spec;
                public ChannelMap channel_map;
                public uint32 owner_module;
                public CVolume volume;
                public int mute;
                public uint32 monitor_source;
                public unowned string monitor_source_name;
                public usec latency;
                public unowned string driver;
                public SinkFlags flags;
                public Proplist proplist;
                public usec configured_latency;
                public Volume base_volume;
                public SinkState state;
                public uint32 n_volume_steps;
                public uint32 card;
                public uint32 n_ports;
                [CCode (array_length_cname="n_ports")]
                public SinkPortInfo*[] ports;
                public SinkPortInfo* active_port;
        }

        [CCode (cname="pa_source_port_info", has_type_id=false)]
        public struct SourcePortInfo {
                public unowned string name;
                public unowned string description;
                public uint32 priority;
        }

        [CCode (cname="pa_source_info", has_type_id=false)]
        public struct SourceInfo {
                public unowned string name;
                public uint32 index;
                public unowned string description;
                public SampleSpec sample_spec;
                public ChannelMap channel_map;
                public uint32 owner_module;
                public CVolume volume;
                public int mute;
                public uint32 monitor_of_sink;
                public unowned string monitor_of_sink_name;
                public usec latency;
                public unowned string driver;
                public SourceFlags flags;
                public Proplist proplist;
                public usec configured_latency;
                public Volume base_volume;
                public SourceState state;
                public uint32 n_volume_steps;
                public uint32 card;
                public uint32 n_ports;
                [CCode (array_length_cname="n_ports")]
                public SourcePortInfo*[] ports;
                public SourcePortInfo* active_port;
        }

        [CCode (cname="pa_server_info", has_type_id=false)]
        public struct ServerInfo {
                public unowned string user_name;
                public unowned string host_name;
                public unowned string server_version;
                public unowned string server_name;
                public SampleSpec sample_spec;
                public unowned string default_sink_name;
                public unowned string default_source_name;
                public uint32 cookie;
                public ChannelMap channel_map;
        }

        [CCode (cname="pa_module_info", has_type_id=false)]
        public struct ModuleInfo {
                public uint32 index;
                public unowned string name;
                public unowned string argument;
                public uint32 n_used;
                public Proplist proplist;
        }

        [CCode (cname="pa_client_info", has_type_id=false)]
        public struct ClientInfo {
                public uint32 index;
                public unowned string name;
                public uint32 owner_module;
                public unowned string driver;
                public Proplist proplist;
        }

        [CCode (cname="pa_card_profile_info", has_type_id=false)]
        public struct CardProfileInfo {
                public unowned string name;
                public unowned string description;
                public uint32 n_sinks;
                public uint32 n_sources;
                public uint32 priority;
        }
        
        [CCode (cname="pa_card_profile_info2", has_type_id=false)]
        public struct CardProfileInfo2 {
                public unowned string name;
                public unowned string description;
                public uint32 n_sinks;
                public uint32 n_sources;
                public uint32 priority;
                int available;
        }
        
        [CCode (cname="pa_card_port_info", has_type_id=false)]
        public struct CardPortInfo {
                public unowned string name;
                public unowned string description;
                public uint32 priority;
                PortAvailable available;
                Direction direction;
                public uint32 n_profiles;
                [CCode (array_length_cname="n_profiles")]
                public CardProfileInfo*[] profiles;
                public Proplist proplist;
                public int64 latency_offset;
                [CCode (array_length_cname="n_profiles")]
                public CardProfileInfo2*[] profiles2;
        }

        [CCode (cname="pa_card_info", has_type_id=false)]
        public struct CardInfo {
                public uint32 index;
                public unowned string name;
                public uint32 owner_module;
                public unowned string driver;
                public uint32 n_profiles;
                [CCode (array_length_cname="n_profiles")]
                public CardProfileInfo[] profiles;
                public CardProfileInfo *active_profile;
                [CCode (array_length_cname="n_profiles")]
                public CardProfileInfo2*[] profiles2;
                public CardProfileInfo2 *active_profile2;
                public Proplist proplist;
                public uint32 n_ports;
                [CCode (array_length_cname="n_ports")]
                public CardPortInfo*[] ports;
        }

        [CCode (cname="pa_sink_input_info", has_type_id=false)]
        public struct SinkInputInfo {
                public uint32 index;
                public unowned string name;
                public uint32 owner_module;
                public uint32 client;
                public uint32 sink;
                public SampleSpec sample_spec;
                public ChannelMap channel_map;
                public CVolume volume;
                public uint32 buffer_usec;
                public uint32 sink_usec;
                public unowned string resample_method;
                public unowned string driver;
                public int mute;
                public Proplist proplist;
                public int corked;
                public int has_volume;
                public int volume_writable;
        }

        [CCode (cname="pa_source_output_info", has_type_id=false)]
        public struct SourceOutputInfo {
                public uint32 index;
                public unowned string name;
                public uint32 owner_module;
                public uint32 client;
                public uint32 source;
                public SampleSpec sample_spec;
                public ChannelMap channel_map;
                public uint32 buffer_usec;
                public uint32 source_usec;
                public unowned string resample_method;
                public unowned string driver;
                public Proplist proplist;
                public int corked;
                public CVolume volume;
                public int mute;
                public int has_volume;
                public int volume_writable;
        }

        [CCode (cname="pa_stat_info", has_type_id=false)]
        public struct StatInfo {
                public uint32 memblock_total;
                public uint32 memblock_total_size;
                public uint32 memblock_allocated;
                public uint32 memblock_allocated_size;
                public uint32 scache_size;
        }

        [CCode (cname="pa_sample_info", has_type_id=false)]
        public struct SampleInfo {
                public uint32 index;
                public unowned string name;
                public CVolume volume;
                public SampleSpec sample_spec;
                public ChannelMap channel_map;
                public usec duration;
                public uint32 bytes;
                public bool lazy;
                public unowned string filename;
                public Proplist proplist;
        }

        [CCode (cname="pa_sink_flags_t", cprefix="PA_SINK_", has_type_id=false)]
        public enum SinkFlags {
                HW_VOLUME_CTRL,
                LATENCY,
                HARDWARE,
                NETWORK,
                HW_MUTE_CTRL,
                DECIBEL_VOLUME,
                FLAT_VOLUME,
                DYNAMIC_LATENCY
        }

        [CCode (cname="pa_source_flags_t", cprefix="PA_SOURCE_", has_type_id=false)]
        public enum SourceFlags {
                HW_VOLUME_CTRL,
                LATENCY,
                HARDWARE,
                NETWORK,
                HW_MUTE_CTRL,
                DECIBEL_VOLUME,
                DYNAMIC_LATENCY
        }

        [CCode (cname="pa_sink_state_t", cprefix="PA_SINK_", has_type_id=false)]
        public enum SinkState {
                INVALID_STATE,
                RUNNING,
                IDLE,
                SUSPENDED;

                [CCode (cname="PA_SINK_IS_OPENED")]
                public bool IS_OPENED();
        }

        [CCode (cname="pa_source_state_t", cprefix="PA_SOURCE_", has_type_id=false)]
        public enum SourceState {
                INVALID_STATE,
                RUNNING,
                IDLE,
                SUSPENDED;

                [CCode (cname="PA_SOURCE_IS_OPENED")]
                public bool IS_OPENED();
        }

        [CCode (cname="pa_gettimeofday")]
        public unowned timeval? gettimeofday(out timeval tv);

        [CCode (cname="pa_timeval_diff")]
        public usec timeval_diff(ref timeval a, ref timeval b);

        [CCode (cname="pa_timeval_cmp")]
        public int timeval_cmp(ref timeval a, ref timeval b);

        [CCode (cname="pa_timeval_age")]
        public usec timeval_age(ref timeval a);

        [CCode (cname="pa_timeval_add")]
        public unowned timeval? timeval_add(ref timeval tv, usec x);

        [CCode (cname="pa_timeval_sub")]
        public unowned timeval? timeval_sub(ref timeval tv, usec x);

        [CCode (cname="pa_timeval_store")]
        public unowned timeval? timeval_store(out timeval tv, usec c);

        [CCode (cname="pa_timeval_load")]
        public usec timeval_load(timeval tv);

        [CCode (cname="PA_USEC_MAX")]
        public const usec USEC_MAX;

        [CCode (cname="PA_USEC_INVALID")]
        public const usec USEC_INVALID;

        [CCode (cname="PA_MSEC_PER_SEC")]
        public const usec MSEC_PER_SEC;

        [CCode (cname="PA_USEC_PER_SEC")]
        public const usec USEC_PER_SEC;

        [CCode (cname="PA_NSEC_PER_SEC")]
        public const uint64 NSEC_PER_SEC;


        [CCode (cname="PA_USEC_PER_MSEC")]
        public const usec USEC_PER_MSEC;

        [CCode (cname="PA_NSEC_PER_MSEC")]
        public const uint64 NSEC_PER_MSEC;


        [CCode (cname="PA_NSEC_PER_USEC")]
        public const uint64 NSEC_PER_USEC;
}
