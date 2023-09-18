using GLib;

namespace PulseAudio {
        [Compact]
        [CCode (cheader_filename="pulse/glib-mainloop.h", cname="pa_glib_mainloop", cprefix="pa_glib_mainloop_", free_function="pa_glib_mainloop_free")]
        public class GLibMainLoop {

                [CCode (cname="pa_glib_mainloop_new")]
                public GLibMainLoop(MainContext? c = null);

                public unowned MainLoopApi get_api();
        }
}
