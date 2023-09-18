set(PULSEAUDIO_FOUND TRUE)

set(PULSEAUDIO_VERSION_MAJOR 16)
set(PULSEAUDIO_VERSION_MINOR 1)
set(PULSEAUDIO_VERSION 16.1)
set(PULSEAUDIO_VERSION_STRING "16.1")

find_path(PULSEAUDIO_INCLUDE_DIR pulse/pulseaudio.h HINTS "/workspaces/scipy/venv/include")
find_library(PULSEAUDIO_LIBRARY NAMES pulse libpulse HINTS "/workspaces/scipy/venv/lib")
find_library(PULSEAUDIO_MAINLOOP_LIBRARY NAMES pulse-mainloop-glib libpulse-mainloop-glib HINTS "/workspaces/scipy/venv/lib")
