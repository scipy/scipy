import gdb
import glib_gdb
import sys

if sys.version_info[0] >= 3:
    long = int
else:
    import itertools

    map = itertools.imap

# FrameDecorator is new in gdb 7.7, so we adapt to its absence.
try:
    import gdb.FrameDecorator

    HAVE_GDB_FRAMEDECORATOR = True
    FrameDecorator = gdb.FrameDecorator.FrameDecorator
except ImportError:
    HAVE_GDB_FRAMEDECORATOR = False


# This is not quite right, as local vars may override symname
def read_global_var(symname):
    return gdb.selected_frame().read_var(symname)


def g_type_to_typenode(gtype):
    def lookup_fundamental_type(typenode):
        if typenode == 0:
            return None
        val = read_global_var("static_fundamental_type_nodes")
        if val is None:
            return None
        return val[typenode >> 2].address

    gtype = long(gtype)
    typenode = gtype - gtype % 4
    if typenode > (255 << 2):
        typenode = gdb.Value(typenode).cast(gdb.lookup_type("TypeNode").pointer())
    else:
        typenode = lookup_fundamental_type(typenode)
    return typenode


def g_type_to_name(gtype):
    typenode = g_type_to_typenode(gtype)
    if typenode is not None:
        return glib_gdb.g_quark_to_string(typenode["qname"])
    return None


def is_g_type_instance(val):
    def is_g_type_instance_helper(type):
        if str(type) == "GTypeInstance":
            return True

        while type.code == gdb.TYPE_CODE_TYPEDEF:
            type = type.target()

        if type.code != gdb.TYPE_CODE_STRUCT:
            return False

        fields = type.fields()
        if len(fields) < 1:
            return False

        first_field = fields[0]
        return is_g_type_instance_helper(first_field.type)

    type = val.type
    if type.code != gdb.TYPE_CODE_PTR:
        return False
    type = type.target()
    return is_g_type_instance_helper(type)


def g_type_name_from_instance(instance):
    if long(instance) != 0:
        try:
            inst = instance.cast(gdb.lookup_type("GTypeInstance").pointer())
            klass = inst["g_class"]
            gtype = klass["g_type"]
            name = g_type_to_name(gtype)
            return name
        except RuntimeError:
            pass
    return None


class GTypePrettyPrinter:
    "Prints a GType instance pointer"

    def __init__(self, val):
        self.val = val

    def to_string(self):
        name = g_type_name_from_instance(self.val)
        if name:
            return ("0x%x [%s]") % (long(self.val), name)
        return ("0x%x") % (long(self.val))


def is_g_type_class_instance(val):
    type = val.type
    if type.code != gdb.TYPE_CODE_PTR:
        return False
    return str(type.target()) == "GTypeClass"


class GTypeHandlePrettyPrinter:
    "Prints a GType instance"

    def __init__(self, val, hint=""):
        self.val = val
        self.hint = hint

    def to_string(self):
        typenode = g_type_to_typenode(self.val)
        if typenode is not None:
            name = glib_gdb.g_quark_to_string(typenode["qname"])
            s = ("0x%x [%s%s") % (long(self.val), self.hint, name)
            for i in range(1, int(typenode["n_supers"])):
                node = g_type_to_typenode(typenode["supers"][i])
                if node:
                    name = glib_gdb.g_quark_to_string(node["qname"])
                else:
                    name = "???"
                s += "/" + name
            return s + "]"
        else:
            return ("0x%x") % (long(self.val))


def pretty_printer_lookup(val):
    if is_g_type_instance(val):
        return GTypePrettyPrinter(val)
    if str(val.type) == "GType":
        return GTypeHandlePrettyPrinter(val)
    if is_g_type_class_instance(val):
        return GTypeHandlePrettyPrinter(val["g_type"], "g_type: ")

    return None


def get_signal_name(id):
    if id is None:
        return None
    id = long(id)
    if id == 0:
        return None
    val = read_global_var("g_signal_nodes")
    max_s = read_global_var("g_n_signal_nodes")
    max_s = long(max_s)
    if id < max_s:
        return val[id]["name"].string()
    return None


def frame_name(frame):
    return str(frame.function())


def frame_var(frame, var):
    return frame.inferior_frame().read_var(var)


class SignalFrame(FrameDecorator):
    def __init__(self, frames):
        FrameDecorator.__init__(self, frames[-1])
        self.frame = frames[-1]
        self.frames = frames

    def name(self):
        return "signal-emission"

    def read_var(self, frame, name, array=None):
        try:
            v = frame_var(frame, name)
            if v is None or v.is_optimized_out:
                return None
            if array is not None:
                array.append(v)
            return v
        except ValueError:
            return None

    def read_object(self, frame, name, array=None):
        try:
            v = frame_var(frame, name)
            if v is None or v.is_optimized_out:
                return None
            v = v.cast(gdb.lookup_type("GObject").pointer())
            # Ensure this is a somewhat correct object pointer
            if v is not None and g_type_name_from_instance(v):
                if array is not None:
                    array.append(v)
                return v
            return None
        except ValueError:
            return None

    def append(self, array, obj):
        if obj is not None:
            array.append(obj)

    def or_join_array(self, array):
        if len(array) == 0:
            return "???"
        else:
            return " or ".join(set(map(str, array)))

    def get_detailed_signal_from_frame(self, frame, signal):
        detail = self.read_var(frame, "detail")
        detail = glib_gdb.g_quark_to_string(detail)
        if detail is not None:
            return signal + ":" + detail
        else:
            return signal

    def function(self):
        instances = []
        signals = []

        for frame in self.frames:
            name = frame_name(frame)
            if name == "signal_emit_unlocked_R":
                self.read_object(frame, "instance", instances)
                node = self.read_var(frame, "node")
                if node:
                    signal = node["name"].string()
                    signal = self.get_detailed_signal_from_frame(frame, signal)
                    self.append(signals, signal)

            if name == "g_signal_emitv":
                instance_and_params = self.read_var(frame, "instance_and_params")
                if instance_and_params:
                    instance = instance_and_params[0]["v_pointer"].cast(
                        gdb.Type("GObject").pointer()
                    )
                    self.append(instances, instance)
                id = self.read_var(frame, "signal_id")
                signal = get_signal_name(id)
                if signal:
                    signal = self.get_detailed_signal_from_frame(frame, signal)
                    self.append(signals, signal)

            if name == "g_signal_emit_valist" or name == "g_signal_emit":
                self.read_object(frame, "instance", instances)
                id = self.read_var(frame, "signal_id")
                signal = get_signal_name(id)
                if signal:
                    signal = self.get_detailed_signal_from_frame(frame, signal)
                    self.append(signals, signal)

            if name == "g_signal_emit_by_name":
                self.read_object(frame, "instance", instances)
                self.read_var(frame, "detailed_signal", signals)
                break

        instance = self.or_join_array(instances)
        signal = self.or_join_array(signals)

        return "<emit signal '%s' on instance %s>" % (signal, instance)

    def elided(self):
        return self.frames[0:-1]

    def describe(self, stream, full):
        stream.write(" " + self.function() + "\n")


class GFrameDecorator:
    def __init__(self, iter):
        self.queue = []
        self.iter = iter

    def __iter__(self):
        return self

    def fill(self):
        while len(self.queue) <= 8:
            try:
                f = next(self.iter)
                self.queue.append(f)
            except StopIteration:
                return

    def find_signal_emission(self):
        for i in range(min(len(self.queue), 3)):
            name = frame_name(self.queue[i])
            if name == "signal_emit_unlocked_R" or name == "_g_closure_invoke_va":
                return i
        return -1

    def next(self):
        # Ensure we have enough frames for a full signal emission
        self.fill()

        # Are we at the end?
        if len(self.queue) == 0:
            raise StopIteration

        emission = self.find_signal_emission()
        if emission > 0:
            start = emission
            while True:
                if start == 0:
                    break
                prev_name = frame_name(self.queue[start - 1])
                if prev_name.find("_marshal") >= 0 or prev_name == "g_closure_invoke":
                    start = start - 1
                else:
                    break
            end = emission + 1
            while end < len(self.queue):
                if frame_name(self.queue[end]) in [
                    "g_signal_emitv",
                    "g_signal_emit_valist",
                    "g_signal_emit",
                    "g_signal_emit_by_name",
                    "signal_emitv_unlocked",
                    "signal_emit_valist_unlocked",
                ]:
                    end = end + 1
                else:
                    break

            signal_frames = self.queue[start:end]
            new_frames = [SignalFrame(signal_frames)]
            self.queue[start:end] = new_frames

        return self.queue.pop(0)

    def __next__(self):
        return self.next()


class GFrameFilter(object):
    name = "glib"
    enabled = True
    priority = 100

    def filter(self, iterator):
        return GFrameDecorator(iterator)


def register(obj):
    if obj is None:
        obj = gdb

    if HAVE_GDB_FRAMEDECORATOR:
        filter = GFrameFilter()
        obj.frame_filters[filter.name] = filter
    obj.pretty_printers.append(pretty_printer_lookup)
