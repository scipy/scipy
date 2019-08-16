
#include <Python.h>
#include <algorithm>
#include <utility>
#include <new>
#include <cstddef>
#include <vector>
#include <unordered_map>
#include <string>

namespace {

/** Handle to a python object that automatically DECREFs */
class py_ref
{
  explicit py_ref(PyObject * object): obj_(object) {}
public:

  py_ref() noexcept: obj_(nullptr) {}
  py_ref(std::nullptr_t) noexcept: py_ref() {}

  py_ref(const py_ref & other) noexcept: obj_(other.obj_) { Py_XINCREF(obj_); }
  py_ref(py_ref && other) noexcept: obj_(other.obj_) { other.obj_ = nullptr; }

  /** Construct from new reference (No INCREF) */
  static py_ref steal(PyObject * object) { return py_ref(object); }

  /** Construct from borrowed reference (and INCREF) */
  static py_ref ref(PyObject * object)
    {
      Py_XINCREF(object);
      return py_ref(object);
    }

  ~py_ref(){ Py_XDECREF(obj_); }

  py_ref & operator = (const py_ref & other) noexcept
    {
      py_ref(other).swap(*this);
      return *this;
    }

  py_ref & operator = (py_ref && other) noexcept
    {
      py_ref(std::move(other)).swap(*this);
      return *this;
    }

  void swap(py_ref & other) noexcept
    {
      std::swap(other.obj_, obj_);
    }

  explicit operator bool () const { return obj_ != nullptr; }

  operator PyObject* () const { return get(); }

  PyObject * get() const { return obj_; }
  PyObject * release()
    {
      PyObject * t = obj_;
      obj_ = nullptr;
      return t;
    }
  void reset()
    {
      Py_CLEAR(obj_);
    }
private:
  PyObject * obj_;
};

/** Make tuple from variadic set of PyObjects */
template <typename ... Ts>
py_ref py_make_tuple(Ts... args)
{
  using py_obj = PyObject *;
  return py_ref::steal(PyTuple_Pack(sizeof...(args), py_obj{args}...));
}

struct backend_options
{
  py_ref backend;
  bool coerce, only;

  bool operator == (const backend_options & other) const
    {
      return (backend == other.backend
              && coerce == other.coerce
              && only == other.only);
    }
};

struct global_backends
{
  backend_options global;
  std::vector<py_ref> registered;
};

struct local_backends
{
  std::vector<py_ref> skipped;
  std::vector<backend_options> preferred;
};


static py_ref BackendNotImplementedError;
static std::unordered_map<std::string, global_backends> global_domain_map;
thread_local std::unordered_map<
  std::string, local_backends> local_domain_map;

/** Constant Python string identifiers

Using these with PyObject_GetAttr is faster than PyObject_GetAttrString which
has to create a new python string internally.
 */
struct
{
  py_ref ua_convert;
  py_ref ua_domain;
  py_ref ua_function;

  bool init()
    {
      ua_convert = py_ref::steal(PyUnicode_InternFromString("__ua_convert__"));
      if (!ua_convert)
        return false;

      ua_domain = py_ref::steal(PyUnicode_InternFromString("__ua_domain__"));
      if (!ua_domain)
        return false;

      ua_function = py_ref::steal(PyUnicode_InternFromString("__ua_function__"));
      if (!ua_function)
        return false;

      return true;
    }

  void clear()
    {
      ua_convert.reset();
      ua_domain.reset();
      ua_function.reset();
    }
} identifiers;


std::string domain_to_string(PyObject * domain)
{
  if (!PyUnicode_Check(domain))
  {
    PyErr_SetString(PyExc_TypeError, "__ua_domain__ must be a string");
    return {};
  }

  Py_ssize_t size;
  const char * str = PyUnicode_AsUTF8AndSize(domain, &size);
  if (!str)
    return {};

  if (size == 0)
  {
    PyErr_SetString(PyExc_ValueError, "__ua_domain__ must be non-empty");
    return {};
  }

  return std::string(str, size);
}

std::string backend_to_domain_string(PyObject * backend)
{
  auto domain = py_ref::steal(
    PyObject_GetAttr(backend, identifiers.ua_domain));
  if (!domain)
    return {};

  return domain_to_string(domain);
}


/** Use to clean up python references before the interpreter is finalized.
 *
 * This must be installed in a python atexit handler. This prevents Py_DECREF
 * being called after the interpreter has already shudown.
 */
PyObject * clear_all_globals(PyObject * /* self */, PyObject * /* args */)
{
  global_domain_map.clear();
  BackendNotImplementedError.reset();
  identifiers.clear();
  Py_RETURN_NONE;
}


PyObject * set_global_backend(PyObject * /* self */, PyObject * args)
{
  PyObject * backend;
  int only = false, coerce = false;
  if (!PyArg_ParseTuple(args, "O|pp",
    &backend,
    &coerce,
    &only))
    return nullptr;

  auto domain = backend_to_domain_string(backend);
  if (domain.empty())
    return nullptr;

  backend_options options;
  options.backend = py_ref::ref(backend);
  options.coerce = coerce;
  options.only = only;

  global_domain_map[domain].global = options;
  Py_RETURN_NONE;
}

PyObject * register_backend(PyObject * /* self */, PyObject * args)
{
  PyObject * backend;
  if (!PyArg_ParseTuple(args, "O", 
    &backend))
    return nullptr;

  auto domain = backend_to_domain_string(backend);
  if (domain.empty())
    return nullptr;

  global_domain_map[domain].registered.push_back(py_ref::ref(backend));
  Py_RETURN_NONE;
}

void clear_single(const std::string& domain, bool registered, bool global)
{
  auto domain_globals = global_domain_map.find(domain);
  if (domain_globals == global_domain_map.end())
    return;

  if (registered && global)
    {
      global_domain_map.erase(domain_globals);
      return;
    }
  
  if (registered)
  {
    domain_globals->second.registered.clear();
  }

  if (global)
  {
    domain_globals->second.global.backend.reset();
  }
}

PyObject * clear_backends(PyObject * /* self */, PyObject * args)
{
  PyObject* domain = nullptr;
  int registered = true, global = false;
  if (!PyArg_ParseTuple(args, "O|pp", 
    &domain, &registered, &global))
    return nullptr;

  if (domain == Py_None && registered && global)
  {
    global_domain_map.clear();
    Py_RETURN_NONE;
  }

  auto domain_str = domain_to_string(domain);
  clear_single(domain_str, registered, global);
  Py_RETURN_NONE;
}


/** Common functionality of set_backend and skip_backend */
template <typename T>
class context_helper
{
  T new_backend_;
  std::vector<T> * backends_;
  size_t enter_size_;
public:

  context_helper():
    backends_(nullptr),
    enter_size_(size_t(-1))
    {}

  bool init(std::vector<T> & backends, T new_backend)
    {
      backends_ = &backends;
      new_backend_ = std::move(new_backend);
      return true;
    }

  bool enter()
    {
      enter_size_ = backends_->size();
      try { backends_->push_back(new_backend_); }
      catch(std::bad_alloc&)
      {
        PyErr_NoMemory();
        return false;
      }
      return true;
    }

  bool exit()
    {
      bool success = (enter_size_ + 1 == backends_->size()
                      && backends_->back() == new_backend_);
      if (enter_size_ < backends_->size())
        backends_->resize(enter_size_);

      if (!success)
        PyErr_SetString(PyExc_RuntimeError,
                        "Found invalid context state while in __exit__");
      return success;
    }
};


struct SetBackendContext
{
  PyObject_HEAD

  context_helper<backend_options> ctx_;

  static void dealloc(SetBackendContext * self)
    {
      self->~SetBackendContext();
      Py_TYPE(self)->tp_free(self);
    }

  static PyObject * new_(PyTypeObject * type, PyObject * args, PyObject * kwargs)
    {
      auto self = reinterpret_cast<SetBackendContext *>(type->tp_alloc(type, 0));
      if (self == nullptr)
        return nullptr;

      // Placement new
      self = new (self) SetBackendContext;
      return reinterpret_cast<PyObject *>(self);
    }

  static int init(
    SetBackendContext * self, PyObject * args, PyObject * kwargs)
    {
      static const char * kwlist[] = {"backend", "coerce", "only", nullptr};
      PyObject * backend = nullptr;
      int coerce = false;
      int only = false;

      if (!PyArg_ParseTupleAndKeywords(
            args, kwargs,
            "O|pp", (char**)kwlist,
            &backend,
            &coerce,
            &only))
        return -1;


      auto domain = backend_to_domain_string(backend);
      if (domain.empty())
        return -1;
      backend_options opt;
      opt.backend = py_ref::ref(backend);
      opt.coerce = coerce;
      opt.only = only;

      try
      {
        if (!self->ctx_.init(local_domain_map[domain].preferred, std::move(opt)))
          return -1;
      }
      catch(std::bad_alloc&)
      {
        PyErr_NoMemory();
        return -1;
      }

      return 0;
    }

  static PyObject * enter__(SetBackendContext * self, PyObject * /* args */)
    {
      if (!self->ctx_.enter())
        return nullptr;
      Py_RETURN_NONE;
    }

  static PyObject * exit__(SetBackendContext * self, PyObject * /*args*/)
    {
      if (!self->ctx_.exit())
        return nullptr;
      Py_RETURN_NONE;
    }
};


struct SkipBackendContext
{
  PyObject_HEAD

  context_helper<py_ref> ctx_;

  static void dealloc(SkipBackendContext * self)
    {
      self->~SkipBackendContext();
      Py_TYPE(self)->tp_free(self);
    }

  static PyObject * new_(PyTypeObject * type, PyObject * args, PyObject * kwargs)
    {
      auto self = reinterpret_cast<SkipBackendContext *>(type->tp_alloc(type, 0));
      if (self == nullptr)
        return nullptr;

      // Placement new
      self = new (self) SkipBackendContext;
      return reinterpret_cast<PyObject *>(self);
    }

  static int init(
    SkipBackendContext * self, PyObject * args, PyObject * kwargs)
    {
      static const char *kwlist[] = {"backend", nullptr};
      PyObject * backend;

      if (!PyArg_ParseTupleAndKeywords(args, kwargs,
                                       "O", (char**)kwlist,
                                       &backend))
        return -1;

      auto domain = backend_to_domain_string(backend);
      if (domain.empty())
        return -1;

      try
      {
        if (!self->ctx_.init(
              local_domain_map[domain].skipped, py_ref::ref(backend)))
          return -1;
      }
      catch(std::bad_alloc&)
      {
        PyErr_NoMemory();
        return -1;
      }

      return 0;
    }

  static PyObject * enter__(SkipBackendContext * self, PyObject * /* args */)
    {
      if (!self->ctx_.enter())
        return nullptr;
      Py_RETURN_NONE;
    }

  static PyObject * exit__(SkipBackendContext * self, PyObject * /*args*/)
    {
      if (!self->ctx_.exit())
        return nullptr;
      Py_RETURN_NONE;
    }
};

enum class LoopReturn { Continue, Break, Error };

template <typename Callback>
LoopReturn for_each_backend(const std::string & domain_key, Callback call)
{
  local_backends * locals = nullptr;
  try
  {
    locals = &local_domain_map[domain_key];
  }
  catch (std::bad_alloc&)
  {
    PyErr_NoMemory();
    return LoopReturn::Error;
  }


  auto & skip = locals->skipped;
  auto & pref = locals->preferred;

  auto should_skip =
    [&](PyObject * backend)
      {
        auto it = std::find_if(
          skip.begin(), skip.end(),
          [&](const py_ref & be) { return be.get() == backend; });

        return (it != skip.end());
      };

  LoopReturn ret = LoopReturn::Continue;
  for (int i = pref.size()-1; i >= 0; --i)
  {
    auto options = pref[i];
    if (should_skip(options.backend))
      continue;

    ret = call(options.backend.get(), options.coerce);
    if (ret != LoopReturn::Continue)
      return ret;

    if (options.only || options.coerce)
      return ret;
  }

  auto& globals = global_domain_map[domain_key];
  auto& global_options = globals.global;
  if (global_options.backend && !should_skip(global_options.backend))
  {
    ret = call(global_options.backend.get(), global_options.coerce);
    if (ret != LoopReturn::Continue)
      return ret;

    if (global_options.only || global_options.coerce)
      return ret;
  }

  for (size_t i = 0; i < globals.registered.size(); ++i)
  {
    py_ref backend = globals.registered[i];
    if (should_skip(backend))
      continue;

    ret = call(backend.get(), false);
    if (ret != LoopReturn::Continue)
      return ret;
  }
  return ret;
}

struct py_func_args { py_ref args, kwargs; };

struct Function
{
  PyObject_HEAD
  py_ref extractor_, replacer_;   // functions to handle dispatchables
  std::string domain_key_;        // associated __ua_domain__ in UTF8
  py_ref def_args_, def_kwargs_;  // default arguments
  py_ref def_impl_;               // default implementation
  py_ref dict_;                   // __dict__

  PyObject * call(PyObject * args, PyObject * kwargs);

  py_func_args replace_dispatchables(
    PyObject * backend, PyObject * args, PyObject * kwargs, PyObject * coerce);

  py_ref canonicalize_args(PyObject * args);
  py_ref canonicalize_kwargs(PyObject * kwargs);

  static void dealloc(Function * self)
    {
      PyObject_GC_UnTrack(self);
      self->~Function();
      Py_TYPE(self)->tp_free(self);
    }

  static PyObject * new_(PyTypeObject * type, PyObject * args, PyObject * kwargs)
    {
      auto self = reinterpret_cast<Function *>(type->tp_alloc(type, 0));
      if (self == nullptr)
        return nullptr;

      // Placement new
      self = new (self) Function;
      return reinterpret_cast<PyObject *>(self);
    }

  static int init(Function * self, PyObject * args, PyObject * /*kwargs*/)
    {
      PyObject * extractor, * replacer;
      PyObject * domain;
      PyObject * def_args, * def_kwargs;
      PyObject * def_impl;

      if (!PyArg_ParseTuple(
            args, "OOO!O!O!O",
            &extractor,
            &replacer,
            &PyUnicode_Type, &domain,
            &PyTuple_Type, &def_args,
            &PyDict_Type, &def_kwargs,
            &def_impl))
      {
        return -1;
      }

      if (!PyCallable_Check(extractor)
          || (replacer != Py_None && !PyCallable_Check(replacer)))
      {
        PyErr_SetString(PyExc_TypeError,
                        "Argument extractor and replacer must be callable");
        return -1;
      }

      if (def_impl != Py_None && !PyCallable_Check(def_impl))
      {
        PyErr_SetString(PyExc_TypeError,
                        "Default implementation must be Callable or None");
        return -1;
      }

      self->domain_key_ = domain_to_string(domain);
      if (PyErr_Occurred())
        return -1;

      self->extractor_ = py_ref::ref(extractor);
      self->replacer_ = py_ref::ref(replacer);
      self->def_args_ = py_ref::ref(def_args);
      self->def_kwargs_ = py_ref::ref(def_kwargs);
      self->def_impl_ = py_ref::ref(def_impl);

      return 0;
    }

  static PyObject * repr(Function * self);
  static PyObject * descr_get(PyObject * self, PyObject * obj, PyObject * type);
  static int traverse(Function * self, visitproc visit, void * arg);
  static int clear(Function * self);
};


bool is_default(PyObject * value, PyObject * def)
{
  // TODO: richer comparison for builtin types? (if cheap)
  return (value == def);
}


py_ref Function::canonicalize_args(PyObject * args)
{
  const auto arg_size = PyTuple_GET_SIZE(args);
  const auto def_size = PyTuple_GET_SIZE(def_args_.get());

  if (arg_size > def_size)
    return py_ref::ref(args);

  Py_ssize_t mismatch = 0;
  for (Py_ssize_t i = arg_size - 1; i >= 0; --i)
  {
    auto val = PyTuple_GET_ITEM(args, i);
    auto def = PyTuple_GET_ITEM(def_args_.get(), i);
    if (!is_default(val, def))
    {
      mismatch = i + 1;
      break;
    }
  }

  return py_ref::steal(PyTuple_GetSlice(args, 0, mismatch));
}


py_ref Function::canonicalize_kwargs(PyObject * kwargs)
{
  if (kwargs == nullptr)
    return py_ref::steal(PyDict_New());

  PyObject * key, * def_value;
  Py_ssize_t pos = 0;
  while (PyDict_Next(def_kwargs_, &pos, &key, &def_value))
  {
    auto val = PyDict_GetItem(kwargs, key);
    if (val && is_default(val, def_value))
    {
      PyDict_DelItem(kwargs, key);
    }
  }
  return py_ref::ref(kwargs);
}


py_func_args Function::replace_dispatchables(
  PyObject * backend, PyObject * args, PyObject * kwargs, PyObject * coerce)
{
  auto ua_convert = py_ref::steal(
    PyObject_GetAttr(backend, identifiers.ua_convert));

  if (!ua_convert)
  {
    PyErr_Clear();
    return {py_ref::ref(args), py_ref::ref(kwargs)};
  }

  auto dispatchables = py_ref::steal(PyObject_Call(extractor_, args, kwargs));
  if (!dispatchables)
    return {};

  auto convert_args = py_make_tuple(dispatchables, coerce);
  auto res = py_ref::steal(PyObject_Call(ua_convert, convert_args, nullptr));
  if (!res)
  {
    return {};
  }

  if (res == Py_NotImplemented)
  {
    return {std::move(res), nullptr};
  }

  auto replaced_args = py_ref::steal(PySequence_Tuple(res));
  if (!replaced_args)
    return {};

  auto replacer_args = py_make_tuple(args, kwargs, replaced_args);
  if (!replacer_args)
    return {};

  res = py_ref::steal(PyObject_Call(replacer_, replacer_args, nullptr));
  if (!res)
    return {};

  if (!PyTuple_Check(res) || PyTuple_Size(res) != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "Argument replacer must return a 2-tuple (args, kwargs)");
    return {};
  }

  auto new_args = py_ref::ref(PyTuple_GET_ITEM(res.get(), 0));
  auto new_kwargs = py_ref::ref(PyTuple_GET_ITEM(res.get(), 1));

  new_kwargs = canonicalize_kwargs(new_kwargs);

  if (!PyTuple_Check(new_args) || !PyDict_Check(new_kwargs))
  {
    PyErr_SetString(PyExc_ValueError, "Invalid return from argument_replacer");
    return {};
  }

  return {std::move(new_args), std::move(new_kwargs)};
}


PyObject * Function_call(
  Function * self, PyObject * args, PyObject * kwargs)
{
  return self->call(args, kwargs);
}


PyObject * Function::call(PyObject * args_, PyObject * kwargs_)
{
  auto args = canonicalize_args(args_);
  auto kwargs = canonicalize_kwargs(kwargs_);

  py_ref result;

  auto ret = for_each_backend(
    domain_key_,
    [&, this](PyObject * backend, bool coerce)
      {
        auto new_args = replace_dispatchables(backend, args, kwargs,
                                              coerce ? Py_True : Py_False);
        if (new_args.args == Py_NotImplemented)
          return LoopReturn::Continue;
        if (new_args.args == nullptr)
          return LoopReturn::Error;

        auto ua_function = py_ref::steal(
          PyObject_GetAttr(backend, identifiers.ua_function));
        if (!ua_function)
          return LoopReturn::Error;

        auto ua_func_args = py_make_tuple(
          reinterpret_cast<PyObject *>(this), new_args.args, new_args.kwargs);
        if (!ua_func_args)
          return LoopReturn::Error;

        result = py_ref::steal(
          PyObject_Call(ua_function, ua_func_args, nullptr));

        // Try the default with this backend
        if (result == Py_NotImplemented && def_impl_ != Py_None)
        {
          backend_options opt;
          opt.backend = py_ref::ref(backend);
          opt.coerce = coerce;
          opt.only = true;
          context_helper<backend_options> ctx;
          try
          {
            if (!ctx.init(
                  local_domain_map[domain_key_].preferred, std::move(opt)))
              return LoopReturn::Error;
          }
          catch(std::bad_alloc&)
          {
            PyErr_NoMemory();
            return LoopReturn::Error;
          }

          if (!ctx.enter())
            return LoopReturn::Error;

          result = py_ref::steal(
            PyObject_Call(def_impl_, new_args.args, new_args.kwargs));

          if (PyErr_Occurred() && PyErr_ExceptionMatches(BackendNotImplementedError))
          {
            PyErr_Clear();  // Suppress exception
            result = py_ref::ref(Py_NotImplemented);
          }

          if (!ctx.exit())
            return LoopReturn::Error;
        }

        if (!result)
          return LoopReturn::Error;

        if (result == Py_NotImplemented)
          return LoopReturn::Continue;

        return LoopReturn::Break;  // Backend called successfully
      }
    );

  if (ret != LoopReturn::Continue)
    return result.release();

  if (def_impl_ == Py_None)
  {
    PyErr_SetString(
      BackendNotImplementedError,
      "No selected backends had an implementation for this function.");
    return nullptr;
  }

  return PyObject_Call(def_impl_, args, kwargs);
}


PyObject * Function::repr(Function * self)
{
  if (self->dict_)
    if (auto name = PyDict_GetItemString(self->dict_, "__name__"))
      return PyUnicode_FromFormat("<uarray multimethod '%S'>", name);

  return PyUnicode_FromString("<uarray multimethod>");
}


/** Implements the descriptor protocol to allow binding to class instances */
PyObject * Function::descr_get(PyObject * self, PyObject * obj, PyObject * type)
{
  if (!obj)
  {
    Py_INCREF(self);
    return self;
  }

  return PyMethod_New(self, obj);
}


/** Make members visible to the garbage collector */
int Function::traverse(Function * self, visitproc visit, void * arg)
{
  Py_VISIT(self->extractor_);
  Py_VISIT(self->replacer_);
  Py_VISIT(self->def_args_);
  Py_VISIT(self->def_kwargs_);
  Py_VISIT(self->def_impl_);
  Py_VISIT(self->dict_);
  return 0;
}


/** Break reference cycles when being GCed */
int Function::clear(Function * self)
{
  self->extractor_.reset();
  self->replacer_.reset();
  self->def_args_.reset();
  self->def_kwargs_.reset();
  self->def_impl_.reset();
  self->dict_.reset();
  return 0;
}


PyGetSetDef Function_getset[] =
{
  {"__dict__", PyObject_GenericGetDict, PyObject_GenericSetDict},
  {NULL} /* Sentinel */
};

PyTypeObject FunctionType = {
  PyVarObject_HEAD_INIT(NULL, 0)
  "uarray._Function",              /* tp_name */
  sizeof(Function),                /* tp_basicsize */
  0,                               /* tp_itemsize */
  (destructor)Function::dealloc,   /* tp_dealloc */
  0,                               /* tp_print */
  0,                               /* tp_getattr */
  0,                               /* tp_setattr */
  0,                               /* tp_reserved */
  (reprfunc)Function::repr,        /* tp_repr */
  0,                               /* tp_as_number */
  0,                               /* tp_as_sequence */
  0,                               /* tp_as_mapping */
  0,                               /* tp_hash  */
  (ternaryfunc)Function_call,      /* tp_call */
  0,                               /* tp_str */
  PyObject_GenericGetAttr,         /* tp_getattro */
  PyObject_GenericSetAttr,         /* tp_setattro */
  0,                               /* tp_as_buffer */
  (Py_TPFLAGS_DEFAULT
   | Py_TPFLAGS_HAVE_GC),          /* tp_flags */
  0,                               /* tp_doc */
  (traverseproc)Function::traverse,/* tp_traverse */
  (inquiry)Function::clear,        /* tp_clear */
  0,                               /* tp_richcompare */
  0,                               /* tp_weaklistoffset */
  0,                               /* tp_iter */
  0,                               /* tp_iternext */
  0,                               /* tp_methods */
  0,                               /* tp_members */
  Function_getset,                 /* tp_getset */
  0,                               /* tp_base */
  0,                               /* tp_dict */
  Function::descr_get,             /* tp_descr_get */
  0,                               /* tp_descr_set */
  offsetof(Function, dict_),       /* tp_dictoffset */
  (initproc)Function::init,        /* tp_init */
  0,                               /* tp_alloc */
  Function::new_,                  /* tp_new */
};


PyMethodDef SetBackendContext_Methods[] = {
  {"__enter__", (PyCFunction)SetBackendContext::enter__, METH_NOARGS, nullptr},
  {"__exit__", (PyCFunction)SetBackendContext::exit__, METH_VARARGS, nullptr},
  {NULL}  /* Sentinel */
};

PyTypeObject SetBackendContextType = {
  PyVarObject_HEAD_INIT(NULL, 0)
  "uarray._SetBackendContext",             /* tp_name */
  sizeof(SetBackendContext),               /* tp_basicsize */
  0,                                       /* tp_itemsize */
  (destructor)SetBackendContext::dealloc,  /* tp_dealloc */
  0,                                       /* tp_print */
  0,                                       /* tp_getattr */
  0,                                       /* tp_setattr */
  0,                                       /* tp_reserved */
  0,                                       /* tp_repr */
  0,                                       /* tp_as_number */
  0,                                       /* tp_as_sequence */
  0,                                       /* tp_as_mapping */
  0,                                       /* tp_hash  */
  0,                                       /* tp_call */
  0,                                       /* tp_str */
  0,                                       /* tp_getattro */
  0,                                       /* tp_setattro */
  0,                                       /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT,                      /* tp_flags */
  0,                                       /* tp_doc */
  0,                                       /* tp_traverse */
  0,                                       /* tp_clear */
  0,                                       /* tp_richcompare */
  0,                                       /* tp_weaklistoffset */
  0,                                       /* tp_iter */
  0,                                       /* tp_iternext */
  SetBackendContext_Methods,               /* tp_methods */
  0,                                       /* tp_members */
  0,                                       /* tp_getset */
  0,                                       /* tp_base */
  0,                                       /* tp_dict */
  0,                                       /* tp_descr_get */
  0,                                       /* tp_descr_set */
  0,                                       /* tp_dictoffset */
  (initproc)SetBackendContext::init,       /* tp_init */
  0,                                       /* tp_alloc */
  SetBackendContext::new_,                 /* tp_new */
};


PyMethodDef SkipBackendContext_Methods[] = {
  {"__enter__", (PyCFunction)SkipBackendContext::enter__, METH_NOARGS, nullptr},
  {"__exit__", (PyCFunction)SkipBackendContext::exit__, METH_VARARGS, nullptr},
  {NULL}  /* Sentinel */
};

PyTypeObject SkipBackendContextType = {
  PyVarObject_HEAD_INIT(NULL, 0)
  "uarray._SkipBackendContext",             /* tp_name */
  sizeof(SkipBackendContext),               /* tp_basicsize */
  0,                                        /* tp_itemsize */
  (destructor)SkipBackendContext::dealloc,  /* tp_dealloc */
  0,                                        /* tp_print */
  0,                                        /* tp_getattr */
  0,                                        /* tp_setattr */
  0,                                        /* tp_reserved */
  0,                                        /* tp_repr */
  0,                                        /* tp_as_number */
  0,                                        /* tp_as_sequence */
  0,                                        /* tp_as_mapping */
  0,                                        /* tp_hash  */
  0,                                        /* tp_call */
  0,                                        /* tp_str */
  0,                                        /* tp_getattro */
  0,                                        /* tp_setattro */
  0,                                        /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT,                       /* tp_flags */
  0,                                        /* tp_doc */
  0,                                        /* tp_traverse */
  0,                                        /* tp_clear */
  0,                                        /* tp_richcompare */
  0,                                        /* tp_weaklistoffset */
  0,                                        /* tp_iter */
  0,                                        /* tp_iternext */
  SkipBackendContext_Methods,               /* tp_methods */
  0,                                        /* tp_members */
  0,                                        /* tp_getset */
  0,                                        /* tp_base */
  0,                                        /* tp_dict */
  0,                                        /* tp_descr_get */
  0,                                        /* tp_descr_set */
  0,                                        /* tp_dictoffset */
  (initproc)SkipBackendContext::init,       /* tp_init */
  0,                                        /* tp_alloc */
  SkipBackendContext::new_,                 /* tp_new */
};


PyMethodDef method_defs[] =
{
  {"set_global_backend", set_global_backend, METH_VARARGS, nullptr},
  {"register_backend", register_backend, METH_VARARGS, nullptr},
  {"clear_all_globals", clear_all_globals, METH_NOARGS, nullptr},
  {"clear_backends", clear_backends, METH_VARARGS, nullptr},
  {NULL} /* Sentinel */
};

PyModuleDef uarray_module =
{
  PyModuleDef_HEAD_INIT,
  "_uarray",
  nullptr,
  -1,
  method_defs,
};

}  // namespace (anonymous)


#if defined(WIN32) || defined(_WIN32)
#  define MODULE_EXPORT __declspec(dllexport)
#else
#  define MODULE_EXPORT __attribute__ ((visibility("default")))
#endif

extern "C" MODULE_EXPORT PyObject *
PyInit__uarray(void)
{

  auto m = py_ref::steal(PyModule_Create(&uarray_module));
  if (!m)
    return nullptr;

  if (PyType_Ready(&FunctionType) < 0)
    return nullptr;
  Py_INCREF(&FunctionType);
  PyModule_AddObject(m, "_Function", (PyObject *)&FunctionType);

  if (PyType_Ready(&SetBackendContextType) < 0)
    return nullptr;
  Py_INCREF(&SetBackendContextType);
  PyModule_AddObject(m, "_SetBackendContext", (PyObject*)&SetBackendContextType);

  if (PyType_Ready(&SkipBackendContextType) < 0)
    return nullptr;
  Py_INCREF(&SkipBackendContextType);
  PyModule_AddObject(
    m, "_SkipBackendContext", (PyObject*)&SkipBackendContextType);

  BackendNotImplementedError = py_ref::steal(
    PyErr_NewExceptionWithDoc(
      "uarray.BackendNotImplementedError",
      "An exception that is thrown when no compatible"
      " backend is found for a method.",
      PyExc_NotImplementedError,
      nullptr));
  if (!BackendNotImplementedError)
    return nullptr;
  Py_INCREF(BackendNotImplementedError.get());
  PyModule_AddObject(
    m, "BackendNotImplementedError", BackendNotImplementedError);

  if (!identifiers.init())
    return nullptr;

  return m.release();
}
