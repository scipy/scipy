#include <Python.h>

#include "small_dynamic_array.h"
#include "vectorcall.h"

#include <algorithm>
#include <cstddef>
#include <new>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>


namespace {

template <typename T>
class immortal {
    alignas(T) std::byte storage[sizeof(T)];

public:
    template <typename... Args>
    immortal(Args&&... args) {
       // Construct new T in storage
       new(&storage) T(std::forward<Args>(args)...);
    }
    ~immortal() {
        // Intentionally don't call destructor
    }

    T* get() { return reinterpret_cast<T*>(&storage); }
    const T* get() const { return reinterpret_cast<const T*>(&storage); }
    const T* get_const() const { return reinterpret_cast<const T*>(&storage); }

    const T* operator ->() const { return get(); }
    T* operator ->() { return get(); }

    T& operator*() { return *get(); }
    const T& operator*() const { return *get(); }

};

/** Handle to a python object that automatically DECREFs */
class py_ref {
  explicit py_ref(PyObject * object): obj_(object) {}

public:
  py_ref() noexcept: obj_(nullptr) {}
  py_ref(std::nullptr_t) noexcept: py_ref() {}

  py_ref(const py_ref & other) noexcept: obj_(other.obj_) { Py_XINCREF(obj_); }
  py_ref(py_ref && other) noexcept: obj_(other.obj_) { other.obj_ = nullptr; }

  /** Construct from new reference (No INCREF) */
  static py_ref steal(PyObject * object) { return py_ref(object); }

  /** Construct from borrowed reference (and INCREF) */
  static py_ref ref(PyObject * object) {
    Py_XINCREF(object);
    return py_ref(object);
  }

  ~py_ref() { Py_XDECREF(obj_); }

  py_ref & operator=(const py_ref & other) noexcept {
    py_ref(other).swap(*this);
    return *this;
  }

  py_ref & operator=(py_ref && other) noexcept {
    py_ref(std::move(other)).swap(*this);
    return *this;
  }

  friend bool operator==(const py_ref & lhs, const py_ref & rhs) {
    return lhs.obj_ == rhs.obj_;
  }
  friend bool operator==(PyObject * lhs, const py_ref & rhs) {
    return lhs == rhs.obj_;
  }
  friend bool operator==(const py_ref & lhs, PyObject * rhs) {
    return lhs.obj_ == rhs;
  }
  friend bool operator!=(const py_ref & lhs, const py_ref & rhs) {
    return lhs.obj_ != rhs.obj_;
  }
  friend bool operator!=(PyObject * lhs, const py_ref & rhs) {
    return lhs != rhs.obj_;
  }
  friend bool operator!=(const py_ref & lhs, PyObject * rhs) {
    return lhs.obj_ != rhs;
  }

  void swap(py_ref & other) noexcept { std::swap(other.obj_, obj_); }

  explicit operator bool() const { return obj_ != nullptr; }

  PyObject * get() const { return obj_; }
  PyObject * release() {
    PyObject * t = obj_;
    obj_ = nullptr;
    return t;
  }
  void reset() { Py_CLEAR(obj_); }

private:
  PyObject * obj_;
};

PyObject * py_get(const py_ref & ref) { return ref.get(); }
PyObject * py_get(PyObject * obj) { return obj; }

/** Make tuple from variadic set of PyObjects */
template <typename... Ts>
py_ref py_make_tuple(const Ts &... args) {
  return py_ref::steal(PyTuple_Pack(sizeof...(args), py_get(args)...));
}

py_ref py_bool(bool input) { return py_ref::ref(input ? Py_True : Py_False); }

template <typename T, size_t N>
constexpr size_t array_size(const T (&array)[N]) {
  return N;
}

struct backend_options {
  py_ref backend;
  bool coerce = false;
  bool only = false;

  bool operator==(const backend_options & other) const {
    return (
        backend == other.backend && coerce == other.coerce &&
        only == other.only);
  }

  bool operator!=(const backend_options & other) const {
    return !(*this == other);
  }
};

struct global_backends {
  backend_options global;
  std::vector<py_ref> registered;
  bool try_global_backend_last = false;
};

struct local_backends {
  std::vector<py_ref> skipped;
  std::vector<backend_options> preferred;
};

using global_state_t = std::unordered_map<std::string, global_backends>;
using local_state_t = std::unordered_map<std::string, local_backends>;

static py_ref BackendNotImplementedError;
static immortal<global_state_t> global_domain_map;
thread_local global_state_t * current_global_state = global_domain_map.get();
thread_local global_state_t thread_local_domain_map;
thread_local local_state_t local_domain_map;

/** Constant Python string identifiers

Using these with PyObject_GetAttr is faster than PyObject_GetAttrString which
has to create a new python string internally.
 */
struct {
  immortal<py_ref> ua_convert;
  immortal<py_ref> ua_domain;
  immortal<py_ref> ua_function;

  bool init() {
    *ua_convert = py_ref::steal(PyUnicode_InternFromString("__ua_convert__"));
    if (!*ua_convert)
      return false;

    *ua_domain = py_ref::steal(PyUnicode_InternFromString("__ua_domain__"));
    if (!*ua_domain)
      return false;

    *ua_function = py_ref::steal(PyUnicode_InternFromString("__ua_function__"));
    if (!*ua_function)
      return false;

    return true;
  }

  void clear() {
    ua_convert->reset();
    ua_domain->reset();
    ua_function->reset();
  }
} identifiers;

bool domain_validate(PyObject * domain) {
  if (!PyUnicode_Check(domain)) {
    PyErr_SetString(PyExc_TypeError, "__ua_domain__ must be a string");
    return false;
  }

  auto size = PyUnicode_GetLength(domain);
  if (size == 0) {
    PyErr_SetString(PyExc_ValueError, "__ua_domain__ must be non-empty");
    return false;
  }

  return true;
}

std::string domain_to_string(PyObject * domain) {
  if (!domain_validate(domain)) {
    return {};
  }

  Py_ssize_t size;
  const char * str = PyUnicode_AsUTF8AndSize(domain, &size);
  if (!str)
    return {};

  if (size == 0) {
    PyErr_SetString(PyExc_ValueError, "__ua_domain__ must be non-empty");
    return {};
  }

  return std::string(str, size);
}

Py_ssize_t backend_get_num_domains(PyObject * backend) {
  auto domain =
      py_ref::steal(PyObject_GetAttr(backend, identifiers.ua_domain->get()));
  if (!domain)
    return -1;

  if (PyUnicode_Check(domain.get())) {
    return 1;
  }

  if (!PySequence_Check(domain.get())) {
    PyErr_SetString(
        PyExc_TypeError,
        "__ua_domain__ must be a string or sequence of strings");
    return -1;
  }

  return PySequence_Size(domain.get());
}

enum class LoopReturn { Continue, Break, Error };

template <typename Func>
LoopReturn backend_for_each_domain(PyObject * backend, Func f) {
  auto domain =
      py_ref::steal(PyObject_GetAttr(backend, identifiers.ua_domain->get()));
  if (!domain)
    return LoopReturn::Error;

  if (PyUnicode_Check(domain.get())) {
    return f(domain.get());
  }

  if (!PySequence_Check(domain.get())) {
    PyErr_SetString(
        PyExc_TypeError,
        "__ua_domain__ must be a string or sequence of strings");
    return LoopReturn::Error;
  }

  auto size = PySequence_Size(domain.get());
  if (size < 0)
    return LoopReturn::Error;
  if (size == 0) {
    PyErr_SetString(PyExc_ValueError, "__ua_domain__ lists must be non-empty");
    return LoopReturn::Error;
  }

  for (Py_ssize_t i = 0; i < size; ++i) {
    auto dom = py_ref::steal(PySequence_GetItem(domain.get(), i));
    if (!dom)
      return LoopReturn::Error;

    auto res = f(dom.get());
    if (res != LoopReturn::Continue) {
      return res;
    }
  }
  return LoopReturn::Continue;
}

template <typename Func>
LoopReturn backend_for_each_domain_string(PyObject * backend, Func f) {
  return backend_for_each_domain(backend, [&](PyObject * domain) {
    auto domain_string = domain_to_string(domain);
    if (domain_string.empty()) {
      return LoopReturn::Error;
    }
    return f(domain_string);
  });
}

bool backend_validate_ua_domain(PyObject * backend) {
  const auto res = backend_for_each_domain(backend, [&](PyObject * domain) {
    return domain_validate(domain) ? LoopReturn::Continue : LoopReturn::Error;
  });
  return (res != LoopReturn::Error);
}

struct BackendState {
  PyObject_HEAD
  global_state_t globals;
  local_state_t locals;
  bool use_thread_local_globals = true;

  static void dealloc(BackendState * self) {
    self->~BackendState();
    Py_TYPE(self)->tp_free(self);
  }

  static PyObject * new_(
      PyTypeObject * type, PyObject * args, PyObject * kwargs) {
    auto self = reinterpret_cast<BackendState *>(type->tp_alloc(type, 0));
    if (self == nullptr)
      return nullptr;

    // Placement new
    self = new (self) BackendState;
    return reinterpret_cast<PyObject *>(self);
  }

  static PyObject * pickle_(BackendState * self) {
    try {
      py_ref py_global = BackendState::convert_py(self->globals);
      py_ref py_locals = BackendState::convert_py(self->locals);
      py_ref py_use_local_globals =
          BackendState::convert_py(self->use_thread_local_globals);

      return py_make_tuple(py_global, py_locals, py_use_local_globals)
          .release();
    } catch (std::runtime_error &) {
      return nullptr;
    }
  }

  static PyObject * unpickle_(PyObject * cls, PyObject * args) {
    try {
      PyObject *py_locals, *py_global;
      py_ref ref =
          py_ref::steal(Q_PyObject_Vectorcall(cls, nullptr, 0, nullptr));
      BackendState * output = reinterpret_cast<BackendState *>(ref.get());
      if (output == nullptr)
        return nullptr;

      int use_thread_local_globals;
      if (!PyArg_ParseTuple(
              args, "OOp", &py_global, &py_locals, &use_thread_local_globals))
        return nullptr;
      local_state_t locals = convert_local_state(py_locals);
      global_state_t globals = convert_global_state(py_global);

      output->locals = std::move(locals);
      output->globals = std::move(globals);
      output->use_thread_local_globals = use_thread_local_globals;

      return ref.release();
    } catch (std::invalid_argument &) {
      return nullptr;
    } catch (std::bad_alloc &) {
      PyErr_NoMemory();
      return nullptr;
    }
  }

  template <typename T, typename Convertor>
  static std::vector<T> convert_iter(
      PyObject * input, Convertor item_convertor) {
    std::vector<T> output;
    py_ref iterator = py_ref::steal(PyObject_GetIter(input));
    if (!iterator)
      throw std::invalid_argument("");

    py_ref item;
    while ((item = py_ref::steal(PyIter_Next(iterator.get())))) {
      output.push_back(item_convertor(item.get()));
    }

    if (PyErr_Occurred())
      throw std::invalid_argument("");

    return output;
  }

  template <
      typename K, typename V, typename KeyConvertor, typename ValueConvertor>
  static std::unordered_map<K, V> convert_dict(
      PyObject * input, KeyConvertor key_convertor,
      ValueConvertor value_convertor) {
    std::unordered_map<K, V> output;

    if (!PyDict_Check(input))
      throw std::invalid_argument("");

    PyObject *key, *value;
    Py_ssize_t pos = 0;

    while (PyDict_Next(input, &pos, &key, &value)) {
      output[key_convertor(key)] = value_convertor(value);
    }

    if (PyErr_Occurred())
      throw std::invalid_argument("");

    return output;
  }

  static std::string convert_domain(PyObject * input) {
    std::string output = domain_to_string(input);
    if (output.empty())
      throw std::invalid_argument("");

    return output;
  }

  static backend_options convert_backend_options(PyObject * input) {
    backend_options output;
    int coerce, only;
    PyObject * py_backend;
    if (!PyArg_ParseTuple(input, "Opp", &py_backend, &coerce, &only))
      throw std::invalid_argument("");

    if (py_backend != Py_None) {
      output.backend = py_ref::ref(py_backend);
    }
    output.coerce = coerce;
    output.only = only;

    return output;
  }

  static py_ref convert_backend(PyObject * input) { return py_ref::ref(input); }

  static local_backends convert_local_backends(PyObject * input) {
    PyObject *py_skipped, *py_preferred;
    if (!PyArg_ParseTuple(input, "OO", &py_skipped, &py_preferred))
      throw std::invalid_argument("");

    local_backends output;
    output.skipped =
        convert_iter<py_ref>(py_skipped, BackendState::convert_backend);
    output.preferred = convert_iter<backend_options>(
        py_preferred, BackendState::convert_backend_options);

    return output;
  }

  static global_backends convert_global_backends(PyObject * input) {
    PyObject *py_global, *py_registered;
    int try_global_backend_last;
    if (!PyArg_ParseTuple(
            input, "OOp", &py_global, &py_registered, &try_global_backend_last))
      throw std::invalid_argument("");

    global_backends output;
    output.global = BackendState::convert_backend_options(py_global);
    output.registered =
        convert_iter<py_ref>(py_registered, BackendState::convert_backend);
    output.try_global_backend_last = try_global_backend_last;

    return output;
  }

  static global_state_t convert_global_state(PyObject * input) {
    return convert_dict<std::string, global_backends>(
        input, BackendState::convert_domain,
        BackendState::convert_global_backends);
  }

  static local_state_t convert_local_state(PyObject * input) {
    return convert_dict<std::string, local_backends>(
        input, BackendState::convert_domain,
        BackendState::convert_local_backends);
  }

  static py_ref convert_py(py_ref input) { return input; }

  static py_ref convert_py(bool input) { return py_bool(input); }

  static py_ref convert_py(backend_options input) {
    if (!input.backend) {
      input.backend = py_ref::ref(Py_None);
    }
    py_ref output = py_make_tuple(
        input.backend, py_bool(input.coerce), py_bool(input.only));
    if (!output)
      throw std::runtime_error("");
    return output;
  }

  static py_ref convert_py(const std::string & input) {
    py_ref output =
        py_ref::steal(PyUnicode_FromStringAndSize(input.c_str(), input.size()));
    if (!output)
      throw std::runtime_error("");
    return output;
  }

  template <typename T>
  static py_ref convert_py(const std::vector<T> & input) {
    py_ref output = py_ref::steal(PyList_New(input.size()));

    if (!output)
      throw std::runtime_error("");

    for (size_t i = 0; i < input.size(); i++) {
      py_ref element = convert_py(input[i]);
      PyList_SET_ITEM(output.get(), i, element.release());
    }

    return output;
  }

  static py_ref convert_py(const local_backends & input) {
    py_ref py_skipped = BackendState::convert_py(input.skipped);
    py_ref py_preferred = BackendState::convert_py(input.preferred);
    py_ref output = py_make_tuple(py_skipped, py_preferred);

    if (!output)
      throw std::runtime_error("");

    return output;
  }

  static py_ref convert_py(const global_backends & input) {
    py_ref py_globals = BackendState::convert_py(input.global);
    py_ref py_registered = BackendState::convert_py(input.registered);
    py_ref output = py_make_tuple(
        py_globals, py_registered, py_bool(input.try_global_backend_last));

    if (!output)
      throw std::runtime_error("");

    return output;
  }

  template <typename K, typename V>
  static py_ref convert_py(const std::unordered_map<K, V> & input) {
    py_ref output = py_ref::steal(PyDict_New());

    if (!output)
      throw std::runtime_error("");

    for (const auto & kv : input) {
      py_ref py_key = convert_py(kv.first);
      py_ref py_value = convert_py(kv.second);

      if (PyDict_SetItem(output.get(), py_key.get(), py_value.get()) < 0) {
        throw std::runtime_error("");
      }
    }

    return output;
  }
};

/** Clean up global python references when the module is finalized. */
void globals_free(void * /* self */) {
  global_domain_map->clear();
  BackendNotImplementedError.reset();
  identifiers.clear();
}

/** Allow GC to break reference cycles between the module and global backends.
 *
 * NOTE: local state and "thread local globals" can't be visited because we
 * can't access locals from another thread. However, since those are only
 * set through context managers they should always be unset before module
 * cleanup.
 */
int globals_traverse(PyObject * self, visitproc visit, void * arg) {
  for (const auto & kv : *global_domain_map) {
    const auto & globals = kv.second;
    PyObject * backend = globals.global.backend.get();
    Py_VISIT(backend);
    for (const auto & reg : globals.registered) {
      backend = reg.get();
      Py_VISIT(backend);
    }
  }
  return 0;
}

int globals_clear(PyObject * /* self */) {
  global_domain_map->clear();
  return 0;
}

PyObject * set_global_backend(PyObject * /* self */, PyObject * args) {
  PyObject * backend;
  int only = false, coerce = false, try_last = false;
  if (!PyArg_ParseTuple(args, "O|ppp", &backend, &coerce, &only, &try_last))
    return nullptr;

  if (!backend_validate_ua_domain(backend)) {
    return nullptr;
  }

  const auto res =
      backend_for_each_domain_string(backend, [&](const std::string & domain) {
        backend_options options;
        options.backend = py_ref::ref(backend);
        options.coerce = coerce;
        options.only = only;

        auto & domain_globals = (*current_global_state)[domain];
        domain_globals.global = options;
        domain_globals.try_global_backend_last = try_last;
        return LoopReturn::Continue;
      });

  if (res == LoopReturn::Error)
    return nullptr;

  Py_RETURN_NONE;
}

PyObject * register_backend(PyObject * /* self */, PyObject * args) {
  PyObject * backend;
  if (!PyArg_ParseTuple(args, "O", &backend))
    return nullptr;

  if (!backend_validate_ua_domain(backend)) {
    return nullptr;
  }

  const auto ret =
      backend_for_each_domain_string(backend, [&](const std::string & domain) {
        (*current_global_state)[domain].registered.push_back(
            py_ref::ref(backend));
        return LoopReturn::Continue;
      });
  if (ret == LoopReturn::Error)
    return nullptr;

  Py_RETURN_NONE;
}

void clear_single(const std::string & domain, bool registered, bool global) {
  auto domain_globals = current_global_state->find(domain);
  if (domain_globals == current_global_state->end())
    return;

  if (registered && global) {
    current_global_state->erase(domain_globals);
    return;
  }

  if (registered) {
    domain_globals->second.registered.clear();
  }

  if (global) {
    domain_globals->second.global.backend.reset();
    domain_globals->second.try_global_backend_last = false;
  }
}

PyObject * clear_backends(PyObject * /* self */, PyObject * args) {
  PyObject * domain = nullptr;
  int registered = true, global = false;
  if (!PyArg_ParseTuple(args, "O|pp", &domain, &registered, &global))
    return nullptr;

  if (domain == Py_None && registered && global) {
    current_global_state->clear();
    Py_RETURN_NONE;
  }

  auto domain_str = domain_to_string(domain);
  clear_single(domain_str, registered, global);
  Py_RETURN_NONE;
}

/** Common functionality of set_backend and skip_backend */
template <typename T>
class context_helper {
public:
  using BackendLists = SmallDynamicArray<std::vector<T> *>;
  // using BackendLists = std::vector<std::vector<T> *>;
private:
  T new_backend_;
  BackendLists backend_lists_;

public:
  const T & get_backend() const { return new_backend_; }

  context_helper() {}

  bool init(BackendLists && backend_lists, T new_backend) {
    static_assert(std::is_nothrow_move_assignable<BackendLists>::value, "");
    backend_lists_ = std::move(backend_lists);
    new_backend_ = std::move(new_backend);
    return true;
  }

  bool init(std::vector<T> & backends, T new_backend) {
    try {
      backend_lists_ = BackendLists(1, &backends);
    } catch (std::bad_alloc &) {
      PyErr_NoMemory();
      return false;
    }
    new_backend_ = std::move(new_backend);
    return true;
  }

  bool enter() {
    auto first = backend_lists_.begin();
    auto last = backend_lists_.end();
    auto cur = first;
    try {
      for (; cur < last; ++cur) {
        (*cur)->push_back(new_backend_);
      }
    } catch (std::bad_alloc &) {
      for (; first < cur; ++first) {
        (*first)->pop_back();
      }
      PyErr_NoMemory();
      return false;
    }
    return true;
  }

  bool exit() {
    bool success = true;

    for (auto * backends : backend_lists_) {
      if (backends->empty()) {
        PyErr_SetString(
            PyExc_SystemExit, "__exit__ call has no matching __enter__");
        success = false;
        continue;
      }

      if (backends->back() != new_backend_) {
        PyErr_SetString(
            PyExc_RuntimeError,
            "Found invalid context state while in __exit__. "
            "__enter__ and __exit__ may be unmatched");
        success = false;
      }

      backends->pop_back();
    }

    return success;
  }
};


struct SetBackendContext {
  PyObject_HEAD

  context_helper<backend_options> ctx_;

  static void dealloc(SetBackendContext * self) {
    PyObject_GC_UnTrack(self);
    self->~SetBackendContext();
    Py_TYPE(self)->tp_free(self);
  }

  static PyObject * new_(
      PyTypeObject * type, PyObject * args, PyObject * kwargs) {
    auto self = reinterpret_cast<SetBackendContext *>(type->tp_alloc(type, 0));
    if (self == nullptr)
      return nullptr;

    // Placement new
    self = new (self) SetBackendContext;
    return reinterpret_cast<PyObject *>(self);
  }

  static int init(
      SetBackendContext * self, PyObject * args, PyObject * kwargs) {
    static const char * kwlist[] = {"backend", "coerce", "only", nullptr};
    PyObject * backend = nullptr;
    int coerce = false;
    int only = false;

    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "O|pp", (char **)kwlist, &backend, &coerce, &only))
      return -1;

    if (!backend_validate_ua_domain(backend)) {
      return -1;
    }

    auto num_domains = backend_get_num_domains(backend);
    if (num_domains < 0) {
      return -1;
    }

    try {
      decltype(ctx_)::BackendLists backend_lists(num_domains);
      int idx = 0;

      const auto ret = backend_for_each_domain_string(
          backend, [&](const std::string & domain) {
            backend_lists[idx] = &local_domain_map[domain].preferred;
            ++idx;
            return LoopReturn::Continue;
          });

      if (ret == LoopReturn::Error) {
        return -1;
      }

      backend_options opt;
      opt.backend = py_ref::ref(backend);
      opt.coerce = coerce;
      opt.only = only;

      if (!self->ctx_.init(std::move(backend_lists), opt)) {
        return -1;
      }
    } catch (std::bad_alloc &) {
      PyErr_NoMemory();
      return -1;
    }

    return 0;
  }

  static PyObject * enter__(SetBackendContext * self, PyObject * /* args */) {
    if (!self->ctx_.enter())
      return nullptr;
    Py_RETURN_NONE;
  }

  static PyObject * exit__(SetBackendContext * self, PyObject * /*args*/) {
    if (!self->ctx_.exit())
      return nullptr;
    Py_RETURN_NONE;
  }

  static int traverse(SetBackendContext * self, visitproc visit, void * arg) {
    Py_VISIT(self->ctx_.get_backend().backend.get());
    return 0;
  }

  static PyObject * pickle_(SetBackendContext * self, PyObject * /*args*/) {
    const backend_options & opt = self->ctx_.get_backend();
    return py_make_tuple(opt.backend, py_bool(opt.coerce), py_bool(opt.only))
        .release();
  }
};


struct SkipBackendContext {
  PyObject_HEAD

  context_helper<py_ref> ctx_;

  static void dealloc(SkipBackendContext * self) {
    PyObject_GC_UnTrack(self);
    self->~SkipBackendContext();
    Py_TYPE(self)->tp_free(self);
  }

  static PyObject * new_(
      PyTypeObject * type, PyObject * args, PyObject * kwargs) {
    auto self = reinterpret_cast<SkipBackendContext *>(type->tp_alloc(type, 0));
    if (self == nullptr)
      return nullptr;

    // Placement new
    self = new (self) SkipBackendContext;
    return reinterpret_cast<PyObject *>(self);
  }

  static int init(
      SkipBackendContext * self, PyObject * args, PyObject * kwargs) {
    static const char * kwlist[] = {"backend", nullptr};
    PyObject * backend;

    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "O", (char **)kwlist, &backend))
      return -1;

    if (!backend_validate_ua_domain(backend)) {
      return -1;
    }

    auto num_domains = backend_get_num_domains(backend);
    if (num_domains < 0) {
      return -1;
    }

    try {
      decltype(ctx_)::BackendLists backend_lists(num_domains);
      int idx = 0;

      const auto ret = backend_for_each_domain_string(
          backend, [&](const std::string & domain) {
            backend_lists[idx] = &local_domain_map[domain].skipped;
            ++idx;
            return LoopReturn::Continue;
          });

      if (ret == LoopReturn::Error) {
        return -1;
      }

      if (!self->ctx_.init(std::move(backend_lists), py_ref::ref(backend))) {
        return -1;
      }
    } catch (std::bad_alloc &) {
      PyErr_NoMemory();
      return -1;
    }

    return 0;
  }

  static PyObject * enter__(SkipBackendContext * self, PyObject * /* args */) {
    if (!self->ctx_.enter())
      return nullptr;
    Py_RETURN_NONE;
  }

  static PyObject * exit__(SkipBackendContext * self, PyObject * /*args*/) {
    if (!self->ctx_.exit())
      return nullptr;
    Py_RETURN_NONE;
  }

  static int traverse(SkipBackendContext * self, visitproc visit, void * arg) {
    Py_VISIT(self->ctx_.get_backend().get());
    return 0;
  }

  static PyObject * pickle_(SkipBackendContext * self, PyObject * /*args*/) {
    return py_make_tuple(self->ctx_.get_backend()).release();
  }
};

const local_backends & get_local_backends(const std::string & domain_key) {
  static const local_backends null_local_backends;
  auto itr = local_domain_map.find(domain_key);
  if (itr == local_domain_map.end()) {
    return null_local_backends;
  }
  return itr->second;
}


const global_backends & get_global_backends(const std::string & domain_key) {
  static const global_backends null_global_backends;
  const auto & cur_globals = *current_global_state;
  auto itr = cur_globals.find(domain_key);
  if (itr == cur_globals.end()) {
    return null_global_backends;
  }
  return itr->second;
}

template <typename Callback>
LoopReturn for_each_backend_in_domain(
    const std::string & domain_key, Callback call) {
  const local_backends & locals = get_local_backends(domain_key);

  auto & skip = locals.skipped;
  auto & pref = locals.preferred;

  auto should_skip = [&](PyObject * backend) -> int {
    bool success = true;
    auto it = std::find_if(skip.begin(), skip.end(), [&](const py_ref & be) {
      auto result = PyObject_RichCompareBool(be.get(), backend, Py_EQ);
      success = (result >= 0);
      return (result != 0);
    });

    if (!success) {
      return -1;
    }

    return (it != skip.end());
  };

  LoopReturn ret = LoopReturn::Continue;
  for (int i = pref.size() - 1; i >= 0; --i) {
    auto options = pref[i];
    int skip_current = should_skip(options.backend.get());
    if (skip_current < 0)
      return LoopReturn::Error;
    if (skip_current)
      continue;

    ret = call(options.backend.get(), options.coerce);
    if (ret != LoopReturn::Continue)
      return ret;

    if (options.only || options.coerce)
      return LoopReturn::Break;
  }

  auto & globals = get_global_backends(domain_key);
  auto try_global_backend = [&] {
    auto & options = globals.global;
    if (!options.backend)
      return LoopReturn::Continue;

    int skip_current = should_skip(options.backend.get());
    if (skip_current < 0)
      return LoopReturn::Error;
    if (skip_current > 0)
      return LoopReturn::Continue;

    return call(options.backend.get(), options.coerce);
  };

  if (!globals.try_global_backend_last) {
    ret = try_global_backend();
    if (ret != LoopReturn::Continue)
      return ret;

    if (globals.global.only || globals.global.coerce)
      return LoopReturn::Break;
  }

  for (size_t i = 0; i < globals.registered.size(); ++i) {
    py_ref backend = globals.registered[i];
    int skip_current = should_skip(backend.get());
    if (skip_current < 0)
      return LoopReturn::Error;
    if (skip_current)
      continue;

    ret = call(backend.get(), false);
    if (ret != LoopReturn::Continue)
      return ret;
  }

  if (!globals.try_global_backend_last) {
    return ret;
  }
  return try_global_backend();
}

template <typename Callback>
LoopReturn for_each_backend(std::string domain, Callback call) {
  do {
    auto ret = for_each_backend_in_domain(domain, call);
    if (ret != LoopReturn::Continue) {
      return ret;
    }

    auto dot_pos = domain.rfind('.');
    if (dot_pos == std::string::npos) {
      return ret;
    }

    domain.resize(dot_pos);
  } while (!domain.empty());
  return LoopReturn::Continue;
}

struct py_func_args {
  py_ref args, kwargs;
};

struct Function {
  PyObject_HEAD
  py_ref extractor_, replacer_;  // functions to handle dispatchables
  std::string domain_key_;       // associated __ua_domain__ in UTF8
  py_ref def_args_, def_kwargs_; // default arguments
  py_ref def_impl_;              // default implementation
  py_ref dict_;                  // __dict__

  PyObject * call(PyObject * args, PyObject * kwargs);

  py_func_args replace_dispatchables(
      PyObject * backend, PyObject * args, PyObject * kwargs,
      PyObject * coerce);

  py_ref canonicalize_args(PyObject * args);
  py_ref canonicalize_kwargs(PyObject * kwargs);

  static void dealloc(Function * self) {
    PyObject_GC_UnTrack(self);
    self->~Function();
    Py_TYPE(self)->tp_free(self);
  }

  static PyObject * new_(
      PyTypeObject * type, PyObject * args, PyObject * kwargs) {
    auto self = reinterpret_cast<Function *>(type->tp_alloc(type, 0));
    if (self == nullptr)
      return nullptr;

    // Placement new
    self = new (self) Function;
    return reinterpret_cast<PyObject *>(self);
  }

  static int init(Function * self, PyObject * args, PyObject * /*kwargs*/) {
    PyObject *extractor, *replacer;
    PyObject * domain;
    PyObject *def_args, *def_kwargs;
    PyObject * def_impl;

    if (!PyArg_ParseTuple(
            args, "OOO!O!O!O", &extractor, &replacer, &PyUnicode_Type, &domain,
            &PyTuple_Type, &def_args, &PyDict_Type, &def_kwargs, &def_impl)) {
      return -1;
    }

    if (!PyCallable_Check(extractor) ||
        (replacer != Py_None && !PyCallable_Check(replacer))) {
      PyErr_SetString(
          PyExc_TypeError, "Argument extractor and replacer must be callable");
      return -1;
    }

    if (def_impl != Py_None && !PyCallable_Check(def_impl)) {
      PyErr_SetString(
          PyExc_TypeError, "Default implementation must be Callable or None");
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
  static PyObject * get_extractor(Function * self);
  static PyObject * get_replacer(Function * self);
  static PyObject * get_domain(Function * self);
  static PyObject * get_default(Function * self);
};


bool is_default(PyObject * value, PyObject * def) {
  // TODO: richer comparison for builtin types? (if cheap)
  return (value == def);
}


py_ref Function::canonicalize_args(PyObject * args) {
  const auto arg_size = PyTuple_GET_SIZE(args);
  const auto def_size = PyTuple_GET_SIZE(def_args_.get());

  if (arg_size > def_size)
    return py_ref::ref(args);

  Py_ssize_t mismatch = 0;
  for (Py_ssize_t i = arg_size - 1; i >= 0; --i) {
    auto val = PyTuple_GET_ITEM(args, i);
    auto def = PyTuple_GET_ITEM(def_args_.get(), i);
    if (!is_default(val, def)) {
      mismatch = i + 1;
      break;
    }
  }

  return py_ref::steal(PyTuple_GetSlice(args, 0, mismatch));
}


py_ref Function::canonicalize_kwargs(PyObject * kwargs) {
  if (kwargs == nullptr)
    return py_ref::steal(PyDict_New());

  PyObject *key, *def_value;
  Py_ssize_t pos = 0;
  while (PyDict_Next(def_kwargs_.get(), &pos, &key, &def_value)) {
    auto val = PyDict_GetItem(kwargs, key);
    if (val && is_default(val, def_value)) {
      PyDict_DelItem(kwargs, key);
    }
  }
  return py_ref::ref(kwargs);
}


py_func_args Function::replace_dispatchables(
    PyObject * backend, PyObject * args, PyObject * kwargs, PyObject * coerce) {
  auto has_ua_convert = PyObject_HasAttr(backend, identifiers.ua_convert->get());
  if (!has_ua_convert) {
    return {py_ref::ref(args), py_ref::ref(kwargs)};
  }

  auto dispatchables =
      py_ref::steal(PyObject_Call(extractor_.get(), args, kwargs));
  if (!dispatchables)
    return {};

  PyObject * convert_args[] = {backend, dispatchables.get(), coerce};
  auto res = py_ref::steal(Q_PyObject_VectorcallMethod(
      identifiers.ua_convert->get(), convert_args,
      array_size(convert_args) | Q_PY_VECTORCALL_ARGUMENTS_OFFSET, nullptr));
  if (!res) {
    return {};
  }

  if (res == Py_NotImplemented) {
    return {std::move(res), nullptr};
  }

  auto replaced_args = py_ref::steal(PySequence_Tuple(res.get()));
  if (!replaced_args)
    return {};

  PyObject * replacer_args[] = {nullptr, args, kwargs, replaced_args.get()};
  res = py_ref::steal(Q_PyObject_Vectorcall(
      replacer_.get(), &replacer_args[1],
      (array_size(replacer_args) - 1) | Q_PY_VECTORCALL_ARGUMENTS_OFFSET,
      nullptr));
  if (!res)
    return {};

  if (!PyTuple_Check(res.get()) || PyTuple_Size(res.get()) != 2) {
    PyErr_SetString(
        PyExc_TypeError,
        "Argument replacer must return a 2-tuple (args, kwargs)");
    return {};
  }

  auto new_args = py_ref::ref(PyTuple_GET_ITEM(res.get(), 0));
  auto new_kwargs = py_ref::ref(PyTuple_GET_ITEM(res.get(), 1));

  new_kwargs = canonicalize_kwargs(new_kwargs.get());

  if (!PyTuple_Check(new_args.get()) || !PyDict_Check(new_kwargs.get())) {
    PyErr_SetString(PyExc_ValueError, "Invalid return from argument_replacer");
    return {};
  }

  return {std::move(new_args), std::move(new_kwargs)};
}


PyObject * Function_call(Function * self, PyObject * args, PyObject * kwargs) {
  return self->call(args, kwargs);
}

class py_errinf {
  py_ref type_, value_, traceback_;

public:
  static py_errinf fetch() {
    PyObject *type, *value, *traceback;
    PyErr_Fetch(&type, &value, &traceback);
    py_errinf err;
    err.set(type, value, traceback);
    return err;
  }

  py_ref get_exception() {
    normalize();
    return value_;
  }

private:
  void set(PyObject * type, PyObject * value, PyObject * traceback) {
    type_ = py_ref::steal(type);
    value_ = py_ref::steal(value);
    traceback_ = py_ref::steal(traceback);
  }

  void normalize() {
    auto type = type_.release();
    auto value = value_.release();
    auto traceback = value_.release();
    PyErr_NormalizeException(&type, &value, &traceback);
    if (traceback) {
      PyException_SetTraceback(value, traceback);
    }
    set(type, value, traceback);
  }
};


PyObject * Function::call(PyObject * args_, PyObject * kwargs_) {
  auto args = canonicalize_args(args_);
  auto kwargs = canonicalize_kwargs(kwargs_);

  py_ref result;
  std::vector<std::pair<py_ref, py_errinf>> errors;


  auto ret =
      for_each_backend(domain_key_, [&, this](PyObject * backend, bool coerce) {
        auto new_args = replace_dispatchables(
            backend, args.get(), kwargs.get(), coerce ? Py_True : Py_False);
        if (new_args.args == Py_NotImplemented)
          return LoopReturn::Continue;
        if (new_args.args == nullptr)
          return LoopReturn::Error;

        PyObject * args[] = {
            backend, reinterpret_cast<PyObject *>(this), new_args.args.get(),
            new_args.kwargs.get()};
        result = py_ref::steal(Q_PyObject_VectorcallMethod(
            identifiers.ua_function->get(), args,
            array_size(args) | Q_PY_VECTORCALL_ARGUMENTS_OFFSET, nullptr));

        // raise BackendNotImplemeted is equivalent to return NotImplemented
        if (!result &&
            PyErr_ExceptionMatches(BackendNotImplementedError.get())) {
          errors.push_back({py_ref::ref(backend), py_errinf::fetch()});
          result = py_ref::ref(Py_NotImplemented);
        }

        // Try the default with this backend
        if (result == Py_NotImplemented && def_impl_ != Py_None) {
          backend_options opt;
          opt.backend = py_ref::ref(backend);
          opt.coerce = coerce;
          opt.only = true;
          context_helper<backend_options> ctx;
          try {
            if (!ctx.init(
                    local_domain_map[domain_key_].preferred, std::move(opt)))
              return LoopReturn::Error;
          } catch (std::bad_alloc &) {
            PyErr_NoMemory();
            return LoopReturn::Error;
          }

          if (!ctx.enter())
            return LoopReturn::Error;

          result = py_ref::steal(PyObject_Call(
              def_impl_.get(), new_args.args.get(), new_args.kwargs.get()));

          if (PyErr_Occurred() &&
              PyErr_ExceptionMatches(BackendNotImplementedError.get())) {
            errors.push_back({py_ref::ref(backend), py_errinf::fetch()});
            result = py_ref::ref(Py_NotImplemented);
          }

          if (!ctx.exit())
            return LoopReturn::Error;
        }

        if (!result)
          return LoopReturn::Error;

        if (result == Py_NotImplemented)
          return LoopReturn::Continue;

        return LoopReturn::Break; // Backend called successfully
      });

  if (ret == LoopReturn::Error)
    return nullptr;

  if (result && result != Py_NotImplemented)
    return result.release();

  // Last resort, try calling default implementation directly
  // Only call if no backend was marked only or coerce
  if (ret == LoopReturn::Continue && def_impl_ != Py_None) {
    result =
        py_ref::steal(PyObject_Call(def_impl_.get(), args.get(), kwargs.get()));
    if (!result) {
      if (!PyErr_ExceptionMatches(BackendNotImplementedError.get()))
        return nullptr;

      errors.push_back({py_ref::ref(Py_None), py_errinf::fetch()});
      result = py_ref::ref(Py_NotImplemented);
    } else if (result != Py_NotImplemented)
      return result.release();
  }


  // All backends and defaults failed, construct the exception
  auto exception_tuple = py_ref::steal(PyTuple_New(errors.size() + 1));
  PyTuple_SET_ITEM(
      exception_tuple.get(), 0,
      PyUnicode_FromString(
          "No selected backends had an implementation for this function."));
  for (size_t i = 0; i < errors.size(); ++i) {
    auto pair =
        py_make_tuple(errors[i].first, errors[i].second.get_exception());
    if (!pair)
      return nullptr;

    PyTuple_SET_ITEM(exception_tuple.get(), i + 1, pair.release());
  }
  PyErr_SetObject(BackendNotImplementedError.get(), exception_tuple.get());
  return nullptr;
}


PyObject * Function::repr(Function * self) {
  if (self->dict_)
    if (auto name = PyDict_GetItemString(self->dict_.get(), "__name__"))
      return PyUnicode_FromFormat("<uarray multimethod '%S'>", name);

  return PyUnicode_FromString("<uarray multimethod>");
}


/** Implements the descriptor protocol to allow binding to class instances */
PyObject * Function::descr_get(
    PyObject * self, PyObject * obj, PyObject * type) {
  if (!obj) {
    Py_INCREF(self);
    return self;
  }

  return PyMethod_New(self, obj);
}


/** Make members visible to the garbage collector */
int Function::traverse(Function * self, visitproc visit, void * arg) {
  Py_VISIT(self->extractor_.get());
  Py_VISIT(self->replacer_.get());
  Py_VISIT(self->def_args_.get());
  Py_VISIT(self->def_kwargs_.get());
  Py_VISIT(self->def_impl_.get());
  Py_VISIT(self->dict_.get());
  return 0;
}


/** Break reference cycles when being GCed */
int Function::clear(Function * self) {
  self->extractor_.reset();
  self->replacer_.reset();
  self->def_args_.reset();
  self->def_kwargs_.reset();
  self->def_impl_.reset();
  self->dict_.reset();
  return 0;
}

PyObject * Function::get_extractor(Function * self) {
  Py_INCREF(self->extractor_.get());
  return self->extractor_.get();
}

PyObject * Function::get_replacer(Function * self) {
  Py_INCREF(self->replacer_.get());
  return self->replacer_.get();
}

PyObject * Function::get_default(Function * self) {
  Py_INCREF(self->def_impl_.get());
  return self->def_impl_.get();
}

PyObject * Function::get_domain(Function * self) {
  return PyUnicode_FromStringAndSize(
      self->domain_key_.c_str(), self->domain_key_.size());
}


PyMethodDef BackendState_Methods[] = {
    {"_pickle", (PyCFunction)BackendState::pickle_, METH_NOARGS, nullptr},
    {"_unpickle", (PyCFunction)BackendState::unpickle_,
     METH_VARARGS | METH_CLASS, nullptr},
    {NULL} /* Sentinel */
};

PyTypeObject BackendStateType = {
    PyVarObject_HEAD_INIT(NULL, 0)     /* boilerplate */
    "uarray._BackendState",            /* tp_name */
    sizeof(BackendState),              /* tp_basicsize */
    0,                                 /* tp_itemsize */
    (destructor)BackendState::dealloc, /* tp_dealloc */
    0,                                 /* tp_print */
    0,                                 /* tp_getattr */
    0,                                 /* tp_setattr */
    0,                                 /* tp_reserved */
    0,                                 /* tp_repr */
    0,                                 /* tp_as_number */
    0,                                 /* tp_as_sequence */
    0,                                 /* tp_as_mapping */
    0,                                 /* tp_hash  */
    0,                                 /* tp_call */
    0,                                 /* tp_str */
    0,                                 /* tp_getattro */
    0,                                 /* tp_setattro */
    0,                                 /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                /* tp_flags */
    0,                                 /* tp_doc */
    0,                                 /* tp_traverse */
    0,                                 /* tp_clear */
    0,                                 /* tp_richcompare */
    0,                                 /* tp_weaklistoffset */
    0,                                 /* tp_iter */
    0,                                 /* tp_iternext */
    BackendState_Methods,              /* tp_methods */
    0,                                 /* tp_members */
    0,                                 /* tp_getset */
    0,                                 /* tp_base */
    0,                                 /* tp_dict */
    0,                                 /* tp_descr_get */
    0,                                 /* tp_descr_set */
    0,                                 /* tp_dictoffset */
    0,                                 /* tp_init */
    0,                                 /* tp_alloc */
    BackendState::new_,                /* tp_new */
};

PyObject * get_state(PyObject * /* self */, PyObject * /* args */) {
  py_ref ref = py_ref::steal(Q_PyObject_Vectorcall(
      reinterpret_cast<PyObject *>(&BackendStateType), nullptr, 0, nullptr));
  BackendState * output = reinterpret_cast<BackendState *>(ref.get());

  output->locals = local_domain_map;
  output->use_thread_local_globals =
      (current_global_state != global_domain_map.get());
  output->globals = *current_global_state;

  return ref.release();
}

PyObject * set_state(PyObject * /* self */, PyObject * args) {
  PyObject * arg;
  int reset_allowed = false;
  if (!PyArg_ParseTuple(args, "O|p", &arg, &reset_allowed))
    return nullptr;

  if (!PyObject_IsInstance(
          arg, reinterpret_cast<PyObject *>(&BackendStateType))) {
    PyErr_SetString(
        PyExc_TypeError, "state must be a uarray._BackendState object.");
    return nullptr;
  }

  BackendState * state = reinterpret_cast<BackendState *>(arg);
  local_domain_map = state->locals;
  bool use_thread_local_globals =
      (!reset_allowed) || state->use_thread_local_globals;
  current_global_state =
      use_thread_local_globals ? &thread_local_domain_map : global_domain_map.get();

  if (use_thread_local_globals)
    thread_local_domain_map = state->globals;
  else
    thread_local_domain_map.clear();


  Py_RETURN_NONE;
}

PyObject * determine_backend(PyObject * /*self*/, PyObject * args) {
  PyObject *domain_object, *dispatchables;
  int coerce;
  if (!PyArg_ParseTuple(
          args, "OOp:determine_backend", &domain_object, &dispatchables,
          &coerce))
    return nullptr;

  auto domain = domain_to_string(domain_object);
  if (domain.empty())
    return nullptr;

  auto dispatchables_tuple = py_ref::steal(PySequence_Tuple(dispatchables));
  if (!dispatchables_tuple)
    return nullptr;

  py_ref selected_backend;
  auto result = for_each_backend_in_domain(
      domain, [&](PyObject * backend, bool coerce_backend) {
        auto has_ua_convert =
            PyObject_HasAttr(backend, identifiers.ua_convert->get());

        if (!has_ua_convert) {
          // If no __ua_convert__, assume it won't accept the type
          return LoopReturn::Continue;
        }

        PyObject * convert_args[] = {
            backend, dispatchables_tuple.get(),
            (coerce && coerce_backend) ? Py_True : Py_False};

        auto res = py_ref::steal(Q_PyObject_VectorcallMethod(
            identifiers.ua_convert->get(), convert_args,
            array_size(convert_args) | Q_PY_VECTORCALL_ARGUMENTS_OFFSET,
            nullptr));
        if (!res) {
          return LoopReturn::Error;
        }

        if (res == Py_NotImplemented) {
          return LoopReturn::Continue;
        }

        // __ua_convert__ succeeded, so select this backend
        selected_backend = py_ref::ref(backend);
        return LoopReturn::Break;
      });

  if (result != LoopReturn::Continue)
    return selected_backend.release();

  // All backends failed, raise an error
  PyErr_SetString(
      BackendNotImplementedError.get(),
      "No backends could accept input of this type.");
  return nullptr;
}


// getset takes mutable char * in python < 3.7
static char dict__[] = "__dict__";
static char arg_extractor[] = "arg_extractor";
static char arg_replacer[] = "arg_replacer";
static char default_[] = "default";
static char domain[] = "domain";
PyGetSetDef Function_getset[] = {
    {dict__, PyObject_GenericGetDict, PyObject_GenericSetDict},
    {arg_extractor, (getter)Function::get_extractor, NULL},
    {arg_replacer, (getter)Function::get_replacer, NULL},
    {default_, (getter)Function::get_default, NULL},
    {domain, (getter)Function::get_domain, NULL},
    {NULL} /* Sentinel */
};

PyTypeObject FunctionType = {
    PyVarObject_HEAD_INIT(NULL, 0) /* boilerplate */
    /* tp_name= */ "uarray._Function",
    /* tp_basicsize= */ sizeof(Function),
    /* tp_itemsize= */ 0,
    /* tp_dealloc= */ (destructor)Function::dealloc,
    /* tp_print= */ 0,
    /* tp_getattr= */ 0,
    /* tp_setattr= */ 0,
    /* tp_reserved= */ 0,
    /* tp_repr= */ (reprfunc)Function::repr,
    /* tp_as_number= */ 0,
    /* tp_as_sequence= */ 0,
    /* tp_as_mapping= */ 0,
    /* tp_hash= */ 0,
    /* tp_call= */ (ternaryfunc)Function_call,
    /* tp_str= */ 0,
    /* tp_getattro= */ PyObject_GenericGetAttr,
    /* tp_setattro= */ PyObject_GenericSetAttr,
    /* tp_as_buffer= */ 0,
    /* tp_flags= */
    (Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC | Q_Py_TPFLAGS_METHOD_DESCRIPTOR),
    /* tp_doc= */ 0,
    /* tp_traverse= */ (traverseproc)Function::traverse,
    /* tp_clear= */ (inquiry)Function::clear,
    /* tp_richcompare= */ 0,
    /* tp_weaklistoffset= */ 0,
    /* tp_iter= */ 0,
    /* tp_iternext= */ 0,
    /* tp_methods= */ 0,
    /* tp_members= */ 0,
    /* tp_getset= */ Function_getset,
    /* tp_base= */ 0,
    /* tp_dict= */ 0,
    /* tp_descr_get= */ Function::descr_get,
    /* tp_descr_set= */ 0,
    /* tp_dictoffset= */ offsetof(Function, dict_),
    /* tp_init= */ (initproc)Function::init,
    /* tp_alloc= */ 0,
    /* tp_new= */ Function::new_,
};


PyMethodDef SetBackendContext_Methods[] = {
    {"__enter__", (PyCFunction)SetBackendContext::enter__, METH_NOARGS,
     nullptr},
    {"__exit__", (PyCFunction)SetBackendContext::exit__, METH_VARARGS, nullptr},
    {"_pickle", (PyCFunction)SetBackendContext::pickle_, METH_NOARGS, nullptr},
    {NULL} /* Sentinel */
};

PyTypeObject SetBackendContextType = {
    PyVarObject_HEAD_INIT(NULL, 0)             /* boilerplate */
    "uarray._SetBackendContext",               /* tp_name */
    sizeof(SetBackendContext),                 /* tp_basicsize */
    0,                                         /* tp_itemsize */
    (destructor)SetBackendContext::dealloc,    /* tp_dealloc */
    0,                                         /* tp_print */
    0,                                         /* tp_getattr */
    0,                                         /* tp_setattr */
    0,                                         /* tp_reserved */
    0,                                         /* tp_repr */
    0,                                         /* tp_as_number */
    0,                                         /* tp_as_sequence */
    0,                                         /* tp_as_mapping */
    0,                                         /* tp_hash  */
    0,                                         /* tp_call */
    0,                                         /* tp_str */
    0,                                         /* tp_getattro */
    0,                                         /* tp_setattro */
    0,                                         /* tp_as_buffer */
    (Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC), /* tp_flags */
    0,                                         /* tp_doc */
    (traverseproc)SetBackendContext::traverse, /* tp_traverse */
    0,                                         /* tp_clear */
    0,                                         /* tp_richcompare */
    0,                                         /* tp_weaklistoffset */
    0,                                         /* tp_iter */
    0,                                         /* tp_iternext */
    SetBackendContext_Methods,                 /* tp_methods */
    0,                                         /* tp_members */
    0,                                         /* tp_getset */
    0,                                         /* tp_base */
    0,                                         /* tp_dict */
    0,                                         /* tp_descr_get */
    0,                                         /* tp_descr_set */
    0,                                         /* tp_dictoffset */
    (initproc)SetBackendContext::init,         /* tp_init */
    0,                                         /* tp_alloc */
    SetBackendContext::new_,                   /* tp_new */
};


PyMethodDef SkipBackendContext_Methods[] = {
    {"__enter__", (PyCFunction)SkipBackendContext::enter__, METH_NOARGS,
     nullptr},
    {"__exit__", (PyCFunction)SkipBackendContext::exit__, METH_VARARGS,
     nullptr},
    {"_pickle", (PyCFunction)SkipBackendContext::pickle_, METH_NOARGS, nullptr},
    {NULL} /* Sentinel */
};

PyTypeObject SkipBackendContextType = {
    PyVarObject_HEAD_INIT(NULL, 0)              /* boilerplate */
    "uarray._SkipBackendContext",               /* tp_name */
    sizeof(SkipBackendContext),                 /* tp_basicsize */
    0,                                          /* tp_itemsize */
    (destructor)SkipBackendContext::dealloc,    /* tp_dealloc */
    0,                                          /* tp_print */
    0,                                          /* tp_getattr */
    0,                                          /* tp_setattr */
    0,                                          /* tp_reserved */
    0,                                          /* tp_repr */
    0,                                          /* tp_as_number */
    0,                                          /* tp_as_sequence */
    0,                                          /* tp_as_mapping */
    0,                                          /* tp_hash  */
    0,                                          /* tp_call */
    0,                                          /* tp_str */
    0,                                          /* tp_getattro */
    0,                                          /* tp_setattro */
    0,                                          /* tp_as_buffer */
    (Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC),  /* tp_flags */
    0,                                          /* tp_doc */
    (traverseproc)SkipBackendContext::traverse, /* tp_traverse */
    0,                                          /* tp_clear */
    0,                                          /* tp_richcompare */
    0,                                          /* tp_weaklistoffset */
    0,                                          /* tp_iter */
    0,                                          /* tp_iternext */
    SkipBackendContext_Methods,                 /* tp_methods */
    0,                                          /* tp_members */
    0,                                          /* tp_getset */
    0,                                          /* tp_base */
    0,                                          /* tp_dict */
    0,                                          /* tp_descr_get */
    0,                                          /* tp_descr_set */
    0,                                          /* tp_dictoffset */
    (initproc)SkipBackendContext::init,         /* tp_init */
    0,                                          /* tp_alloc */
    SkipBackendContext::new_,                   /* tp_new */
};


PyMethodDef method_defs[] = {
    {"set_global_backend", set_global_backend, METH_VARARGS, nullptr},
    {"register_backend", register_backend, METH_VARARGS, nullptr},
    {"clear_backends", clear_backends, METH_VARARGS, nullptr},
    {"determine_backend", determine_backend, METH_VARARGS, nullptr},
    {"get_state", get_state, METH_NOARGS, nullptr},
    {"set_state", set_state, METH_VARARGS, nullptr},
    {NULL} /* Sentinel */
};

PyModuleDef uarray_module = {
    PyModuleDef_HEAD_INIT,
    /* m_name= */ "uarray._uarray",
    /* m_doc= */ nullptr,
    /* m_size= */ -1,
    /* m_methods= */ method_defs,
    /* m_slots= */ nullptr,
    /* m_traverse= */ globals_traverse,
    /* m_clear= */ globals_clear,
    /* m_free= */ globals_free};

} // namespace


PyMODINIT_FUNC PyInit__uarray(void) {

  auto m = py_ref::steal(PyModule_Create(&uarray_module));
  if (!m)
    return nullptr;

  if (PyType_Ready(&FunctionType) < 0)
    return nullptr;
  Py_INCREF(&FunctionType);
  PyModule_AddObject(m.get(), "_Function", (PyObject *)&FunctionType);

  if (PyType_Ready(&SetBackendContextType) < 0)
    return nullptr;
  Py_INCREF(&SetBackendContextType);
  PyModule_AddObject(
      m.get(), "_SetBackendContext", (PyObject *)&SetBackendContextType);

  if (PyType_Ready(&SkipBackendContextType) < 0)
    return nullptr;
  Py_INCREF(&SkipBackendContextType);
  PyModule_AddObject(
      m.get(), "_SkipBackendContext", (PyObject *)&SkipBackendContextType);

  if (PyType_Ready(&BackendStateType) < 0)
    return nullptr;
  Py_INCREF(&BackendStateType);
  PyModule_AddObject(m.get(), "_BackendState", (PyObject *)&BackendStateType);

  BackendNotImplementedError = py_ref::steal(PyErr_NewExceptionWithDoc(
      "uarray.BackendNotImplementedError",
      "An exception that is thrown when no compatible"
      " backend is found for a method.",
      PyExc_NotImplementedError, nullptr));
  if (!BackendNotImplementedError)
    return nullptr;
  Py_INCREF(BackendNotImplementedError.get());
  PyModule_AddObject(
      m.get(), "BackendNotImplementedError", BackendNotImplementedError.get());

  if (!identifiers.init())
    return nullptr;

#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(m.get(), Py_MOD_GIL_NOT_USED);
#endif

  return m.release();
}
