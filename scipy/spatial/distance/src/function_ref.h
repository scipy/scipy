#pragma once
#include <type_traits>

// Type-erased reference to a function object.
template <class Signature>
class FunctionRef;

template <class Ret, class... Args>
class FunctionRef<Ret(Args...)> {
public:
    template <class FunctionObject,
        typename std::enable_if<  // Don't disable the default copy constructor
            !std::is_same<typename std::decay<FunctionObject>::type, FunctionRef>::value,
        int>::type = 0
    >
    FunctionRef(FunctionObject && a_FunctionObject) {
        data_ = &a_FunctionObject;
        call_function_ = &ObjectFunctionCaller<FunctionObject>;
    }

    Ret operator () (Args... args) {
        return call_function_(data_, std::forward<Args>(args)...);
    }

private:

    template <class ObjectType>
    static Ret ObjectFunctionCaller(void * callable, Args...  args) {
        // Convert opaque reference to the concrete type.
        using ObjectPtr = typename std::add_pointer<ObjectType>::type;
        auto & Object = *static_cast<ObjectPtr>(callable);

        // Forward the call down to the object.
        return Object(std::forward<Args>(args)...);
    }

    using CallFunction = Ret(*)(void *, Args...);

    void* data_;
    CallFunction call_function_;
};
