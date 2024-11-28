#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <map>
#include <mutex>

#include "lp_data/HighsOptions.h"

namespace py = pybind11;

class HighsOptionsManager {
public:
  HighsOptionsManager() {
    initialize_log_options();
    for (const auto &record : highs_options_.records) {
      record_type_lookup_.emplace(record->name, record->type);
    }
  }

  const HighsOptions &get_highs_options() const { return highs_options_; }

  const std::map<std::string, HighsOptionType> &get_record_type_lookup() const {
    return record_type_lookup_;
  }

  template <typename OptionRecordType, typename T>
  bool check_option(const std::string &name, const T value) {
    std::lock_guard<std::mutex> guard(highs_options_mutex);
    HighsInt idx = 0;
    const OptionStatus idx_status = getOptionIndex(
        highs_log_options, name.c_str(), highs_options_.records, idx);

    if (OptionStatus::kOk != idx_status) {
      return false;
    }

    OptionRecordType &record =
        static_cast<OptionRecordType &>(*highs_options_.records.at(idx));
    const OptionStatus check_status =
        checkOptionValue(highs_log_options, record, value);
    if (OptionStatus::kIllegalValue == check_status) {
      return false;
    }

    return true;
  }

private:
  HighsOptions highs_options_;
  std::mutex highs_options_mutex;
  std::map<std::string, HighsOptionType> record_type_lookup_;
  HighsLogOptions highs_log_options;

  static constexpr bool log_to_console = false;
  static constexpr bool output_flag = true;

  void initialize_log_options() {
    highs_log_options.log_stream = nullptr;
    highs_log_options.output_flag = const_cast<bool *>(&output_flag);
    highs_log_options.log_to_console = const_cast<bool *>(&log_to_console);
    highs_log_options.log_dev_level = nullptr;
    highs_log_options.user_log_callback = nullptr;
    highs_log_options.user_log_callback_data = nullptr;
    highs_log_options.user_callback_data = nullptr;
    highs_log_options.user_callback_active = false;
  }
};

PYBIND11_MODULE(_highs_options, m, py::mod_gil_not_used()) {
  py::class_<HighsOptionsManager>(m, "HighsOptionsManager")
      .def(py::init<>())
      .def("get_option_type",
           [](const HighsOptionsManager &manager, const std::string &name) {
             const auto &lookup = manager.get_record_type_lookup().find(name);
             if (manager.get_record_type_lookup().end() == lookup) {
               return -1;
             }
             return static_cast<int>(lookup->second);
           })
      .def("get_all_option_types", &HighsOptionsManager::get_record_type_lookup)
      .def("get_highs_options_records",
           [](const HighsOptionsManager &manager) {
             std::vector<std::string> records_names;
             for (const auto &record : manager.get_highs_options().records) {
               records_names.push_back(record->name);
             }
             return records_names;
           })
      .def("check_int_option",
           [](HighsOptionsManager &self, const std::string &name, int value) {
             try {
               return self.check_option<OptionRecordInt, int>(name, value);
             } catch (const std::exception &e) {
               py::print("Exception caught in check_int_option:", e.what());
               return false;
             }
           })
      .def(
          "check_double_option",
          [](HighsOptionsManager &self, const std::string &name, double value) {
            try {
              return self.check_option<OptionRecordDouble, double>(name, value);
            } catch (const std::exception &e) {
              py::print("Exception caught in check_double_option:", e.what());
              return false;
            }
          })
      .def("check_string_option",
           [](HighsOptionsManager &self, const std::string &name,
              const std::string &value) {
             try {
               return self.check_option<OptionRecordString, std::string>(name,
                                                                         value);
             } catch (const std::exception &e) {
               py::print("Exception caught in check_string_option:", e.what());
               return false;
             }
           });
}
