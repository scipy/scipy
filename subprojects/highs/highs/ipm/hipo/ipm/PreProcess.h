#ifndef HIPO_PRE_POST_PROCESS
#define HIPO_PRE_POST_PROCESS

#include <cassert>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

class Model;
struct Iterate;

struct PreprocessorPoint {
  std::vector<double>& x;
  std::vector<double>& xl;
  std::vector<double>& xu;
  std::vector<double>& slack;
  std::vector<double>& y;
  std::vector<double>& zl;
  std::vector<double>& zu;

  void assertConsistency(Int n, Int m) const;
};

struct PreprocessAction {
  Int n_pre, m_pre, n_post, m_post;

  virtual ~PreprocessAction() = default;
  virtual void apply(Model& model) = 0;
  virtual void undo(PreprocessorPoint& point, const Model& model,
                    const Iterate& it) const = 0;
  virtual void print(std::stringstream& stream) const = 0;
};

struct PreprocessEmptyRows : public PreprocessAction {
  std::vector<Int> rows_shift;
  std::vector<double> b_pre;
  Int empty_rows{};

  void apply(Model& model) override;
  void undo(PreprocessorPoint& point, const Model& model,
            const Iterate& it) const override;
  void print(std::stringstream& stream) const override;
};

struct PreprocessFixedVars : public PreprocessAction {
  Int fixed_vars{};
  std::vector<double> fixed_at;

  // information about the columns that get removed
  struct FixedVarsData {
    double c;
    std::vector<Int> indA, indQ;
    std::vector<double> valA, valQ;
  };
  std::map<Int, FixedVarsData> data;

  void apply(Model& model) override;
  void undo(PreprocessorPoint& point, const Model& model,
            const Iterate& it) const override;
  void print(std::stringstream& stream) const override;
};

struct PreprocessScaling : public PreprocessAction {
  bool scaled = false;

  void apply(Model& model) override;
  void undo(PreprocessorPoint& point, const Model& model,
            const Iterate& it) const override;
  void print(std::stringstream& stream) const override;
};

struct PreprocessFormulation : public PreprocessAction {
  void apply(Model& model) override;
  void undo(PreprocessorPoint& point, const Model& model,
            const Iterate& it) const override;
  void print(std::stringstream& stream) const override;
};

struct PreprocessFreeVars : public PreprocessAction {
  Int free_vars_count{};

  void apply(Model& model) override;
  void undo(PreprocessorPoint& point, const Model& model,
            const Iterate& it) const override;
  void print(std::stringstream& stream) const override;
};

struct Preprocessor {
  std::vector<std::unique_ptr<PreprocessAction>> stack;

  void apply(Model& model);
  void undo(PreprocessorPoint& point, const Model& model,
            const Iterate& it) const;
  void print(std::stringstream& log_stream) const;
};

}  // namespace hipo

#endif
