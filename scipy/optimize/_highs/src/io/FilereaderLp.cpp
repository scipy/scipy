/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/FilereaderLp.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "io/FilereaderLp.h"

#include <cstdarg>
#include <exception>
#include <map>

#include "../external/filereaderlp/reader.hpp"

FilereaderRetcode FilereaderLp::readModelFromFile(const HighsOptions& options,
                                                  HighsLp& model) {
  try {
    Model m = readinstance(options.model_file);

    // build variable index and gather variable information
    std::map<std::string, unsigned int> varindex;

    model.numCol_ = m.variables.size();
    model.numRow_ = m.constraints.size();
    for (unsigned int i = 0; i < m.variables.size(); i++) {
      varindex[m.variables[i]->name] = i;
      model.colLower_.push_back(m.variables[i]->lowerbound);
      model.colUpper_.push_back(m.variables[i]->upperbound);
      model.col_names_.push_back(m.variables[i]->name);
    }

    // get objective
    model.offset_ = m.objective->offset;
    model.colCost_.resize(model.numCol_, 0.0);
    for (unsigned int i = 0; i < m.objective->linterms.size(); i++) {
      std::shared_ptr<LinTerm> lt = m.objective->linterms[i];
      model.colCost_[varindex[lt->var->name]] = lt->coef;
    }

    // handle constraints
    std::map<std::shared_ptr<Variable>, std::vector<unsigned int>>
        consofvarmap_index;
    std::map<std::shared_ptr<Variable>, std::vector<double>> consofvarmap_value;
    for (unsigned int i = 0; i < m.constraints.size(); i++) {
      std::shared_ptr<Constraint> con = m.constraints[i];
      for (unsigned int j = 0; j < con->expr->linterms.size(); j++) {
        std::shared_ptr<LinTerm> lt = con->expr->linterms[j];
        if (consofvarmap_index.count(lt->var) == 0) {
          consofvarmap_index[lt->var] = std::vector<unsigned int>();
          consofvarmap_value[lt->var] = std::vector<double>();
        }
        consofvarmap_index[lt->var].push_back(i);
        consofvarmap_value[lt->var].push_back(lt->coef);
      }

      model.rowLower_.push_back(con->lowerbound);
      model.rowUpper_.push_back(con->upperbound);
    }

    int nz = 0;
    for (int i = 0; i < model.numCol_; i++) {
      std::shared_ptr<Variable> var = m.variables[i];
      model.Astart_.push_back(nz);
      for (unsigned int j = 0; j < consofvarmap_index[var].size(); j++) {
        model.Aindex_.push_back(consofvarmap_index[var][j]);
        model.Avalue_.push_back(consofvarmap_value[var][j]);
        nz++;
      }
    }
    model.Astart_.push_back(nz);
    model.sense_ = m.sense == ObjectiveSense::MIN ? ObjSense::MINIMIZE
                                                  : ObjSense::MAXIMIZE;
  } catch (std::invalid_argument& ex) {
    return FilereaderRetcode::PARSERERROR;
  }
  return FilereaderRetcode::OK;
}

void FilereaderLp::writeToFile(FILE* file, const char* format, ...) {
  va_list argptr;
  va_start(argptr, format);
  char stringbuffer[LP_MAX_LINE_LENGTH + 1];
  int tokenlength = vsprintf(stringbuffer, format, argptr);
  if (this->linelength + tokenlength >= LP_MAX_LINE_LENGTH) {
    fprintf(file, "\n");
    fprintf(file, "%s", stringbuffer);
    this->linelength = tokenlength;
  } else {
    fprintf(file, "%s", stringbuffer);
    this->linelength += tokenlength;
  }
}

void FilereaderLp::writeToFileLineend(FILE* file) {
  fprintf(file, "\n");
  this->linelength = 0;
}

HighsStatus FilereaderLp::writeModelToFile(const HighsOptions& options,
                                           const std::string filename,
                                           HighsLp& model) {
  FILE* file = fopen(filename.c_str(), "w");

  // write comment at the start of the file
  this->writeToFile(file, "\\ %s", LP_COMMENT_FILESTART);
  this->writeToFileLineend(file);

  // write objective
  this->writeToFile(file, "%s",
                    model.sense_ == ObjSense::MINIMIZE ? "min" : "max");
  this->writeToFileLineend(file);
  this->writeToFile(file, " obj: ");
  for (int i = 0; i < model.numCol_; i++) {
    this->writeToFile(file, "%+g x%d ", model.colCost_[i], (i + 1));
  }
  this->writeToFileLineend(file);

  // write constraint section, lower & upper bounds are one constraint each
  this->writeToFile(file, "st");
  this->writeToFileLineend(file);
  for (int row = 0; row < model.numRow_; row++) {
    if (model.rowLower_[row] == model.rowUpper_[row]) {
      // equality constraint
      this->writeToFile(file, " con%d: ", row + 1);
      for (int var = 0; var < model.numCol_; var++) {
        for (int idx = model.Astart_[var]; idx < model.Astart_[var + 1];
             idx++) {
          if (model.Aindex_[idx] == row) {
            this->writeToFile(file, "%+g x%d ", model.Avalue_[idx], var + 1);
          }
        }
      }
      this->writeToFile(file, "= %+g", model.rowLower_[row]);
      this->writeToFileLineend(file);
    } else {
      if (model.rowLower_[row] > -HIGHS_CONST_INF) {
        // has a lower bounds
        this->writeToFile(file, " con%dlo: ", row + 1);
        for (int var = 0; var < model.numCol_; var++) {
          for (int idx = model.Astart_[var]; idx < model.Astart_[var + 1];
               idx++) {
            if (model.Aindex_[idx] == row) {
              this->writeToFile(file, "%+g x%d ", model.Avalue_[idx], var + 1);
            }
          }
        }
        this->writeToFile(file, ">= %+g", model.rowLower_[row]);
        this->writeToFileLineend(file);
      } else if (model.rowUpper_[row] < HIGHS_CONST_INF) {
        // has an upper bounds
        this->writeToFile(file, " con%dup: ", row + 1);
        for (int var = 0; var < model.numCol_; var++) {
          for (int idx = model.Astart_[var]; idx < model.Astart_[var + 1];
               idx++) {
            if (model.Aindex_[idx] == row) {
              this->writeToFile(file, "%+g x%d ", model.Avalue_[idx], var + 1);
            }
          }
        }
        this->writeToFile(file, "<= %+g", model.rowUpper_[row]);
        this->writeToFileLineend(file);
      } else {
        // constraint has infinite lower & upper bounds so not a proper
        // constraint, does not get written
      }
    }
  }

  // write bounds section
  this->writeToFile(file, "bounds");
  this->writeToFileLineend(file);
  for (int i = 0; i < model.numCol_; i++) {
    // if both lower/upper bound are +/-infinite: [name] free
    if (model.colLower_[i] > -HIGHS_CONST_INF &&
        model.colUpper_[i] < HIGHS_CONST_INF) {
      this->writeToFile(file, " %+g <= x%d <= %+g", model.colLower_[i], i + 1,
                        model.colUpper_[i]);
      this->writeToFileLineend(file);
    } else if (model.colLower_[i] <= -HIGHS_CONST_INF &&
               model.colUpper_[i] < HIGHS_CONST_INF) {
      this->writeToFile(file, " -inf <= x%d <= %+g", i + 1, model.colUpper_[i]);
      this->writeToFileLineend(file);

    } else if (model.colLower_[i] > -HIGHS_CONST_INF &&
               model.colUpper_[i] >= HIGHS_CONST_INF) {
      this->writeToFile(file, " %+g <= x%d <= +inf", model.colLower_[i], i + 1);
      this->writeToFileLineend(file);
    } else {
      this->writeToFile(file, " x%d free", i + 1);
      this->writeToFileLineend(file);
    }
  }

  // write binary section
  this->writeToFile(file, "bin");
  this->writeToFileLineend(file);

  // write general section
  this->writeToFile(file, "gen");
  this->writeToFileLineend(file);

  // write semi section
  this->writeToFile(file, "semi");
  this->writeToFileLineend(file);

  // write end
  this->writeToFile(file, "end");
  this->writeToFileLineend(file);

  fclose(file);
  return HighsStatus::OK;
}
