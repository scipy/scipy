/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../extern/pdqsort/pdqsort.h"
#include "presolve/HPresolve.h"
#include "presolve/HPresolveUtils.h"

namespace presolve {

HPresolve::Result HPresolve::strengthenInequalities(
    HighsPostsolveStack& postsolve_stack, HighsInt& num_strengthened) {
  std::vector<int8_t> complementation;
  std::vector<HighsCDouble> reducedcost;
  std::vector<HighsCDouble> upper;
  std::vector<HighsInt> indices;
  std::vector<HighsInt> positions;
  std::vector<HighsInt> stack;
  std::vector<double> coefs;
  std::vector<HighsInt> cover;

  num_strengthened = 0;
  // Check for timeout according to this frequency
  const HighsInt check_time_frequency = 100;

  for (HighsInt row = 0; row != model->num_row_; ++row) {
    if (rowsize[row] <= 1) continue;
    if (isRanged(row)) continue;

    // do not run on very dense rows as this could get expensive
    HighsInt rowsize_limit =
        std::max(HighsInt{1000}, (model->num_col_ - numDeletedCols) / 20);
    if (rowsize[row] > rowsize_limit) continue;

    // printf("strengthening knapsack of %" HIGHSINT_FORMAT " vars\n",
    // rowsize[row]);

    HighsCDouble maxviolation;
    HighsCDouble continuouscontribution = 0.0;
    // maximum contribution of negligible terms not used in covers; a term is
    // negligible if its maximum contribution is at most primal_feastol
    HighsCDouble smallcontribution = 0.0;
    double scale;

    if (model->row_lower_[row] != -kHighsInf) {
      // ax >= lb, treat as -ax <= -lb --> -ax + lb <= 0
      maxviolation = model->row_lower_[row];
      scale = -1.0;
    } else {
      // ax <= ub --> ax - ub <= 0
      maxviolation = -model->row_upper_[row];
      scale = 1.0;
    }

    complementation.clear();
    reducedcost.clear();
    upper.clear();
    indices.clear();
    positions.clear();
    complementation.reserve(rowsize[row]);
    reducedcost.reserve(rowsize[row]);
    upper.reserve(rowsize[row]);
    indices.reserve(rowsize[row]);
    positions.reserve(rowsize[row]);
    stack.reserve(rowsize[row]);
    stack.push_back(rowroot[row]);

    bool skiprow = false;

    while (!stack.empty()) {
      // pop element from stack
      HighsInt pos = stack.back();
      stack.pop_back();

      // add non-zeros to stack
      if (ARright[pos] != -1) stack.push_back(ARright[pos]);
      if (ARleft[pos] != -1) stack.push_back(ARleft[pos]);

      // get column index
      HighsInt col = Acol[pos];

      // skip row if a column bound is not finite
      skiprow = model->col_lower_[col] == -kHighsInf ||
                model->col_upper_[col] == kHighsInf;
      if (skiprow) break;

      // compute maximum violation
      // scale =  1:  ax <=  ub --> violation =  ax - ub > 0
      // scale = -1: -ax <= -lb --> violation = -ax + lb > 0
      // this means that for scale = 1 we sum up an upper bound on constraint
      // activity, and for scale = -1 we sum up a lower bound on constraint
      // activity.
      int8_t comp;
      double weight = Avalue[pos] * scale;
      HighsCDouble ub = static_cast<HighsCDouble>(model->col_upper_[col]) -
                        static_cast<HighsCDouble>(model->col_lower_[col]);
      if (weight > 0) {
        comp = 1;
        maxviolation +=
            static_cast<HighsCDouble>(model->col_upper_[col]) * weight;
      } else {
        comp = -1;
        maxviolation +=
            static_cast<HighsCDouble>(model->col_lower_[col]) * weight;
        weight = -weight;
      }

      if (model->integrality_[col] == HighsVarType::kContinuous) {
        continuouscontribution += weight * ub;
        continue;
      }

      if (weight * ub <= primal_feastol) {
        smallcontribution += weight * ub;
        continue;
      }

      indices.push_back(static_cast<HighsInt>(reducedcost.size()));
      positions.push_back(pos);
      reducedcost.push_back(weight);
      complementation.push_back(comp);
      upper.push_back(ub);
    }

    // Check for timeout according to frequency, unless a particularly
    // dense row has just been analysed
    if ((row & check_time_frequency) == 0 || 10 * rowsize[row] > rowsize_limit)
      HPRESOLVE_CHECKED_CALL(checkTimeLimit());

    if (skiprow) {
      stack.clear();
      continue;
    }

    // maxviolation <= 0 implies that the constraint is redundant:
    // scale =  1: upper bound on activity <= model->row_upper_[row]
    // scale = -1: model->row_lower_[row]  <= lower bound on activity
    if (maxviolation <= primal_feastol) {
      HPRESOLVE_CHECKED_CALL(rowPresolve(postsolve_stack, row));
      continue;
    }

    const double smallVal =
        std::max(100 * primal_feastol,
                 primal_feastol * static_cast<double>(maxviolation));
    while (true) {
      if (maxviolation - continuouscontribution - smallcontribution <=
              smallVal ||
          indices.empty())
        break;

      pdqsort(indices.begin(), indices.end(), [&](HighsInt i1, HighsInt i2) {
        return std::make_pair(reducedcost[i1], i1) >
               std::make_pair(reducedcost[i2], i2);
      });

      HighsCDouble lambda =
          maxviolation - continuouscontribution - smallcontribution;

      cover.clear();
      cover.reserve(indices.size());

      for (size_t i = indices.size(); i > 0; --i) {
        HighsInt index = indices[i - 1];
        HighsCDouble delta = upper[index] * reducedcost[index];

        if (upper[index] <= 1000.0 && reducedcost[index] > smallVal &&
            lambda - delta <= smallVal)
          cover.push_back(index);
        else
          lambda -= delta;
      }

      if (cover.empty() || lambda <= smallVal) break;

      HighsInt alpos = *std::min_element(
          cover.begin(), cover.end(), [&](HighsInt i1, HighsInt i2) {
            if (reducedcost[i1] <= 1e-3 || reducedcost[i2] <= 1e-3)
              return reducedcost[i1] > reducedcost[i2];
            return reducedcost[i1] < reducedcost[i2];
          });

      HighsCDouble al = reducedcost[alpos];
      coefs.resize(cover.size());
      HighsCDouble coverrhs = max(ceil(lambda / al - primal_feastol), 1.0);
      HighsCDouble slackupper = -coverrhs;

      HighsCDouble step = kHighsInf;
      for (size_t i = 0; i != cover.size(); ++i) {
        coefs[i] =
            static_cast<double>(ceil(min(reducedcost[cover[i]], lambda) / al -
                                     options->small_matrix_value));
        slackupper += upper[cover[i]] * coefs[i];
        step = min(step, reducedcost[cover[i]] / coefs[i]);
      }
      assert(slackupper >= 0.0);
      step = min(step, maxviolation / coverrhs);
      maxviolation -= step * coverrhs;

      // only the cover members' reduced costs change, so only they can become
      // negligible; move their maximum contribution to smallcontribution and
      // set their reduced cost to zero, so that they are removed below
      for (size_t i = 0; i != cover.size(); ++i) {
        HighsInt index = cover[i];
        reducedcost[index] -= step * coefs[i];
        if (reducedcost[index] * upper[index] <= primal_feastol) {
          if (reducedcost[index] > 0.0)
            smallcontribution += reducedcost[index] * upper[index];
          reducedcost[index] = 0.0;
        }
      }

      // terms with non-positive reduced cost contribute nothing
      indices.erase(
          std::remove_if(indices.begin(), indices.end(),
                         [&](HighsInt i) { return reducedcost[i] <= 0.0; }),
          indices.end());
      if (step * slackupper <= primal_feastol) {
        smallcontribution += step * slackupper;
      } else {
        indices.push_back(static_cast<HighsInt>(reducedcost.size()));
        reducedcost.push_back(step);
        upper.push_back(slackupper);
      }
    }

    double threshold = static_cast<double>(maxviolation + primal_feastol);

    indices.erase(std::remove_if(indices.begin(), indices.end(),
                                 [&](HighsInt i) {
                                   return static_cast<size_t>(i) >=
                                              positions.size() ||
                                          abs(reducedcost[i]) <= threshold;
                                 }),
                  indices.end());
    if (indices.empty()) continue;

    auto updateNonZeros = [&](HighsInt row, HighsCDouble& rhs,
                              HighsInt direction) {
      for (HighsInt i : indices) {
        assert(Arow[positions[i]] == row);
        HighsCDouble coefdelta = direction * (reducedcost[i] - maxviolation);
        HighsInt col = Acol[positions[i]];

        if (complementation[i] == -1) {
          rhs += coefdelta * model->col_lower_[col];
          addToMatrix(row, col, static_cast<double>(coefdelta));
        } else {
          rhs -= coefdelta * model->col_upper_[col];
          addToMatrix(row, col, static_cast<double>(-coefdelta));
        }
      }
    };

    // update / add non-zeros
    if (scale < 0) {
      HighsCDouble lhs = model->row_lower_[row];
      updateNonZeros(row, lhs, HighsInt{-1});
      changeRowLower(row, static_cast<double>(lhs));
    } else {
      HighsCDouble rhs = model->row_upper_[row];
      updateNonZeros(row, rhs, HighsInt{1});
      changeRowUpper(row, static_cast<double>(rhs));
    }

    num_strengthened += indices.size();
  }

  return Result::kOk;
}

}  // namespace presolve
