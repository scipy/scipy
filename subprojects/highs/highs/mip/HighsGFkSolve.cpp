/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsGFkSolve.h"

#include <cassert>

#include "util/HighsSplay.h"

void HighsGFkSolve::link(HighsInt pos) {
  Anext[pos] = colhead[Acol[pos]];
  Aprev[pos] = -1;
  colhead[Acol[pos]] = pos;
  if (Anext[pos] != -1) Aprev[Anext[pos]] = pos;

  ++colsize[Acol[pos]];

  auto get_row_left = [&](HighsInt pos) -> HighsInt& { return ARleft[pos]; };
  auto get_row_right = [&](HighsInt pos) -> HighsInt& { return ARright[pos]; };
  auto get_row_key = [&](HighsInt pos) { return Acol[pos]; };
  highs_splay_link(pos, rowroot[Arow[pos]], get_row_left, get_row_right,
                   get_row_key);
  ++rowsize[Arow[pos]];
}

void HighsGFkSolve::unlink(HighsInt pos) {
  HighsInt next = Anext[pos];
  HighsInt prev = Aprev[pos];

  if (next != -1) Aprev[next] = prev;

  if (prev != -1)
    Anext[prev] = next;
  else
    colhead[Acol[pos]] = next;
  --colsize[Acol[pos]];

  auto get_row_left = [&](HighsInt pos) -> HighsInt& { return ARleft[pos]; };
  auto get_row_right = [&](HighsInt pos) -> HighsInt& { return ARright[pos]; };
  auto get_row_key = [&](HighsInt pos) { return Acol[pos]; };
  highs_splay_unlink(pos, rowroot[Arow[pos]], get_row_left, get_row_right,
                     get_row_key);
  --rowsize[Arow[pos]];

  Avalue[pos] = 0;
  freeslots.push(pos);
}

void HighsGFkSolve::storeRowPositions(HighsInt pos) {
  if (pos == -1) return;

  assert(iterstack.empty());
  iterstack.push_back(pos);
  do {
    pos = iterstack.back();
    iterstack.pop_back();

    rowpositions.push_back(pos);
    rowposColsizes.push_back(colsize[Acol[pos]]);

    if (ARleft[pos] != -1) iterstack.push_back(ARleft[pos]);
    if (ARright[pos] != -1) iterstack.push_back(ARright[pos]);
  } while (!iterstack.empty());
}

HighsInt HighsGFkSolve::findNonzero(HighsInt row, HighsInt col) {
  if (rowroot[row] == -1) return -1;

  auto get_row_left = [&](HighsInt pos) -> HighsInt& { return ARleft[pos]; };
  auto get_row_right = [&](HighsInt pos) -> HighsInt& { return ARright[pos]; };
  auto get_row_key = [&](HighsInt pos) { return Acol[pos]; };
  rowroot[row] =
      highs_splay(col, rowroot[row], get_row_left, get_row_right, get_row_key);

  if (Acol[rowroot[row]] == col) return rowroot[row];

  return -1;
}

void HighsGFkSolve::addNonzero(HighsInt row, HighsInt col, unsigned int val) {
  assert(findNonzero(row, col) == -1);
  HighsInt pos;
  if (freeslots.empty()) {
    pos = Avalue.size();
    Avalue.push_back(val);
    Arow.push_back(row);
    Acol.push_back(col);
    Anext.push_back(-1);
    Aprev.push_back(-1);
    ARleft.push_back(-1);
    ARright.push_back(-1);
  } else {
    pos = freeslots.top();
    freeslots.pop();
    Avalue[pos] = val;
    Arow[pos] = row;
    Acol[pos] = col;
    Aprev[pos] = -1;
  }

  link(pos);
}
