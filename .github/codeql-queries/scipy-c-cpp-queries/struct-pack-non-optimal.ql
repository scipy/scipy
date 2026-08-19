/**
 * @name Suboptimal type definition
 * @description Highlights structs whose members are not laid out optimally, in the sense
 *    that by reordering them one could reduce the amount of internal padding on a 64-bit architecture.
 * @kind problem
 * @id cpp/suboptimal-64-bit-type
 * @problem.severity recommendation
 * @tags efficiency
 */

import semmle.code.cpp.padding.Padding

from
  PaddedType t, Architecture arch, WideCharType wc, int holes, int size, int percentage, int optimum
where
  arch.pointerSize() = 64 and // Select 64-bit architecture
  arch.wideCharSize() = (wc.getSize() * 8) and // Select Windows(sizeof(wchar_t == 2)) or non-Windows(sizeof(wchar_t == 4))
  t.isPrecise() and
  optimum = t.optimalSize(arch) and
  size = arch.paddedSize(t) and
  holes = size - optimum and
  holes > 0 and
  percentage = (holes * 100.0 / size.(float)).ceil()
select t,
  t.getName() + " could be optimized to save " + holes + "/" + t.wastedSpace(arch) +
    " bits of padding (or " + percentage + "% of its size)."
