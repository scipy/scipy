/* Analyze RTL for GNU compiler.
   Copyright (C) 2020-2022 Free Software Foundation, Inc.

This file is part of GCC.

GCC is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 3, or (at your option) any later
version.

GCC is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with GCC; see the file COPYING3.  If not see
<http://www.gnu.org/licenses/>.  */

/* Note that for historical reasons, many rtlanal.cc functions are
   declared in rtl.h rather than here.  */

#ifndef GCC_RTLANAL_H
#define GCC_RTLANAL_H

/* A dummy register value that represents the whole of variable memory.
   Using ~0U means that arrays that track both registers and memory can
   be indexed by regno + 1.  */
const unsigned int MEM_REGNO = ~0U;

/* Bitmasks of flags describing an rtx_obj_reference.  See the accessors
   in the class for details.  */
namespace rtx_obj_flags
{
  const uint16_t IS_READ = 1U << 0;
  const uint16_t IS_WRITE = 1U << 1;
  const uint16_t IS_CLOBBER = 1U << 2;
  const uint16_t IS_PRE_POST_MODIFY = 1U << 3;
  const uint16_t IS_MULTIREG = 1U << 4;
  const uint16_t IN_MEM_LOAD = 1U << 5;
  const uint16_t IN_MEM_STORE = 1U << 6;
  const uint16_t IN_SUBREG = 1U << 7;
  const uint16_t IN_NOTE = 1U << 8;

  /* Flags that apply to all subrtxes of the rtx they were originally
     added for.  */
  static const uint16_t STICKY_FLAGS = IN_NOTE;
}

/* Contains information about a reference to a register or variable memory.  */
class rtx_obj_reference
{
public:
  rtx_obj_reference () = default;
  rtx_obj_reference (unsigned int regno, uint16_t flags,
		     machine_mode mode, unsigned int multireg_offset = 0);

  bool is_reg () const { return regno != MEM_REGNO; }
  bool is_mem () const { return regno == MEM_REGNO; }

  /* True if the reference is a read or a write respectively.
     Both flags are set in a read-modify-write context, such as
     for read_modify_subreg_p.  */
  bool is_read () const { return flags & rtx_obj_flags::IS_READ; }
  bool is_write () const { return flags & rtx_obj_flags::IS_WRITE; }

  /* True if IS_WRITE and if the write is a clobber rather than a set.  */
  bool is_clobber () const { return flags & rtx_obj_flags::IS_CLOBBER; }

  /* True if the reference is updated by an RTX_AUTOINC.  Both IS_READ
     and IS_WRITE are also true if so.  */
  bool is_pre_post_modify () const
  {
    return flags & rtx_obj_flags::IS_PRE_POST_MODIFY;
  }

  /* True if the register is part of a multi-register hard REG.  */
  bool is_multireg () const { return flags & rtx_obj_flags::IS_MULTIREG; }

  /* True if the reference occurs in the address of a load MEM.  */
  bool in_mem_load () const { return flags & rtx_obj_flags::IN_MEM_LOAD; }

  /* True if the reference occurs in the address of a store MEM.  */
  bool in_mem_store () const { return flags & rtx_obj_flags::IN_MEM_STORE; }

  /* True if the reference occurs in any kind of MEM address.  */
  bool in_address () const { return in_mem_load () || in_mem_store (); }

  /* True if the reference occurs in a SUBREG.  */
  bool in_subreg () const { return flags & rtx_obj_flags::IN_SUBREG; }

  /* True if the reference occurs in a REG_EQUAL or REG_EQUIV note.  */
  bool in_note () const { return flags & rtx_obj_flags::IN_NOTE; }

  /* The referenced register, or MEM_REGNO for variable memory.  */
  unsigned int regno;

  /* A bitmask of rtx_obj_flags.  */
  unsigned int flags : 16;

  /* The mode of the reference.  If IS_MULTIREG, this is the mode of
     REGNO - MULTIREG_OFFSET.  */
  machine_mode mode : 8;

  /* If IS_MULTIREG, the offset of REGNO from the start of the register.  */
  unsigned int multireg_offset : 8;
};

/* Construct a reference with the given fields.  */

inline rtx_obj_reference::rtx_obj_reference (unsigned int regno, uint16_t flags,
					     machine_mode mode,
					     unsigned int multireg_offset)
  : regno (regno),
    flags (flags),
    mode (mode),
    multireg_offset (multireg_offset)
{
}

/* Contains information about an rtx or an instruction, including a
   list of rtx_obj_references.  The storage backing the list needs
   to be filled in by assigning to REF_BEGIN and REF_END.  */

class rtx_properties
{
public:
  rtx_properties ();

  void try_to_add_reg (const_rtx x, unsigned int flags = 0);
  void try_to_add_dest (const_rtx x, unsigned int flags = 0);
  void try_to_add_src (const_rtx x, unsigned int flags = 0);
  void try_to_add_pattern (const_rtx pat);
  void try_to_add_note (const_rtx x);
  void try_to_add_insn (const rtx_insn *insn, bool include_notes);

  iterator_range<rtx_obj_reference *> refs () const;

  /* Return the number of rtx_obj_references that have been recorded.  */
  size_t num_refs () const { return ref_iter - ref_begin; }

  bool has_side_effects () const;

  /* [REF_BEGIN, REF_END) is the maximum extent of the memory available
     for recording references.  REG_ITER is the first unused entry.  */
  rtx_obj_reference *ref_begin;
  rtx_obj_reference *ref_iter;
  rtx_obj_reference *ref_end;

  /* True if the rtx includes an asm.  */
  unsigned int has_asm : 1;

  /* True if the rtx includes a call.  */
  unsigned int has_call : 1;

  /* True if the rtx includes an RTX_AUTOINC expression.  */
  unsigned int has_pre_post_modify : 1;

  /* True if the rtx contains volatile references, in the sense of
     volatile_refs_p.  */
  unsigned int has_volatile_refs : 1;

  /* For future expansion.  */
  unsigned int spare : 28;
};

inline rtx_properties::rtx_properties ()
  : ref_begin (nullptr),
    ref_iter (nullptr),
    ref_end (nullptr),
    has_asm (false),
    has_call (false),
    has_pre_post_modify (false),
    has_volatile_refs (false),
    spare (0)
{
}

/* Like add_src, but treat X has being part of a REG_EQUAL or
   REG_EQUIV note.  */

inline void
rtx_properties::try_to_add_note (const_rtx x)
{
  try_to_add_src (x, rtx_obj_flags::IN_NOTE);
}

/* Return true if the rtx has side effects, in the sense of
   side_effects_p (except for side_effects_p's special handling
   of combine.cc clobbers).  */

inline bool
rtx_properties::has_side_effects () const
{
  return has_volatile_refs || has_pre_post_modify || has_call;
}

/* Return an iterator range for all the references, suitable for
   range-based for loops.  */

inline iterator_range<rtx_obj_reference *>
rtx_properties::refs () const
{
  return { ref_begin, ref_iter };
}

/* BASE is derived from rtx_properties and provides backing storage
   for REF_BEGIN.  It has a grow () method that increases the amount
   of memory available if the initial allocation was too small.  */

template<typename Base>
class growing_rtx_properties : public Base
{
public:
  template<typename... Args>
  growing_rtx_properties (Args...);

  template<typename AddFn>
  void repeat (AddFn add);

  /* Wrappers around the try_to_* functions that always succeed.  */
  void add_dest (const_rtx x, unsigned int flags = 0);
  void add_src (const_rtx x, unsigned int flags = 0);
  void add_pattern (const_rtx pat);
  void add_note (const_rtx x);
  void add_insn (const rtx_insn *insn, bool include_notes);
};

template<typename Base>
template<typename... Args>
growing_rtx_properties<Base>::growing_rtx_properties (Args... args)
  : Base (std::forward<Args> (args)...)
{
}

/* Perform ADD until there is enough room to hold the result.  */

template<typename Base>
template<typename AddFn>
inline void
growing_rtx_properties<Base>::repeat (AddFn add)
{
  ptrdiff_t count = this->num_refs ();
  for (;;)
    {
      add ();
      /* This retries if the storage happened to be exactly the right size,
	 but that's expected to be a rare case and so isn't worth
	 optimizing for.  */
      if (__builtin_expect (this->ref_iter != this->ref_end, 1))
	break;
      this->grow (count);
    }
}

template<typename Base>
inline void
growing_rtx_properties<Base>::add_dest (const_rtx x, unsigned int flags)
{
  repeat ([&]() { this->try_to_add_dest (x, flags); });
}

template<typename Base>
inline void
growing_rtx_properties<Base>::add_src (const_rtx x, unsigned int flags)
{
  repeat ([&]() { this->try_to_add_src (x, flags); });
}

template<typename Base>
inline void
growing_rtx_properties<Base>::add_pattern (const_rtx pat)
{
  repeat ([&]() { this->try_to_add_pattern (pat); });
}

template<typename Base>
inline void
growing_rtx_properties<Base>::add_note (const_rtx x)
{
  repeat ([&]() { this->try_to_add_note (x); });
}

template<typename Base>
inline void
growing_rtx_properties<Base>::add_insn (const rtx_insn *insn, bool include_notes)
{
  repeat ([&]() { this->try_to_add_insn (insn, include_notes); });
}

/* A base class for vec_rtx_properties; see there for details.  */

class vec_rtx_properties_base : public rtx_properties
{
  static const size_t SIZE = 32;

public:
  vec_rtx_properties_base ();
  ~vec_rtx_properties_base ();

protected:
  void grow (ptrdiff_t);

private:
  rtx_obj_reference m_storage[SIZE];
};

inline vec_rtx_properties_base::vec_rtx_properties_base ()
{
  ref_begin = ref_iter = m_storage;
  ref_end = m_storage + SIZE;
}

inline vec_rtx_properties_base::~vec_rtx_properties_base ()
{
  if (__builtin_expect (ref_begin != m_storage, 0))
    free (ref_begin);
}

/* A rtx_properties that stores its references in a temporary array.
   Like auto_vec, the array is initially on the stack, but can switch
   to the heap if necessary.

   The reason for implementing this as a derived class is that the
   default on-stack size should be enough for the vast majority of
   expressions and instructions.  It's therefore not worth paying
   the cost of conditionally calling grow code at every site that
   records a new reference.  Instead, the rtx_properties code can use
   trivial iterator updates for the common case, and in the rare case
   that the vector needs to be resized, we can pay the cost of
   collecting the references a second time.  */
using vec_rtx_properties = growing_rtx_properties<vec_rtx_properties_base>;

bool
vec_series_highpart_p (machine_mode result_mode, machine_mode op_mode,
		       rtx sel);

bool
vec_series_lowpart_p (machine_mode result_mode, machine_mode op_mode, rtx sel);

#endif
