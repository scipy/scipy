.. _deciding-on-new-features:

========================
Deciding on new features
========================
The general decision rule to accept a proposed new feature has so far
been conditional on:

1. The method is applicable in many fields and "generally agreed" to
   be useful,
2. It fits the topic of the submodule, and does not require extensive
   support frameworks to operate,
3. The implementation looks sound and unlikely to need much tweaking in
   the future (e.g., limited expected maintenance burden),
4. Someone wants to contribute it, and
5. Someone wants to review it.

The last criterion is often a sticking point for proposed features. Code cannot
be merged until it has been thoroughly reviewed, and there is always a backlog
of maintenance tasks that compete for reviewers' time. Ideally, contributors
should line up a reviewer with suitable domain expertise before beginning
work.

Although it's difficult to give hard rules on what "generally useful
and generally agreed to work" means, it may help to weigh the following
against each other:

- Is the method used/useful in different domains in practice?
  How much domain-specific background knowledge is needed to use it
  properly?
- Consider the code already in the module.  Is what you are adding
  an omission?  Does it solve a problem that you'd expect the module
  be able to solve?  Does it supplement an existing feature in
  a significant way?
- Consider the equivalence class of similar methods / features usually
  expected. Among them, what would in principle be the minimal set so
  that there's not a glaring omission in the offered features remaining?
  How much stuff would that be? Does including a representative one of
  them cover most use cases? Would it in principle sound reasonable to
  include everything from the minimal set in the module?
- Is what you are adding something that is well understood in the
  literature? If not, how sure are you that it will turn out well?
  Does the method perform well compared to other similar ones?
- Note that the twice-a-year release cycle and backward-compatibility
  policy makes correcting things later on more difficult.

The scopes of the submodules also vary, so it's probably best to consider
each as if it's a separate project - "numerical evaluation of special
functions" is relatively well-defined, but "commonly needed optimization
algorithms" less so.
