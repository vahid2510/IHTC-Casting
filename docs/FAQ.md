# FAQ

**Why 1D?**  
Local identification near an instrumented spot is well approximated by a 1D normal heat path, which improves identifiability and runtime.

**Can I separate radiation from IHTC?**  
You can, but early-time conduction dominates. If you insist, add a radiative term at the interface and solve for conductive `h(t)` only.

**How do I choose sensor depths?**  
Place at least two sensors close to the interface on the casting side and one deeper anchor. See the methodology doc for sensitivity notes.
