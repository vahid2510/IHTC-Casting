# Methodology

This package implements a transient 1D conduction model for the casting and the sand mold coupled by a time-dependent interfacial heat transfer coefficient h(t).

- Governing equations: transient conduction with an enthalpy method in steel over the mushy interval.
- Interface: Robin coupling `q = h(t) (T_cast - T_mold)` enforced implicitly.
- Inverse problem: regularized least squares with early-time weighting.
- IHTC model: physics-aware spike–decay–plateau `h(t)=h_inf+(h_peak-h_inf)exp(-(t-t0)/tau)`.
- Uncertainty: simple residual bootstrap to get 95% bands for h(t).

For a deeper literature context and parameter guidelines, see the references.
