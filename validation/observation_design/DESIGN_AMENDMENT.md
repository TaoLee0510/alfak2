# Shared-anchor composition arm, before primary accuracy comparison

Before launching the primary accuracy batch or inspecting its design-risk
comparisons, a support-connectivity audit identified a weak link in the original
composition design. At beta=0, noiseless within-culture composition changes
identify physical birth contrasts only inside that culture's initial support.
The original A and B supports are disjoint: their 15-state union has two
contrast components (rank 13). Adding the common 2.2.2 state to B joins the
same 15-state union into one component (rank 14). At the actual beta=.005,
mutation can weakly bridge the original supports, so this is a limiting
connectivity diagnostic, not an impossibility claim about the count model.

Add `composition_anchor600`: A1 is unchanged; B_anchor has .15 at 2.2.2,
.35 at 1.1.1, .35 at 3.3.3, and .15 over the same six corner-neighbor states.
Both cultures still start at 600 cells and receive t0/T counts of 150 each.
The original `composition600` arm remains a disclosed diagnostic control;
none of its generated data or results is removed. The primary three 600-cell
layouts are depth600, time600 and composition_anchor600. All results from both
composition variants are reported. This adds 12 development fits to the original
48 (60 total), without opening independent calibration/test seeds.

Each new B_anchor path uses the previously unused process seed A1_seed+3 and
the same independent B count uniforms. The original truth, A1/A2/B paths,
reference and control-based distance groups remain unchanged. This amendment
is motivated by normalized-composition invariance, not observed accuracy.
