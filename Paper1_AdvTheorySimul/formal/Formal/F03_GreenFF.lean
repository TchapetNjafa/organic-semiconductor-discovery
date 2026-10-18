import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity

namespace Formal.F03

noncomputable def FF (v : ℝ) : ℝ :=
  (v - Real.log (v + 0.72)) / (v + 1)

theorem FF_bounds {v : ℝ} (hv : 0 ≤ v) : 0 < FF v ∧ FF v < 1 := by
  have hy : 0 < v + (0.72 : ℝ) := by linarith
  have hd : 0 < v + 1 := by linarith
  have hu := Real.log_le_sub_one_of_pos hy
  have hl := Real.one_sub_inv_le_log_of_pos hy
  have hi : (v + (0.72 : ℝ))⁻¹ < 2 := by
    rw [inv_eq_one_div, div_lt_iff₀ hy]
    linarith
  constructor
  · exact div_pos (by linarith) hd
  · unfold FF
    apply (div_lt_one hd).2
    linarith

/-- The quotient rule, with its derivative kept in a form suitable for bounds. -/
theorem hasDerivAt_FF {v : ℝ} (hv : 0 ≤ v) :
    HasDerivAt FF
      (((1 - 1 / (v + 0.72)) * (v + 1) - (v - Real.log (v + 0.72))) /
        (v + 1) ^ 2) v := by
  have hy : v + (0.72 : ℝ) ≠ 0 := by linarith
  have hd : v + (1 : ℝ) ≠ 0 := by linarith
  unfold FF
  simpa only [id_eq, mul_one] using
    ((hasDerivAt_id v).fun_sub (((hasDerivAt_id v).add_const (0.72 : ℝ)).log hy)).fun_div
      ((hasDerivAt_id v).add_const (1 : ℝ)) hd

theorem FF_strictMonoOn : StrictMonoOn FF (Set.Ici 1) := by
  apply strictMonoOn_of_deriv_pos (convex_Ici 1)
  · intro v hv
    exact (hasDerivAt_FF (by have := hv; simp only [Set.mem_Ici] at this; linarith)).continuousAt.continuousWithinAt
  · intro v hv
    have hv1 : 1 < v := by simpa only [interior_Ici, Set.mem_Ioi] using hv
    have hv0 : 0 ≤ v := by linarith
    have hy : 0 < v + (0.72 : ℝ) := by linarith
    have hl := Real.one_sub_inv_le_log_of_pos hy
    have hi : 1 / (v + (0.72 : ℝ)) < (3 / 4 : ℝ) := by
      apply (div_lt_iff₀ hy).2
      linarith
    have hmul : (1 / (v + (0.72 : ℝ))) * (v + 0.72) = 1 :=
      div_mul_cancel₀ 1 (ne_of_gt hy)
    rw [(hasDerivAt_FF hv0).deriv]
    apply div_pos _ (sq_pos_of_pos (by linarith))
    rw [inv_eq_one_div] at hl
    nlinarith

end Formal.F03
