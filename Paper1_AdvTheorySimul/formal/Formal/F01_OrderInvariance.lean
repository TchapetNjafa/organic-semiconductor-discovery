import Mathlib.Data.Finset.Filter
import Mathlib.Data.List.Perm.Basic
import Mathlib.Tactic.Tauto

open scoped Classical

namespace Formal.F01

variable {α : Type*}

/-- Apply each predicate in the list, in order. Classical decidability also
allows lists of arbitrary predicates, hence in particular decidable ones. -/
noncomputable def successiveFilter (s : Finset α) : List (α → Prop) → Finset α
  | [] => s
  | p :: ps => successiveFilter (s.filter p) ps

theorem successiveFilter_eq_filter (s : Finset α) (ps : List (α → Prop)) :
    successiveFilter s ps = s.filter (fun x => ∀ p ∈ ps, p x) := by
  induction ps generalizing s with
  | nil => ext x; simp [successiveFilter]
  | cons p ps ih =>
    rw [successiveFilter, ih]
    ext x
    simp only [Finset.mem_filter, List.mem_cons, forall_eq_or_imp]
    tauto

theorem successiveFilter_perm (s : Finset α) {ps qs : List (α → Prop)}
    (h : ps.Perm qs) : successiveFilter s ps = successiveFilter s qs := by
  rw [successiveFilter_eq_filter, successiveFilter_eq_filter]
  ext x
  simp only [Finset.mem_filter]
  constructor
  · rintro ⟨hx, hp⟩
    exact ⟨hx, fun p hq => hp p (h.mem_iff.mpr hq)⟩
  · rintro ⟨hx, hq⟩
    exact ⟨hx, fun p hp => hq p (h.mem_iff.mp hp)⟩

end Formal.F01
