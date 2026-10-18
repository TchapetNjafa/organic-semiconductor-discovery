import Mathlib.Data.Finset.Basic
import Mathlib.Data.Fintype.Basic
import Mathlib.Tactic.NormNum

namespace Formal.F02

variable {ι α : Type*}

/-- The admissible elements common to every convention in the finite index set. -/
def Robust (C : Finset ι) (Adm : ι → Finset α) : Set α :=
  {x | ∀ c ∈ C, x ∈ Adm c}

theorem mem_Robust (C : Finset ι) (Adm : ι → Finset α) (x : α) :
    x ∈ Robust C Adm ↔ ∀ c ∈ C, x ∈ Adm c := Iff.rfl

theorem Robust_subset (C : Finset ι) (Adm : ι → Finset α) {c : ι} (hc : c ∈ C) :
    Robust C Adm ⊆ (Adm c : Set α) := by
  intro x hx
  exact hx c hc

/-- Even two nonempty admissible sets can have empty robust intersection. -/
theorem Robust_empty_example :
    Robust (Finset.univ : Finset Bool)
      (fun c => if c then ({0} : Finset ℕ) else {1}) = ∅ := by
  ext x
  simp only [Robust, Set.mem_setOf_eq, Set.mem_empty_iff_false, iff_false]
  intro h
  have h0 := h true (Finset.mem_univ true)
  have h1 := h false (Finset.mem_univ false)
  norm_num at h0 h1
  exact Nat.zero_ne_one (h0.symm.trans h1)

end Formal.F02
