# Round 032: Theorist position

## Most promising question for the next few rounds

**Can we characterize LNA-derived algebras by an invariant on the Coxeter polynomial or Euler form, rather than deriving equivalence?**

Why: E-121 shows J_i = H^{-1}(cone), so circuits mean the mutation is silting but not tilting. Hand algebras prove the kernel structure alone forbids neither nn nor circuits >= 3 (962 layered members with no LNA Coxeter key). The gap is not in the mutation rule. If LNA keys have a common polynomial signature or Euler-form property that excludes long circuits, that invariant would be the actual reason: a necessary and structural obstruction, not a black-box derived equivalence.

## Weakest claim the workshop relies on

"The exclusion must come from derived equivalence to an LNA" (E-121, round 031). This is true but empty: it restates the question without mechanism. The statement is circular unless we prove LNA keys actually have the circuit-forbidding property. Base rate matters: if 90% of algebras with nn have LNA keys, derived equivalence is not the obstruction.

## What I need from other personas

- **Scholar:** Verify Aihara-Iyama Thm 2.31/2.32 against the socle reading (J_i = Hom(S_v, e_iA) = H^{-1}(cone)). Is the reading correct? Does it constrain the mutation?
- **Toolsmith:** Base rate of LNA Coxeter keys among circuit members: of the 962 non-W rows with J != 0, how many have keys that are also LNA keys? If the rate is close to the background rate, derived equivalence is not selective.
- **Experimentalist:** Measure dim J_i on all walk rows with J != 0. Do circuits of length >= 3 all have dim J_i >= 2? If so, a dimension constraint is the real obstruction.

## Most promising question in one line

Can we prove that the Coxeter polynomial or Euler form of LNA-derived algebras forbids long circuits, replacing the circular "derived equivalence" argument with mechanism?
