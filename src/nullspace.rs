use core::convert::identity;
use core::iter;

use num_bigint::BigUint;

use crate::bitvec::BitVec;
use crate::gaussian_elimination::is_in_row_echelon_form;
use crate::util::transpose;

/// Find the solution with index `solution_index` within some enumeration of the nullspace of the
/// given matrix.
///
/// Note that `solution_index = 0` will give the trivial solution. To enumerate nontrivial
/// solutions, start with `solution_index = 1` and increment until `None` is returned.
pub(crate) fn nullspace_member(rows: &[BitVec], solution_index: &BigUint) -> Option<BitVec> {
    assert!(is_in_row_echelon_form(rows), "not RE: {:?}", rows);

    let num_vars = rows[0].len();
    rows.iter().for_each(|row| assert_eq!(row.len(), num_vars));

    // Back substitution.
    let mut assignment = BitVec::new(num_vars);
    let mut num_assigned = 0usize;

    let mut next_solution_index_bit = 0;
    let mut bits = iter::from_fn(|| {
        let b = solution_index.bit(next_solution_index_bit);
        next_solution_index_bit += 1;
        Some(b)
    });

    for r in (0..rows.len()).rev() {
        let row = &rows[r];
        if let Some(pivot_index) = row.iter().position(identity) {
            // We're going to set the pivot to whichever bool satisfies this constraint,
            // so any unset variable after the pivot can be freely chosen.
            for c in pivot_index + 1..row.len() - num_assigned {
                assignment.set(c, bits.next().unwrap());
                num_assigned += 1;
            }

            // Set the pivot to whichever bool satisfies this constraint.
            let mut target = false;
            for c in pivot_index + 1..row.len() {
                target ^= row.get(c) & assignment.get(c);
            }
            assignment.set(pivot_index, target);
            num_assigned += 1;
        }
    }

    // Set any remaining free variables.
    for c in 0..num_vars - num_assigned {
        assignment.set(c, bits.next().unwrap());
    }

    if solution_index.bits() > next_solution_index_bit {
        // The requested index exceeds the size of the nullspace. The solution we found is just a
        // repeat of a previous one, since we didn't use the highest bits of the index.
        None
    } else {
        validate(rows, &assignment);
        Some(assignment)
    }
}

fn validate(rows: &[BitVec], assignment: &BitVec) {
    let selected_cols = transpose(rows.to_vec())
        .into_iter()
        .enumerate()
        .filter_map(|(i, col)| assignment.get(i).then_some(col));

    let mut acc = BitVec::new(rows.len());
    for col in selected_cols {
        acc ^= &col;
    }
    assert_eq!(
        acc.count_ones(),
        0,
        "{:?} isn't in the nullspace of {:?}",
        assignment,
        rows
    );
}

#[cfg(test)]
mod tests {
    use alloc::vec;
    use alloc::vec::Vec;

    use num_bigint::BigUint;
    use num_integer::Integer;
    use num_traits::{One, Zero};

    use crate::bitvec::BitVec;
    use crate::nullspace::nullspace_member;

    #[test]
    fn test_nullspace_members() {
        // Solving for (x, y, z) subject to
        //     x + y = 0
        //         y = 0
        // There are two solutions: (0, 0, 0) and (0, 0, 1).
        test_nullspace_member(
            vec![
                BitVec::from([true, true, false]),
                BitVec::from([false, true, false]),
            ],
            BigUint::from(2u8),
        );

        // Solving for (x, y, z) subject to
        //     y + z = 0
        // There are four solutions: (0, 0, 0), (0, 1, 1), (1, 0, 0), (1, 1, 1).
        test_nullspace_member(vec![BitVec::from([false, true, true])], BigUint::from(4u8));

        // Solving for (x, y, z) subject to
        // x         = 0
        //     y + z = 0
        // There are two solutions: (0, 0, 0) and (0, 1, 1).
        test_nullspace_member(
            vec![
                BitVec::from([true, false, false]),
                BitVec::from([false, true, true]),
            ],
            BigUint::from(2u8),
        );

        // Sam as above, but adding a trivial 0 = 0 constraint which shouldn't change anything.
        test_nullspace_member(
            vec![
                BitVec::from([true, false, false]),
                BitVec::from([false, true, true]),
                BitVec::from([false, false, false]),
            ],
            BigUint::from(2u8),
        );

        // Solving for (x, y, z) subject to
        //     x + y     = 0
        //         y     = 0
        //             z = 0
        // There is only the trivial solution.
        test_nullspace_member(
            vec![
                BitVec::from([true, true, false]),
                BitVec::from([false, true, false]),
                BitVec::from([false, false, true]),
            ],
            BigUint::one(),
        );
    }

    fn test_nullspace_member(rows: Vec<BitVec>, num_solutions: BigUint) {
        let width = rows[0].len();

        // Every matrix should have a trivial solution.
        assert_eq!(
            nullspace_member(&rows, &BigUint::zero()),
            Some(BitVec::new(width))
        );

        let mut solution_index = BigUint::one();
        while solution_index < num_solutions {
            let assignment = nullspace_member(&rows, &solution_index).expect("No solutions found");
            assert_ne!(
                assignment.count_ones(),
                0,
                "another trivial solution was returned"
            );
            solution_index.inc();
        }

        assert_eq!(nullspace_member(&rows, &solution_index), None);
    }
}
