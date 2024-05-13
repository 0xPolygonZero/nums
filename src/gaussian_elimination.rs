use crate::bitvec::BitVec;
use alloc::vec::Vec;

pub(crate) fn gaussian_elimination(rows: &mut [BitVec]) {
    let height = rows.len();

    let mut row = 0;

    while row < height {
        let next_row: usize = (row..height)
            .min_by_key(|&r| rows[r].leading_zeros())
            .unwrap();
        rows.swap(row, next_row);
        let pivot = rows[row].leading_zeros();

        for r in row + 1..height {
            if rows[r].get(pivot) {
                // We want to do
                //     rows[r] ^= &rows[row];
                // Since the borrow checker doesn't know that these are disjoint pieces of data,
                // we need split_at to borrow both rows at the same time.
                let (upper_rows, lower_rows) = rows.split_at_mut(r);
                lower_rows[0] ^= &upper_rows[row];
            }
        }

        row += 1;
    }

    debug_assert!(is_in_row_echelon_form(rows), "not RE: {:?}", rows);
}

pub(crate) fn is_in_row_echelon_form(rows: &[BitVec]) -> bool {
    let width = rows[0].len();
    let pivots: Vec<usize> = rows.iter().map(BitVec::leading_zeros).collect();
    for i in 1..pivots.len() {
        let p1 = pivots[i - 1];
        let p2 = pivots[i];
        let both_zero_rows = p1 == width && p2 == width;
        if p1 >= p2 && !both_zero_rows {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use crate::bitvec::BitVec;
    use crate::gaussian_elimination::gaussian_elimination;
    use alloc::vec;

    #[test]
    fn test_gaussian_elimination_1() {
        let mut mat = vec![
            BitVec::from([true, true, true, true]),
            BitVec::from([true, false, true, false]),
            BitVec::from([false, true, true, true]),
        ];
        gaussian_elimination(&mut mat);
        assert_eq!(
            mat,
            vec![
                BitVec::from([true, true, true, true]),
                BitVec::from([false, true, false, true]),
                BitVec::from([false, false, true, false]),
            ]
        );
    }

    #[test]
    fn test_gaussian_elimination_2() {
        let mut mat = vec![
            BitVec::from([false, false, false]),
            BitVec::from([true, false, false]),
            BitVec::from([true, true, false]),
            BitVec::from([false, true, false]),
        ];
        gaussian_elimination(&mut mat);
        assert_eq!(
            mat,
            vec![
                BitVec::from([true, false, false]),
                BitVec::from([false, true, false]),
                BitVec::from([false, false, false]),
                BitVec::from([false, false, false]),
            ]
        );
    }

    #[test]
    fn test_gaussian_elimination_3() {
        let mut mat = vec![
            BitVec::from([false, true, true, false]),
            BitVec::from([false, false, true, true]),
            BitVec::from([false, false, false, false]),
        ];
        let original = mat.clone();
        gaussian_elimination(&mut mat);
        assert_eq!(mat, original);
    }

    #[test]
    fn test_gaussian_elimination_4() {
        let mut mat = vec![
            BitVec::from([false, false]),
            BitVec::from([false, true]),
            BitVec::from([false, false]),
            BitVec::from([true, true]),
        ];
        gaussian_elimination(&mut mat);
        // Zero rows should be at the bottom.
        assert_eq!(
            mat,
            vec![
                BitVec::from([true, true]),
                BitVec::from([false, true]),
                BitVec::from([false, false]),
                BitVec::from([false, false]),
            ]
        );
    }
}
