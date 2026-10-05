from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from checking.conventions import reverse_qubit_bits, reverse_pauli_label, symplectic, msb_x_high_to_shared_label


class ConventionTests(unittest.TestCase):
    def test_basis_and_inverse(self):
        for n in (3,4,5):
            for x in range(1<<n):
                self.assertEqual(reverse_qubit_bits(reverse_qubit_bits(x,n),n),x)
            self.assertEqual(reverse_qubit_bits(1,n),1<<(n-1))

    def test_symplectic_and_hermitian_phase(self):
        n=3
        for a in range(64):
            ar=reverse_pauli_label(a,n)
            self.assertEqual(((a&7)&(a>>3)).bit_count(),((ar&7)&(ar>>3)).bit_count())
            for b in range(64):
                self.assertEqual(symplectic(a,b,n),symplectic(ar,reverse_pauli_label(b,n),n))

    def test_wrong_dimension_rejected(self):
        with self.assertRaises(ValueError):
            reverse_pauli_label(64,3)

    def test_c16_pauli_matrices_equal_after_named_conversion(self):
        from checking.c16 import pauli
        for v in range(64):
            shared=msb_x_high_to_shared_label(v,3)
            x,z=shared&7,shared>>3
            p,f=pauli(v)
            self.assertEqual(msb_x_high_to_shared_label(shared,3),v)
            for j in range(8):
                old=reverse_qubit_bits(j,3)
                self.assertEqual(reverse_qubit_bits(p[old],3),j^x)
                self.assertEqual(f[old],((x&z).bit_count()+2*(z&j).bit_count())%4)


if __name__=='__main__':
    unittest.main()
