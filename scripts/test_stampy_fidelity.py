import tempfile
import unittest
from pathlib import Path
from check_stampy_fidelity import compare


class FidelityTest(unittest.TestCase):
    def compare_text(self, actual):
        base = '@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:1000\n'
        record = 'r\t99\tchr1\t21\t99\t4M\t=\t40\t23\tACGT\tIIII\tNM:i:0\tSM:i:80\n'
        with tempfile.TemporaryDirectory() as directory:
            a, b = Path(directory) / 'a.sam', Path(directory) / 'b.sam'
            a.write_text(base + record)
            b.write_text(actual(base + record))
            return compare(a, b)

    def test_equivalent_notation(self):
        self.assertTrue(self.compare_text(lambda s: s.replace('\t=\t', '\tchr1\t')
                        .replace('NM:i:0\tSM:i:80', 'SM:i:80\tNM:i:0'))['passed'])

    def test_one_base_difference_fails(self):
        self.assertFalse(self.compare_text(lambda s: s.replace('\t21\t', '\t22\t'))['passed'])

    def test_missing_tag_fails(self):
        self.assertFalse(self.compare_text(lambda s: s.replace('\tNM:i:0', ''))['passed'])

    def test_duplicate_record_fails(self):
        self.assertFalse(self.compare_text(lambda s: s + s.splitlines(True)[-1])['passed'])

    def test_mapq_zero_is_not_ignored(self):
        self.assertFalse(self.compare_text(lambda s: s.replace('\t99\t4M', '\t0\t4M'))['passed'])

    def test_header_difference_fails(self):
        self.assertFalse(self.compare_text(lambda s: s.replace('VN:1.0', 'VN:1.6'))['passed'])


if __name__ == '__main__':
    unittest.main()
