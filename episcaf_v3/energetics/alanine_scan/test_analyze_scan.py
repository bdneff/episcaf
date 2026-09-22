"""Synthetic format tests, not a claim of validation against a real Gemini run."""
import csv
from pathlib import Path
import tempfile
import unittest
from analyze_scan import parse


class PairedFrames(unittest.TestCase):
    def fixture(self, directory, mutant_frames=(1,21)):
        path = Path(directory) / 'energies.csv'
        with path.open('w', newline='') as f:
            out = csv.writer(f)
            for label, frames, offset in [('GENERALIZED BORN:', (1,21), 0),
                    ('LYS Mutant \nGENERALIZED BORN:\n\n', mutant_frames, 2)]:
                out.writerow([label])
                out.writerow(['Complex Energy Terms'])
                out.writerow(['Frame #','TOTAL'])
                out.writerow([1, -999])  # must not mistake complex energy for binding
                out.writerow([])
                out.writerow(['Delta Energy Terms'])
                out.writerow(['Frame #','VDWAALS','TOTAL'])
                for frame in frames:
                    out.writerow([frame,-3,-10+offset])
                out.writerow([])
        return path

    def test_select_binding_and_pair_sign(self):
        with tempfile.TemporaryDirectory() as tmp:
            data = parse(self.fixture(tmp))
            self.assertEqual(data['mutant'][1]['TOTAL']-data['wt'][1]['TOTAL'],2)

    def test_missing_mutant_frame_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                parse(self.fixture(tmp,(1,)))

    def test_duplicate_frame_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                parse(self.fixture(tmp,(1,1)))


if __name__ == '__main__':
    unittest.main()
