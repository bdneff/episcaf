"""Unit tests for tile_fasta.tile_sequence (synthetic; no data)."""
import sys; from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from episcaf_pipeline.tile_fasta import tile_sequence

def test_basic_windows_exact_end():
    seq = "ACDEFGHIKL"          # len 10
    out = list(tile_sequence(seq, 4, 2))   # 1,3,5,7 -> ends at 10 exactly, no cterm
    assert [p for p, _ in out] == [1, 3, 5, 7]
    assert out[0] == (1, "ACDE") and out[-1] == (7, "GHIK".replace("GHIK","GHIK")) or out[-1][1] == seq[6:10]

def test_cterm_added_when_grid_misses_end():
    seq = "ACDEFGHIKLMNO"       # len 13
    out = list(tile_sequence(seq, 4, 4))   # starts 0,4,8 (cover to 12) + cterm at 13-4=9 (start10)
    pos = [p for p, _ in out]
    assert pos == [1, 5, 9, 10], pos
    assert out[-1][1] == seq[-4:]          # final tile reaches the C-terminus

def test_cterm_disabled():
    seq = "ACDEFGHIKLMNO"       # len 13
    out = list(tile_sequence(seq, 4, 4, include_cterm=False))
    assert [p for p, _ in out] == [1, 5, 9]

def test_too_short():
    assert list(tile_sequence("ACD", 4, 2)) == []

if __name__ == "__main__":
    for fn in [test_basic_windows_exact_end, test_cterm_added_when_grid_misses_end,
               test_cterm_disabled, test_too_short]:
        fn(); print("ok", fn.__name__)
    print("ALL TILING TESTS PASS")
