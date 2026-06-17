from __future__ import annotations
import logging
from pathlib import Path
from typing import Iterable, List, Any
import re

_INT_RE = re.compile(r"-?\d+")

def setup_logging(verbosity: int = 0) -> None:
    level = logging.WARNING
    if verbosity == 1: level = logging.INFO
    if verbosity >= 2: level = logging.DEBUG
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )

def abs_path(p: str | Path) -> str:
    return str(Path(p).expanduser().resolve())

def parse_index_list(x: Any) -> List[int]:
    """Parse list-like parquet values into a list[int]."""
    try:
        import pandas as pd
        if x is None or (isinstance(x, float) and pd.isna(x)):
            return []
    except Exception:
        if x is None:
            return []

    if isinstance(x, (list, tuple)):
        out=[]
        for v in x:
            try: out.append(int(v))
            except Exception: pass
        return out
    try:
        if hasattr(x, "tolist"):
            return parse_index_list(x.tolist())
    except Exception:
        pass

    s=str(x).strip()
    if not s: return []
    return [int(n) for n in _INT_RE.findall(s)]

def contig_to_rfd3(contig_string: str) -> str:
    return ",".join([seg for seg in str(contig_string).split("/") if seg != ""])
