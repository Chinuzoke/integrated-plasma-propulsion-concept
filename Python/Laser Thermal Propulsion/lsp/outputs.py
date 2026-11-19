"""
outputs.py -- output-directory resolution and verified file writing.

Why this exists: output directories used to be given relative to the *current working directory* and the applications printed "saved ..." with
that relative path. Run from another directory the files appeared somewhere unexpected, and nothing confirmed that a file existed after saving.
Now:  * default output directories are anchored to the PROJECT ROOT (independent of the current directory);
      * a user-supplied --outdir is resolved to an ABSOLUTE path;
      * every file is checked (exists, non-zero size) after it is written, and the absolute path is reported;
      * a failure raises OutputError instead of printing a success message.
"""
from __future__ import annotations

import os

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class OutputError(RuntimeError):
    pass


def resolve_outdir(outdir: str | None, default_subdir: str) -> str:
    """Absolute output directory (created). None -> <project root>/outputs/<default_subdir>."""
    path = os.path.abspath(outdir) if outdir else os.path.join(PROJECT_ROOT, "outputs", default_subdir)
    os.makedirs(path, exist_ok=True)
    if not os.path.isdir(path) or not os.access(path, os.W_OK):
        raise OutputError(f"output directory is not writable: {path}")
    return path


def verify_file(path: str) -> dict:
    ap = os.path.abspath(path)
    exists = os.path.isfile(ap)
    size = os.path.getsize(ap) if exists else 0
    return {"path": ap, "exists": exists, "bytes": size, "ok": exists and size > 0}


def save_figure(fig, path: str, dpi: int = 140) -> str:
    """Save and VERIFY a matplotlib figure; returns the absolute path or raises OutputError."""
    import matplotlib.pyplot as plt
    ap = os.path.abspath(path)
    try:
        fig.tight_layout()
        fig.savefig(ap, dpi=dpi)
    except Exception as e:                       # never swallow: report and re-raise as an output failure
        raise OutputError(f"could not save figure {ap}: {e!r}") from e
    finally:
        plt.close(fig)
    v = verify_file(ap)
    if not v["ok"]:
        raise OutputError(f"figure was not written (missing or empty): {ap}")
    return ap


def write_text(path: str, text: str) -> str:
    ap = os.path.abspath(path)
    with open(ap, "w", encoding="utf-8") as f:
        f.write(text)
    if not verify_file(ap)["ok"]:
        raise OutputError(f"text file was not written: {ap}")
    return ap


def verification_table(paths) -> list:
    return [verify_file(p) for p in paths]
