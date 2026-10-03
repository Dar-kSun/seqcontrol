"""Download-once file cache under data/raw/."""

from __future__ import annotations

import hashlib
import urllib.request
from collections.abc import Sequence
from pathlib import Path

USER_AGENT = "seqcontrol/0.1 (research)"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fetch(urls: str | Sequence[str], dest: Path) -> Path:
    """Download the first URL that works to `dest`, unless `dest` already exists."""
    if dest.exists():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    errors = []
    for url in [urls] if isinstance(urls, str) else urls:
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as f:
                while block := r.read(1 << 20):
                    f.write(block)
            tmp.replace(dest)
            return dest
        except OSError as e:  # URLError and HTTPError are OSErrors
            tmp.unlink(missing_ok=True)
            errors.append(f"{url}: {e}")
    raise RuntimeError("all downloads failed:\n" + "\n".join(errors))
