#!/usr/bin/env python3
"""Write a deterministic SHA-256 manifest for an analysis-product directory."""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output-name", default="hashes.sha256")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    directory = args.directory.resolve()
    output = directory / args.output_name
    staged = output.with_name(output.name + ".inprogress")
    if not directory.is_dir():
        raise NotADirectoryError(directory)
    if staged.exists():
        raise FileExistsError(staged)
    if output.exists() and not args.overwrite:
        raise FileExistsError(f"refusing to overwrite {output}")

    files = sorted(
        path
        for path in directory.rglob("*")
        if path.is_file()
        and path not in {output, staged}
        and not path.name.endswith(".inprogress")
    )
    if not files:
        raise RuntimeError(f"no files to hash under {directory}")
    lines = [f"{sha256(path)}  {path.relative_to(directory).as_posix()}" for path in files]
    staged.write_text("\n".join(lines) + "\n", encoding="utf-8")
    os.replace(staged, output)
    print(f"[done] {len(files)} files -> {output}")


if __name__ == "__main__":
    main()
