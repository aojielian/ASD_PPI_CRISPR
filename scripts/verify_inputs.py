#!/usr/bin/env python3
"""Verify files from the Zenodo v1.0.0 archive or a legacy manifest-based bundle.

Only data-file integrity is checked; this does not rerun statistical models.
"""
import argparse
import csv
import hashlib
import re
from pathlib import Path, PurePosixPath

HEX64 = re.compile(r"^[0-9a-fA-F]{64}$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rows_from_zenodo_sha256s(path: Path):
    for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        fields = line.split("  ", 1)
        if len(fields) != 2 or not HEX64.fullmatch(fields[0]):
            raise ValueError(f"Invalid SHA256SUMS.txt line {lineno}")
        yield fields[1], fields[0].lower(), None


def rows_from_legacy_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"relative_path", "sha256", "bytes"} <= set(reader.fieldnames or []):
            raise ValueError("Legacy manifest must have relative_path, sha256, bytes")
        for row in reader:
            if not HEX64.fullmatch(row["sha256"]):
                raise ValueError(f"Invalid SHA-256 in legacy manifest: {row['relative_path']!r}")
            yield row["relative_path"], row["sha256"].lower(), int(row["bytes"])


def verify(root: Path) -> int:
    root = root.resolve(strict=True)
    zenodo = root / "SHA256SUMS.txt"
    legacy = root / "metadata" / "MANIFEST_SHA256.tsv"
    if zenodo.is_file():
        entries = list(rows_from_zenodo_sha256s(zenodo))
        excluded = {"SHA256SUMS.txt"}
    elif legacy.is_file():
        entries = list(rows_from_legacy_tsv(legacy))
        excluded = {"metadata/MANIFEST_SHA256.tsv"}
    else:
        raise ValueError("Missing SHA256SUMS.txt or metadata/MANIFEST_SHA256.tsv in extracted input root")
    if not entries:
        raise ValueError("Empty checksum manifest")
    seen = set()
    failures = []
    for rel_text, expected_hash, expected_bytes in entries:
        rel = PurePosixPath(rel_text)
        if rel.is_absolute() or any(part in {'.', '..'} for part in rel.parts) or '\\' in rel_text or not rel.parts:
            failures.append(f"{rel_text}: unsafe relative path")
            continue
        rel_str = rel.as_posix()
        if rel_str in seen:
            failures.append(f"{rel_text}: duplicate manifest path")
            continue
        seen.add(rel_str)
        file = root.joinpath(*rel.parts)
        if file.is_symlink() or any(parent.is_symlink() for parent in file.parents if parent != root and root in parent.parents):
            failures.append(f"{rel_text}: symlink not allowed")
            continue
        if not file.is_file() or not file.resolve().is_relative_to(root):
            failures.append(f"{rel_text}: missing or outside root")
            continue
        if expected_bytes is not None and file.stat().st_size != expected_bytes:
            failures.append(f"{rel_text}: byte-count mismatch")
            continue
        if sha256(file) != expected_hash:
            failures.append(f"{rel_text}: SHA-256 mismatch")
    actual = {p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file()}
    extras = actual - seen - excluded
    if extras:
        failures.extend(f"{name}: unlisted file" for name in sorted(extras))
    if failures:
        raise ValueError("Input verification failed: " + "; ".join(failures))
    print(f"PASS: {len(entries)} files verified by SHA-256 under {root}")
    return len(entries)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', required=True, type=Path)
    args = parser.parse_args()
    try:
        verify(args.data_dir)
    except (ValueError, OSError) as exc:
        parser.exit(1, f"FAIL: {exc}\n")

if __name__ == '__main__':
    main()
