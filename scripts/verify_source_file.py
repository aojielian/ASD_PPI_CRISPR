#!/usr/bin/env python3
"""Check one locally acquired source file without downloading it."""
import argparse, hashlib
from pathlib import Path

def sha256(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--path",required=True,type=Path)
    ap.add_argument("--expected-sha256",required=True)
    ap.add_argument("--expected-bytes",required=True,type=int)
    a=ap.parse_args(); p=a.path.resolve(strict=True)
    size=p.stat().st_size; digest=sha256(p)
    if size!=a.expected_bytes or digest.lower()!=a.expected_sha256.lower():
        raise SystemExit(f"FAIL size={size} sha256={digest}")
    print(f"PASS bytes={size} sha256={digest}")
if __name__=="__main__":main()
