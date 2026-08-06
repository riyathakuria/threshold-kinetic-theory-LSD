#!/usr/bin/env python3
"""
Phase 1 - Step 1: Acquire HuMicA v2.0.0 (.h5ad) from Zenodo record 18458280.

- Queries the Zenodo API for the record's file list, sizes, and MD5 checksums.
- Downloads HuMicA_v2.0.0.h5ad into data/raw/ ONLY if absent or integrity fails.
- Streams to disk while computing SHA256 and MD5 on the fly (no full-file re-read).
- Verifies the downloaded MD5 against the Zenodo-published checksum.
- Writes/updates data/download_manifest.csv.

Reproducible: safe to re-run. An already-verified file is reused, not re-downloaded.
"""
import os, sys, json, hashlib, csv, datetime, time, re
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = os.path.join(ROOT, "data", "raw")
MANIFEST = os.path.join(ROOT, "data", "download_manifest.csv")

RECORD_ID = "18458280"
TARGET_KEY = "HuMicA_v2.0.0.h5ad"
API_URL = f"https://zenodo.org/api/records/{RECORD_ID}"


def get_record():
    with urllib.request.urlopen(API_URL, timeout=60) as resp:
        return json.load(resp)


def file_md5(path, chunk=8 * 1024 * 1024):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(chunk), b""):
            h.update(blk)
    return h.hexdigest()


def resumable_download(url, dest, expected_size, max_retries=40):
    """
    Download url -> dest with POSITION-AUTHORITATIVE HTTP Range resume.

    part = dest + ".part"
    # ensure file exists for r+b seeking
    if not os.path.exists(part):
        open(part, "wb").close()
    attempt = 0
    while True:
        have = os.path.getsize(part)
        if expected_size and have >= expected_size:
            break
        req = urllib.request.Request(url)
        req.add_header("Range", f"bytes={have}-")
        try:
            with urllib.request.urlopen(req, timeout=180) as resp:
                # Determine the true starting offset the server is sending from.
                cr = resp.headers.get("Content-Range")  # e.g. "bytes 2859485756-6803451710/6803451711"
                if resp.status == 200:
                    start = 0            # server ignored Range -> full body from 0
                elif cr:
                    m = re.search(r"bytes\s+(\d+)-", cr)
                    start = int(m.group(1)) if m else have
                else:
                    start = have
                print(f"  resuming from {have/1e9:.2f} GB; server sends from {start/1e9:.2f} GB "
                      f"(attempt {attempt+1})", flush=True)
                last_pct = -1
                with open(part, "r+b") as out:
                    out.seek(start)
                    pos = start
                    while True:
                        buf = resp.read(8 * 1024 * 1024)
                        if not buf:
                            break
                        # never write past expected_size
                        if expected_size and pos + len(buf) > expected_size:
                            buf = buf[: expected_size - pos]
                            if buf:
                                out.write(buf); pos += len(buf)
                            break
                        out.write(buf); pos += len(buf)
                        if expected_size:
                            pct = int(pos * 100 / expected_size)
                            if pct != last_pct and pct % 5 == 0:
                                print(f"  ... {pct}%  ({pos/1e9:.2f} / {expected_size/1e9:.2f} GB)", flush=True)
                                last_pct = pct
                # if we wrote a contiguous prefix shorter than the current file
                # size (shouldn't happen), truncate stray tail
                with open(part, "r+b") as out:
                    out.seek(0, os.SEEK_END)
                    if expected_size and out.tell() > expected_size:
                        out.truncate(expected_size)
        except Exception as e:
            attempt += 1
            if attempt > max_retries:
                raise
            print(f"  connection error ({type(e).__name__}: {e}); retry {attempt}/{max_retries}", flush=True)
            time.sleep(min(30, 2 ** min(attempt, 5)))
            continue
    final_size = os.path.getsize(part)
    if expected_size and final_size != expected_size:
        raise RuntimeError(f"size mismatch after download: {final_size} != {expected_size}")
    sha, md5 = hashlib.sha256(), hashlib.md5()
    with open(part, "rb") as f:
        for blk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            sha.update(blk); md5.update(blk)
    os.replace(part, dest)
    return sha.hexdigest(), md5.hexdigest(), final_size


# Backwards-compatible alias
def stream_download(url, dest, expected_size):
    return resumable_download(url, dest, expected_size)


def main():
    os.makedirs(RAW_DIR, exist_ok=True)
    print(f"Querying Zenodo record {RECORD_ID} ...", flush=True)
    rec = get_record()
    files = {f["key"]: f for f in rec.get("files", [])}
    if TARGET_KEY not in files:
        sys.exit(f"ERROR: {TARGET_KEY} not in record. Available: {list(files)}")

    fmeta = files[TARGET_KEY]
    size = fmeta.get("size", 0)
    zchecksum = fmeta.get("checksum", "")          # e.g. 'md5:0630...'
    zmd5 = zchecksum.split(":", 1)[-1] if ":" in zchecksum else zchecksum
    url = f"https://zenodo.org/records/{RECORD_ID}/files/{TARGET_KEY}?download=1"
    dest = os.path.join(RAW_DIR, TARGET_KEY)

    print(f"Target: {TARGET_KEY}  size={size/1e9:.3f} GB  zenodo_md5={zmd5}", flush=True)

    # Reuse if present and MD5 matches
    if os.path.exists(dest) and os.path.getsize(dest) == size:
        print("File already present with matching size; verifying MD5 ...", flush=True)
        local_md5 = file_md5(dest)
        if local_md5 == zmd5:
            print("MD5 MATCH -> reusing existing file, skipping download.", flush=True)
            sha256 = None
            md5 = local_md5
        else:
            print(f"MD5 MISMATCH (local {local_md5} != zenodo {zmd5}); re-downloading.", flush=True)
            sha256, md5, _ = stream_download(url, dest, size)
    else:
        print("Downloading ...", flush=True)
        sha256, md5, written = stream_download(url, dest, size)
        print(f"Downloaded {written/1e9:.3f} GB", flush=True)

    ok = (md5 == zmd5)
    print(f"MD5 verification: {'PASS' if ok else 'FAIL'} (local={md5}, zenodo={zmd5})", flush=True)

    # If we reused the file we never computed sha256; compute it now for the manifest
    if sha256 is None:
        print("Computing SHA256 for manifest ...", flush=True)
        h = hashlib.sha256()
        with open(dest, "rb") as f:
            for blk in iter(lambda: f.read(8 * 1024 * 1024), b""):
                h.update(blk)
        sha256 = h.hexdigest()

    row = {
        "Dataset": "HuMicA v2.0.0 (.h5ad)",
        "URL": url,
        "Local Path": dest,
        "Size": os.path.getsize(dest),
        "Download Date": datetime.date.today().isoformat(),
        "SHA256": sha256,
        "Zenodo_MD5": zmd5,
        "MD5_verified": str(ok),
    }
    write_header = not os.path.exists(MANIFEST)
    # Rebuild manifest row-by-key so re-runs update rather than duplicate
    rows = []
    if os.path.exists(MANIFEST):
        with open(MANIFEST) as f:
            rows = [r for r in csv.DictReader(f) if r.get("Dataset") != row["Dataset"]]
    rows.append(row)
    with open(MANIFEST, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(row.keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"Manifest updated: {MANIFEST}", flush=True)
    print("DONE.", flush=True)


if __name__ == "__main__":
    main()
