import csv
import hashlib
import sys
import time
import urllib.request
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = PROJECT_ROOT / "data" / "raw" / "npc"
MANIFEST = RAW_DIR / "download_manifest.csv"

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
FILES = [
    # (GSE, filename, expected_size_bytes)
    ("GSE221nnn/GSE221609", "GSE221609_barcodes.tsv.gz",     378555),
    ("GSE221nnn/GSE221609", "GSE221609_features.tsv.gz",     290896),
    ("GSE221nnn/GSE221609", "GSE221609_matrix.mtx.gz",   730628877),
    ("GSE152nnn/GSE152158", "GSE152158_NormCountData.csv.gz", 4816848),
]

# Genotype key for GSE221609 (barcode suffix -> label), recorded for step 07.
GSE221609_SAMPLE_KEY = {
    "-1": ("GSM6890285", "WT",      "WT forebrain 1"),
    "-2": ("GSM6890286", "WT",      "WT forebrain 2"),
    "-3": ("GSM6890287", "WT",      "WT forebrain 3"),
    "-4": ("GSM6890288", "Npc1_KO", "Npc1-/- forebrain 1"),
    "-5": ("GSM6890289", "Npc1_KO", "Npc1-/- forebrain 2"),
    "-6": ("GSM6890290", "Npc1_KO", "Npc1-/- forebrain 3"),
}


def _log(msg: str) -> None:
    print(f"[07a] {msg}", flush=True)


def md5_of(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def download(url: str, dest: Path, expected_size: int | None) -> None:
    if dest.exists() and expected_size and dest.stat().st_size == expected_size:
        _log(f"skip (present, size ok): {dest.name}")
        return
    tmp = dest.with_suffix(dest.suffix + ".part")
    _log(f"download: {dest.name}  <- {url}")
    t0 = time.time()
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=120) as resp, open(tmp, "wb") as out:
        total = 0
        while True:
            buf = resp.read(1 << 20)
            if not buf:
                break
            out.write(buf)
            total += len(buf)
    tmp.rename(dest)
    _log(f"  done {total:,} bytes in {time.time()-t0:.1f}s")


def main() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for gse, fname, size in FILES:
        url = f"{BASE}/{gse}/suppl/{fname}"
        dest = RAW_DIR / fname
        try:
            download(url, dest, size)
        except Exception as e:  # noqa: BLE001
            _log(f"ERROR downloading {fname}: {e}")
            sys.exit(1)
        digest = md5_of(dest)
        rows.append({
            "gse": gse.split("/")[-1],
            "filename": fname,
            "bytes": dest.stat().st_size,
            "md5": digest,
            "url": url,
        })
        _log(f"  md5={digest}")

    with open(MANIFEST, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["gse", "filename", "bytes", "md5", "url"])
        w.writeheader()
        w.writerows(rows)
    _log(f"manifest -> {MANIFEST}")

    # Persist the GSE221609 genotype key for step 07.
    key_path = RAW_DIR / "GSE221609_sample_key.csv"
    with open(key_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["barcode_suffix", "gsm", "genotype", "title"])
        for suf, (gsm, geno, title) in GSE221609_SAMPLE_KEY.items():
            w.writerow([suf, gsm, geno, title])
    _log(f"sample key -> {key_path}")
    _log("acquisition complete")


if __name__ == "__main__":
    main()
