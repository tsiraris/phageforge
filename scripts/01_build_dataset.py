from __future__ import annotations 
import os
from pathlib import Path
from phageforge.data.download import download_file, VIRUS_HOST_DAILY_URL
from phageforge.data.build_dataset import build_rbp_dataset


def main():
    """ What this main does in detail is:
    1. Download Virus-Host DB TSV from NCBI and filter to ESKAPEE host genera
    2. Fetch GenBank for bacteriophage accessions
    3. Extract CDS with RBP-like keywords
    4. Write rbp_dataset.csv
    """
    # Paths
    root = Path(__file__).resolve().parents[1] # root of the repo
    raw = root / "data" / "raw"
    processed = root / "data" / "processed"
    cache = root / "data" / "cache"
    
    virushost_path = raw / "virushostdb.daily.tsv"
    out_csv = processed / "rbp_dataset_eskapee_stage1.csv"

    # Download Virus-Host DB TSV
    download_file(VIRUS_HOST_DAILY_URL, virushost_path)

    # Build dataset
    email = os.environ.get("NCBI_EMAIL", "")
    if not email:
        raise RuntimeError("Set NCBI_EMAIL env var (required by NCBI Entrez usage policies).")

    api_key = os.environ.get("NCBI_API_KEY")  # optional but recommended

    build_rbp_dataset(
        virushost_tsv=virushost_path,
        out_csv=out_csv,
        cache_dir=cache,
        email=email,
        api_key=api_key,
        max_viruses_per_genus=80,  
    )

    print(f"✅ Wrote dataset: {out_csv}")


if __name__ == "__main__":
    main()