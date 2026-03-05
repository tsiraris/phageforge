from __future__ import annotations

from pathlib import Path
import requests

# tsv files
VIRUS_HOST_DAILY_URL = "https://www.genome.jp/ftp/db/virushostdb/virushostdb.daily.tsv"


def download_file(url: str, out_path: Path, chunk_size: int = 1024 * 1024) -> Path:
    """ Download a file from the given url to the given path."""
    # Create the directory if it doesn't exist
    out_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Check if the file already exists
    if out_path.exists() and out_path.stat().st_size > 0:
        return out_path

    # Download the file to the given path 
    with requests.get(url, stream=True, timeout=60) as r:
        # Check if the request was successful
        r.raise_for_status()
        # Write the content of the request to the file
        with open(out_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
    return out_path # Return the path to the file that was downloaded
