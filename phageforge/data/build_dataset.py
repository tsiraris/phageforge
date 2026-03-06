from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Optional
import pandas as pd
from tqdm import tqdm
from Bio import Entrez, SeqIO
from Bio.SeqFeature import SeqFeature
import re
import time
import socket
from urllib.error import HTTPError, URLError
from io import StringIO

# ESKAPEE host genera
ESKAPEE = {
    "Enterococcus",
    "Staphylococcus",
    "Klebsiella",
    "Acinetobacter",
    "Pseudomonas",
    "Enterobacter",
    "Escherichia",
}

# RBP keywords
RBP_KEYWORDS = [
    "fiber",
    "fibre",
    "adhesin",
    "tail fiber",
    "tail fibre",
    "tail spike",
    "tailspike",
    "receptor binding",
    "receptor-binding",
    "rbp",
    "baseplate",
    "adsorption",
    "host recognition",
    "host-recognition"
]

@dataclass(frozen=True) # Immutable class
class RBPRecord:
    virus_accession: str
    host_genus: str
    protein_id: str
    product: str
    aa_sequence: str


def _host_genus_from_name(host_name: str) -> Optional[str]:
    """ 
    Extract host genus from host name. 
    Takes as input the "host name" field from the Virus-Host DB, keeps the first word and returns it if it is in ESKAPEE. 
    """
    # Virus-Host database's "host name" can contain multiple entries separated by commas.
    if not isinstance(host_name, str) or not host_name.strip(): 
        return None # if host_name is None or empty
    # Keep first word of the host name (for example for Echerichia coli, keep "Echerichia")
    genus = host_name.strip().split()[0]
    return genus if genus in ESKAPEE else None  # return None if not in ESKAPEE

# Accept typical nucleotide accessions:
# NC_#####, NZ_#####, MN#####, etc. (we keep it permissive but sane)
ACC_RE = re.compile(r"^[A-Z]{1,3}_?\d+(\.\d+)?$")

def _extract_accessions(refseq_field: str) -> list[str]:
    """ Extracts GenBank accessions from the "refseq id" field. """
    
    # If the "refseq id" field is empty or not a string, return an empty list
    if not isinstance(refseq_field, str) or not refseq_field.strip():  
        return []
    
    # Replace all semicolons with commas and then split on commas
    parts = [p.strip() for p in refseq_field.replace(";", ",").split(",")] 
    out = []
    
    # For each part
    for p in parts:
        if not p:
            continue
        p = p.split()[0]           # remove trailing junk
        p = p.split("|")[-1]       # handle pipe-formatted fields
        p = p.strip()
        if ACC_RE.match(p):
            out.append(p.split(".")[0])  # strip version
    return out

def _is_candidate_rbp(feature: SeqFeature) -> bool:
    """ Check if the sequence feature is a CDS (code DNA sequence) and if it contains any of the keywords ("product", "gene", "note", "function") in RBP_KEYWORDS, and returns True if so. """
    # Check if the feature is a CDS
    if feature.type != "CDS":
        return False
    # Check if the CDS feature contains any of the keywords: "product", "gene", "note", "function"
    quals = feature.qualifiers
    text = " ".join(
        str(x).lower()
        for k in ("product", "gene", "note", "function") # Try these keys
        for x in quals.get(k, [])
    )
    return any(kw in text for kw in RBP_KEYWORDS) # Return True if any keyword is found

def _get_protein_id(feature: SeqFeature, fallback_prefix: str) -> str:
    """ Returns the protein ID from the feature sequence, or a fallback if not present. """
    quals = feature.qualifiers
    # Check if the feature has a "protein_id", "locus_tag", or "gene" qualifier
    for key in ("protein_id", "locus_tag", "gene"):
        if key in quals and len(quals[key]) > 0: # If the key is present and has a value
            return str(quals[key][0])            # Return the value
    return f"{fallback_prefix}_cds"              # Return a fallback

def _get_product(feature: SeqFeature) -> str:
    """ Returns the product (protein type/name) from the feature sequence, or "unknown" if not present. """
    quals = feature.qualifiers
    # Check if the feature has a "product" or "note" qualifier and return the first value
    if "product" in quals and len(quals["product"]) > 0:
        return str(quals["product"][0])
    if "note" in quals and len(quals["note"]) > 0:
        return str(quals["note"][0])
    return "unknown" # else return "unknown"

def _get_translation(feature: SeqFeature) -> Optional[str]:
    """ Returns the translation (amino acid sequence) from the feature sequence, or None if not present. """
    quals = feature.qualifiers
    # Check if the feature has a "translation" qualifier and return the first value
    if "translation" in quals and len(quals["translation"]) > 0:
        return str(quals["translation"][0]).replace(" ", "").replace("\n", "")
    return None

def batched(items: list[str], batch_size: int) -> Iterator[list[str]]:
    """ Yields batches of items. """
    for i in range(0, len(items), batch_size):
        yield items[i : i + batch_size] # Returns a generator of lists of length batch_size

def _efetch_genbank(ids_csv: str):
    """ Fetches GenBank records for a comma-separated list of accessions. """
    # db = database, id = comma-separated list of accessions, rettype = record type, retmode = return mode
    return Entrez.efetch(db="nuccore", id=ids_csv, rettype="gb", retmode="text")

def fetch_genbank_records( accessions: list[str], email: str, api_key: Optional[str] = None, cache_dir: Optional[Path] = None, sleep_s: float = 0.34, max_retries: int = 4) -> Iterator[tuple[str, object]]:
    """
    Fetches GenBank records for a list of accessions. 
    
    """
    # Configure Entrez 
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    # Set timeout
    socket.setdefaulttimeout(30)

    # Create cache directory if it doesn't exist
    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)

    # Create an empty set to keep track of seen accessions and remove duplicates
    seen = set() 
    accessions = [a for a in accessions if not (a in seen or seen.add(a))] 

    # Fetch GenBank records
    for acc in accessions:
        # Get the first part of the accession
        acc = acc.split(".")[0]

        cache_path = None
        
        # Check cache: if not empty, try to read from cache the record for the given accessions
        if cache_dir is not None:
            cache_path = cache_dir / f"{acc}.gb"
            if cache_path.exists() and cache_path.stat().st_size > 0:
                try:
                    with open(cache_path, "r", encoding="utf-8", errors="ignore") as f:
                        recs = list(SeqIO.parse(f, "genbank"))
                    if recs:
                        yield acc, recs[0]
                        continue
                except Exception:
                    pass
        
        # Sleep used for rate limiting
        time.sleep(sleep_s)

        last_err = None
        
        # Try to fetch the record
        for attempt in range(max_retries):  # Retry up to max_retries times
            try:                            # Try to fetch the record
                handle = Entrez.efetch(
                    db="nuccore",
                    id=acc,
                    rettype="gb",
                    retmode="text"
                )
                text = handle.read()    # Read the record text
                handle.close()

                # Write the record text to cache
                if cache_path is not None:
                    with open(cache_path, "w", encoding="utf-8") as f:
                        f.write(text)

                # 
                recs = list(SeqIO.parse(StringIO(text), "genbank"))
                if recs:
                    yield acc, recs[0]
                break
            
            # Handle errors: if the error is a 400 (bad request), skip the record
            except HTTPError as e:
                if getattr(e, "code", None) == 400:
                    print(f"[WARN] skipping accession (HTTP 400): {acc}")
                    break
                last_err = e # Save the error as a variable  

            # Handle other errors
            except (URLError, TimeoutError, ConnectionError, OSError) as e: 
                last_err = e

            time.sleep(1.5 * (attempt + 1)) # Sleep for 1.5 seconds

        else:
            print(f"[WARN] giving up after retries: {acc} | last_err={last_err}")


def build_rbp_dataset(virushost_tsv: Path, out_csv: Path, cache_dir: Path, email: str, api_key: Optional[str] = None, max_viruses_per_genus: int = 300) -> Path:
    """
    Stage 1 dataset:
      - reads Virus-Host DB TSV (Which phage infects which bacterium)
      - filters to ESKAPEE host (bacteria) genera 
      - fetches GenBank for virus (bacteriophage) accessions (What proteins does that phage encode)
      - extracts CDS with RBP-like keywords
      - writes rbp_dataset.csv
    """
    # Create output directory
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    #---- Virus-Host TSV Filtering: Only for ESKAPEE host genera, keep only the refseq id of the viruses and the corresponding bacteria genus names first word----
    
    # Read Virus-Host (phage-host mapping) database as pandas dataframe
    df = pd.read_csv(virushost_tsv, sep="\t", dtype=str)

    # Only for ESKAPEE host genera, add host_genus column and drop rows where host_genus is None
    df["host_genus"] = df["host name"].apply(_host_genus_from_name) if "host name" in df.columns else None  
    df = df.dropna(subset=["host_genus"])   
    
    # Keep only viruses with refseq id not empty rows, add accessions column and drop rows where accessions is empty
    if "refseq id" not in df.columns:
        raise ValueError("Expected column 'refseq id' not found in Virus-Host DB TSV.")
    df["accessions"] = df["refseq id"].apply(_extract_accessions) 
    df = df[df["accessions"].map(len) > 0]  

    # Add virus_accession column, strip the version from each row, and drop duplicates.
    df["virus_accession"] = df["accessions"].apply(lambda xs: xs[0].split(".")[0]) 
    df = df.drop_duplicates(subset=["virus_accession", "host_genus"])

    # Keep different max_viruses_per_genus accessions for each host_genus
    CAPS = {
    "Klebsiella": 300,
    "Acinetobacter": 300,
    "Pseudomonas": 200,
    "Escherichia": 200,
    "Enterobacter": 150,
    "Staphylococcus": 150,
    "Enterococcus": 150,
    }

    def cap_group(g):
        """ Cap the number of accessions for each host_genus """
        cap = CAPS.get(g.name, max_viruses_per_genus)
        return g.head(cap) # Return the first cap rows

    # Group by host_genus, keep only up to first cap max_viruses_per_genus accessions for each host_genus and reset index
    df = df.groupby("host_genus", group_keys=False).apply(cap_group).reset_index(drop=True)

    # Create a host map: a dictionary including only the virus genome accession (unversioned refseq id) to host_genus
    accessions = df["virus_accession"].tolist()  # Convert the virus_accession column to a list
    host_map = dict(zip(df["virus_accession"], df["host_genus"]))   

    #---- GenBank Fetching: Use the virus accessions to fetch GenBank records, filter CDS features with RBP-like keywords eventually keeping "product" = protein type, "translation"= amico acid sequence and "protein_id" for ESKAPEE genuses ----
    
    rows: list[dict] = [] # Create an empty list to store the rows

    # Fetch GenBank records
    it = fetch_genbank_records(accessions, email=email, api_key=api_key, cache_dir=cache_dir / "genbank")
    
    # Iterate over the GenBank records
    for acc, rec in tqdm(it, total=len(accessions)):    # for each tuple of (accession, record)
        host_genus = host_map.get(acc)  # Get the host_genus for the current accessions
        
        # Skip if host_genus is None
        if host_genus is None:
            continue

        fallback_prefix = acc
        for idx, feat in enumerate(rec.features): # for each feature
            # Skip if the feature is not a CDS
            if not _is_candidate_rbp(feat): 
                continue
            # Skip if the feature does not have a translation
            aa = _get_translation(feat)
            if not aa or len(aa) < 60:
                continue
            
            # Get the protein ID and product
            pid = _get_protein_id(feat, fallback_prefix=f"{fallback_prefix}_{idx}")
            product = _get_product(feat)

            # Append the row
            rows.append(
                {
                    "virus_accession": acc,
                    "host_genus": host_genus,
                    "protein_id": pid,
                    "product": product,
                    "aa_sequence": aa,
                }
            )

    # Write the dataset as a CSV
    out_df = pd.DataFrame(rows).drop_duplicates(subset=["aa_sequence"])
    out_df.to_csv(out_csv, index=False)
    return out_csv