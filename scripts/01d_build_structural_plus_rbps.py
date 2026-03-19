from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import re


def normalize_product_string(product: str) -> str:
    """ Normalize a product annotation string to simplify regex matching. """
    if pd.isna(product):
        return ""
    return re.sub(r"\s+", " ", str(product).strip().lower())                             # Lowercase, trim, and collapse repeated whitespace


def compile_patterns(pattern_strings: list[str]) -> list[re.Pattern]:
    """ Compile a list of regex pattern strings for repeated matching. """
    return [re.compile(p, flags=re.IGNORECASE) for p in pattern_strings]                 # Ignore case


def matches_any_pattern(text: str, patterns: list[re.Pattern]) -> bool:
    """ Return True if the text matches any compiled regex pattern. """
    return any(p.search(text) is not None for p in patterns)


def main():
    # Create an ArgumentParser object to run the same script with different parameters without editing code.
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_csv", type=str, default="data/processed/rbp_dataset_eskapee_stage1.csv")
    ap.add_argument("--output_csv", type=str, default="data/processed/rbp_dataset_eskapee_structural_plus.csv")
    ap.add_argument("--review_csv", type=str, default="data/processed/rbp_dataset_eskapee_structural_plus_review.csv")
    args = ap.parse_args()

    # Load data
    input_csv = Path(args.input_csv)
    output_csv = Path(args.output_csv)
    review_csv = Path(args.review_csv)

    # Raise error if input CSV does not exist   
    if not input_csv.exists():
        raise FileNotFoundError(f"Missing file: {input_csv}")

    # Load the stage-1 dataset
    df = pd.read_csv(input_csv)                                                          # Read the broader stage-1 dataset so structural_plus can recover useful adsorption-related proteins that the strict structural filter may miss
    
    # Raise error if input CSV does not contain 'product' column
    if "product" not in df.columns:
        raise ValueError("Expected a 'product' column in the input CSV.")

    # Clean product strings
    df["product_norm"] = df["product"].apply(normalize_product_string)                   # Normalize product strings once so all downstream regex checks use the same cleaned annotation text

    # High-confidence keep patterns inherited from the structural-style host-recognition / adsorption set
    keep_patterns = compile_patterns([
        r"\btail fiber protein\b",
        r"\btail fibre protein\b",
        r"\btail spike protein\b",
        r"\btailspike protein\b",
        r"\breceptor[- ]?binding protein\b",
        r"\breceptor[- ]?binding\b",
        r"\bhost[- ]?range and adsorption protein\b",
        r"\bhost[- ]?recognition\b",
        r"\badsorption protein\b",
        r"\badsorption associated\b",
        r"\battachment protein\b",
        r"\btail fiber receptor[- ]?binding protein\b",
        r"\bputative tail fiber protein\b",
        r"\bputative adsorption protein\b",
        r"\bputative adsorption associated tail protein\b",
        r"\bnon[- ]contractile tail fiber protein\b",
        r"\bshort tail fiber protein\b",
        r"\blong[- ]?tail fiber protein\b",
        r"\bl[- ]shaped tail fiber protein\b",
        r"\blytic tail fiber protein\b",
        r"\bhead fiber protein\b",
    ])                                                                                   

    # Expanded keep patterns that recover adsorption-apparatus and host-interaction neighborhood annotations for a broader but still curated dataset
    expanded_keep_patterns = compile_patterns([
        r"\bdepolymerase\b",
        r"\bcapsular polysaccharide depolymerase\b",
        r"\btail tip\b",
        r"\bdistal tail\b",
        r"\breceptor recognition\b",
        r"\bhost specificity\b",
        r"\badsorption apparatus\b",
        r"\battachment catalyst\b",
        r"\btail fiber protein attachment catalyst\b",
        r"\btail fiber adhesin\b",
        r"\bcentral tail fiber\b",
        r"\bproximal tail fiber\b",
        r"\blong tail fiber proximal\b",
        r"\blong tail fiber distal\b",
        r"\bhinge connector of long tail fiber\b",
        r"\bcollar tail protein for l[- ]shaped tail fiber attachment\b",
        r"\blong tail fiber protein proximal connector\b",
        r"\blong tail fiber protein distal subunit\b",
        r"\btail collar fiber protein\b",
        r"\bbaseplate protein\b",
        r"\bputative baseplate protein\b",
        r"\bbaseplate spike\b",
        r"\bbaseplate spike protein\b",
        r"\bbaseplate j[- ]like protein\b",
        r"\bbaseplate j family protein\b",
        r"\bbppu family baseplate upper protein\b",
        r"\bbaseplate upper protein\b",
        r"\bbaseplate central spike complex protein\b",
        r"\bbaseplate component\b",
    ])                                              
    
    # Borderline patterns that are useful to inspect but are broad enough that they should stay in the review file by default                                
    review_patterns = compile_patterns([
        r"\btail fiber\b",
        r"\btail fibers protein\b",
        r"\bshort tail fiber\b",
        r"\blong tail fiber\b",
        r"\bbaseplate\b",
        r"\bbaseplate subunit\b",
        r"\bbaseplate tail tube cap\b",
        r"\bbaseplate tail tube initiator\b",
        r"\bbaseplate tail[- ]tube junction protein\b",
        r"\bputative tail fiber\b",
        r"\btail protein\b",
    ])                                                                                   
    
    # Strong exclusion patterns for proteins unlikely to be the primary host-recognition determinant or too noisy for this curated adsorption-enriched panel
    exclude_patterns = compile_patterns([
        r"\bchaperone\b",
        r"\bassembly\b",
        r"\bhead[- ]tail connector\b",
        r"\bwedge\b",
        r"\bhub\b",
        r"\blysozyme\b",
        r"\bmuramidase\b",
        r"\bpeptidase\b",
        r"\bpeptidoglycan hydrolase\b",
        r"\bendolysin\b",
        r"\bhypothetical protein\b",
        r"\bportal protein\b",
        r"\bmajor tail protein\b",
        r"\btail tubular protein\b",
        r"\btail tube protein\b",
        r"\btape measure protein\b",
        r"\bterminase\b",
        r"\bcapsid\b",
        r"\bneck protein\b",
        r"\bRNA ligase\b",
    ])                                                                                   

    # Create masks for each dataset
    keep_mask = df["product_norm"].apply(lambda s: matches_any_pattern(s, keep_patterns))                                   # High-confidence structural-style keep set
    expanded_keep_mask = df["product_norm"].apply(lambda s: matches_any_pattern(s, expanded_keep_patterns))                 # Expanded adsorption / host-interaction keep set
    review_mask = df["product_norm"].apply(lambda s: matches_any_pattern(s, review_patterns))                               # Review-only set for borderline generic annotations
    exclude_mask = df["product_norm"].apply(lambda s: matches_any_pattern(s, exclude_patterns))                             # Exclusion set for clearly non-target annotations

    # Structural plus dataset will keep the original keep patterns and the expanded keep patterns and not the review patterns
    structural_plus_mask = (keep_mask | expanded_keep_mask) & (~exclude_mask)                                               
    
    # Review dataset will keep the review patterns and not the structural_plus patterns or the exclusion patterns
    review_output_mask = review_mask & (~structural_plus_mask) & (~exclude_mask)                                            

    # Create the final datasets
    structural_plus_df = df[structural_plus_mask].copy().reset_index(drop=True)                 # Final structural_plus dataset: broader and higher-recall than structural while still curated
    review_df = df[review_output_mask].copy().reset_index(drop=True)                            # Review dataset: borderline proteins worth inspecting manually but not included in the main structural_plus set

    # Drop the helper normalization column before saving the final datasets
    if "product_norm" in structural_plus_df.columns:
        structural_plus_df = structural_plus_df.drop(columns=["product_norm"])                  
    if "product_norm" in review_df.columns:
        review_df = review_df.drop(columns=["product_norm"])                                    

    # Create output directories if they don't exist
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    review_csv.parent.mkdir(parents=True, exist_ok=True)

    # Save the final datasets to CSV
    structural_plus_df.to_csv(output_csv, index=False)                                          
    review_df.to_csv(review_csv, index=False)                                                   

    # Print success messages
    print(f"✅ Saved structural-plus dataset: {output_csv}")
    print(f"✅ Saved review dataset:          {review_csv}")
    print("")
    print("Structural-plus dataset size:")
    print(structural_plus_df.shape)
    print("")
    print("Structural-plus dataset host distribution:")
    if "host_genus" in structural_plus_df.columns:
        print(structural_plus_df["host_genus"].value_counts())
    else:
        print("Missing host_genus column")
    print("")
    print("Top structural-plus product strings:")
    if "product" in structural_plus_df.columns:
        print(structural_plus_df["product"].value_counts().head(30))
    else:
        print("Missing product column")
    print("")
    print("Review dataset size:")
    print(review_df.shape)
    print("")
    print("Top review product strings:")
    if "product" in review_df.columns:
        print(review_df["product"].value_counts().head(30))
    else:
        print("Missing product column")


if __name__ == "__main__":
    main()
