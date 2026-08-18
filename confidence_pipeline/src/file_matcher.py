"""Core dataframe filtering and substring matching logic."""

from __future__ import annotations

# Import logging so per-file match statistics can be reported
import logging
# Import Any because configuration objects are passed as generic mappings
from typing import Any

# Import pandas for row filtering and sequence column operations
import pandas as pd

# Import normalization helpers from the local string utilities module
from src.string_normalizer import is_valid_peptide_sequence, normalize_sequence

# Create a module-level logger for matcher messages
logger = logging.getLogger(__name__)


class ProteinSubstringIndex:
    """In-memory canonical protein sequences used for substring exclusion."""

    def __init__(self, proteins: list[str], *, il_equivalent: bool = False) -> None:
        """Store protein sequences and optional I/L-normalized copies."""
        # Save the original uppercase protein sequence tuple
        self.proteins = tuple(proteins)
        # Save whether I/L equivalence should be used during substring checks
        self.il_equivalent = il_equivalent
        # Build an optional second tuple with isoleucine mapped to leucine
        if il_equivalent:
            # Create I/L-normalized protein copies for isobaric substring matching
            self.proteins_il = tuple(protein.replace("I", "L") for protein in self.proteins)
        else:
            # Use None when I/L equivalence is disabled
            self.proteins_il = None


def build_protein_substring_index(
    protein_sequences: list[str],
    *,
    il_equivalent: bool = False,
) -> ProteinSubstringIndex:
    """Build a protein substring index from canonical proteome sequences."""
    # Construct and return the in-memory protein index object
    return ProteinSubstringIndex(protein_sequences, il_equivalent=il_equivalent)


def is_known_fragment(
    peptide: str,
    index: ProteinSubstringIndex,
    *,
    strip_mods: bool = True,
    il_equivalent: bool = True,
) -> bool:
    """Return whether a peptide is an exact substring of any canonical protein."""
    # Normalize the predicted peptide for substring comparison
    normalized_peptide = normalize_sequence(
        peptide,
        strip_mods=strip_mods,
        il_equivalent=il_equivalent,
    )
    # Reject empty or invalid peptide strings without searching the proteome
    if not is_valid_peptide_sequence(normalized_peptide):
        # Return false so invalid sequences are kept for downstream review
        return False

    # Choose the protein list that corresponds to the configured I/L behavior
    protein_list = index.proteins_il if il_equivalent and index.proteins_il else index.proteins
    # Scan every canonical protein for the normalized peptide substring
    for protein in protein_list:
        # Return true immediately when the peptide is found inside a protein
        if normalized_peptide in protein:
            return True

    # Return false when no canonical protein contains the peptide substring
    return False


def build_sequence_set(
    dataframe: pd.DataFrame,
    sequence_column: str,
    *,
    strip_mods: bool = True,
    il_equivalent: bool = False,
) -> set[str]:
    """Build a normalized set of sequence strings from a DataFrame column."""
    # Raise a helpful error when the requested sequence column is missing
    if sequence_column not in dataframe.columns:
        # Stop because sequence set construction cannot proceed without the column
        raise KeyError(
            f"Column '{sequence_column}' not found. Available: {list(dataframe.columns)}"
        )

    # Initialize the set that will store normalized sequence values
    sequences: set[str] = set()
    # Normalize every non-null sequence value in the target column
    for value in dataframe[sequence_column].dropna().astype(str):
        # Normalize the current sequence using the configured options
        normalized = normalize_sequence(
            value,
            strip_mods=strip_mods,
            il_equivalent=il_equivalent,
        )
        # Skip empty strings that may appear after normalization
        if normalized:
            # Add the normalized sequence to the aggregate set
            sequences.add(normalized)

    # Return the completed normalized sequence set
    return sequences


def filter_rows_by_membership(
    dataframe: pd.DataFrame,
    sequence_column: str,
    reference_sequences: set[str],
    *,
    strip_mods: bool = True,
    il_equivalent: bool = False,
) -> pd.DataFrame:
    """Keep rows whose normalized sequence is exactly present in a reference set."""
    # Raise a helpful error when the requested sequence column is missing
    if sequence_column not in dataframe.columns:
        # Stop because exact membership filtering cannot proceed without the column
        raise KeyError(
            f"Column '{sequence_column}' not found. Available: {list(dataframe.columns)}"
        )

    # Normalize each row sequence for exact set membership comparison
    normalized_series = dataframe[sequence_column].astype(str).apply(
        lambda value: normalize_sequence(
            value,
            strip_mods=strip_mods,
            il_equivalent=il_equivalent,
        )
    )
    # Build a boolean mask for rows whose normalized sequence is in the reference set
    match_mask = normalized_series.isin(reference_sequences)
    # Return only the rows that exactly match a reference sequence
    return dataframe.loc[match_mask].copy()


def filter_rows_by_containment(
    dataframe: pd.DataFrame,
    sequence_column: str,
    reference_sequences: list[str],
    *,
    direction: str = "query_in_reference",
    strip_mods: bool = True,
    il_equivalent: bool = False,
) -> pd.DataFrame:
    """Keep rows based on substring containment between query and reference sequences."""
    # Raise a helpful error when the requested sequence column is missing
    if sequence_column not in dataframe.columns:
        # Stop because substring containment filtering cannot proceed without the column
        raise KeyError(
            f"Column '{sequence_column}' not found. Available: {list(dataframe.columns)}"
        )

    # Normalize every reference sequence once before row-by-row matching
    normalized_references = [
        normalize_sequence(
            sequence,
            strip_mods=strip_mods,
            il_equivalent=il_equivalent,
        )
        for sequence in reference_sequences
    ]
    # Remove empty reference sequences that cannot participate in substring checks
    normalized_references = [sequence for sequence in normalized_references if sequence]

    def row_matches(value: object) -> bool:
        """Return whether one query row satisfies the configured containment rule."""
        # Normalize the current row sequence for comparison
        query = normalize_sequence(
            str(value),
            strip_mods=strip_mods,
            il_equivalent=il_equivalent,
        )
        # Reject empty query sequences immediately
        if not query:
            # Return false because an empty query cannot match anything
            return False
        # Check whether the query is contained in any reference sequence
        if direction == "query_in_reference":
            # Return true when the query is a substring of any reference sequence
            return any(query in reference for reference in normalized_references)
        # Check whether any reference sequence is contained in the query
        if direction == "reference_in_query":
            # Return true when any reference sequence is a substring of the query
            return any(reference in query for reference in normalized_references)
        # Raise an error when an unsupported containment direction is requested
        raise ValueError(f"Unsupported containment direction: {direction}")

    # Build a boolean mask by applying the row matcher to every sequence value
    match_mask = dataframe[sequence_column].apply(row_matches)
    # Return only the rows that satisfy the substring containment rule
    return dataframe.loc[match_mask].copy()


def filter_table_against_proteome(
    dataframe: pd.DataFrame,
    sequence_column: str,
    index: ProteinSubstringIndex,
    config: Any,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Remove rows whose sequences are known canonical protein fragments."""
    # Count how many rows were present before proteome filtering
    input_rows = len(dataframe)
    # Build a boolean mask marking rows whose sequences are known protein fragments
    known_mask = dataframe[sequence_column].apply(
        lambda value: is_known_fragment(
            str(value),
            index,
            strip_mods=bool(getattr(config, "strip_modifications", True)),
            il_equivalent=bool(getattr(config, "treat_leucine_isoleucine_as_identical", True)),
        )
    )
    # Keep only rows that are not exact substring matches to canonical proteins
    filtered_df = dataframe.loc[~known_mask].copy()
    # Count how many rows were removed by the proteome exclusion step
    removed_rows = int(known_mask.sum())
    # Count how many novel rows remain after proteome exclusion
    retained_rows = len(filtered_df)
    # Package the per-table count summary for the run summary JSON
    summary = {
        "input_rows": input_rows,
        "removed_rows": removed_rows,
        "retained_rows": retained_rows,
    }
    # Return the filtered DataFrame and its count summary together
    return filtered_df, summary


def match_table(
    dataframe: pd.DataFrame,
    sequence_column: str,
    reference_sequences: set[str] | list[str],
    *,
    mode: str,
    direction: str = "query_in_reference",
    strip_mods: bool = True,
    il_equivalent: bool = False,
) -> pd.DataFrame:
    """Filter a table using either exact membership or substring containment."""
    # Use exact set membership when the caller requests exact matching
    if mode == "exact":
        # Ensure the reference is a set for fast membership testing
        reference_set = set(reference_sequences)
        # Return rows whose normalized sequence is exactly in the reference set
        return filter_rows_by_membership(
            dataframe,
            sequence_column,
            reference_set,
            strip_mods=strip_mods,
            il_equivalent=il_equivalent,
        )
    # Use substring containment when the caller requests substring matching
    if mode == "substring":
        # Ensure the reference is a list for ordered substring scanning
        reference_list = list(reference_sequences)
        # Return rows that satisfy the configured substring containment rule
        return filter_rows_by_containment(
            dataframe,
            sequence_column,
            reference_list,
            direction=direction,
            strip_mods=strip_mods,
            il_equivalent=il_equivalent,
        )
    # Raise an error when an unsupported matching mode is requested
    raise ValueError(f"Unsupported match mode: {mode}")
