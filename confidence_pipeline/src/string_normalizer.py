"""Peptide sequence normalization utilities for confidence-tier matching."""

from __future__ import annotations

# Import the regular expression module for modification stripping
import re

# Compile a pattern that matches bracketed modification tokens like [UNIMOD:4]
MODIFICATION_PATTERN = re.compile(r"\[[^\]]+\]")

# Compile a pattern that removes all whitespace characters including newlines
WHITESPACE_PATTERN = re.compile(r"[\s\r\n\t]+")

# Compile a pattern that validates standard uppercase amino acid letters
AMINO_ACID_PATTERN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


def strip_modifications(sequence: str) -> str:
    """Remove bracketed mass-spectrometry modification annotations from a sequence."""
    # Convert the input value to a sanitized raw string
    text = sanitize_raw_sequence(sequence)
    # Remove every bracketed modification token from the peptide text
    stripped = MODIFICATION_PATTERN.sub("", text)
    # Return the modification-free sequence string
    return stripped


def sanitize_raw_sequence(sequence: str) -> str:
    """Remove surrounding quotes and all whitespace from a raw peptide string."""
    # Convert the input value to a plain string for defensive cleaning
    text = str(sequence)
    # Remove leading and trailing ASCII whitespace from the string
    text = text.strip()
    # Remove one layer of surrounding double quotes when present
    if len(text) >= 2 and text[0] == '"' and text[-1] == '"':
        # Strip the outer double-quote characters from the peptide text
        text = text[1:-1]
    # Remove one layer of surrounding single quotes when present
    if len(text) >= 2 and text[0] == "'" and text[-1] == "'":
        # Strip the outer single-quote characters from the peptide text
        text = text[1:-1]
    # Delete every internal whitespace, tab, carriage-return, and newline character
    text = WHITESPACE_PATTERN.sub("", text)
    # Return the sanitized raw peptide string
    return text


def unify_il(sequence: str) -> str:
    """Map isoleucine to leucine for isobaric mass-spectrometry matching."""
    # Convert the sequence to uppercase text before amino acid replacement
    upper_text = str(sequence).upper()
    # Replace every isoleucine with leucine for isobaric comparison
    normalized = upper_text.replace("I", "L")
    # Return the I/L-normalized sequence
    return normalized


def normalize_sequence(
    sequence: str,
    *,
    strip_mods: bool = True,
    il_equivalent: bool = False,
) -> str:
    """Normalize a peptide sequence for consistent string comparison."""
    # Start with the raw sequence converted to a sanitized string
    text = sanitize_raw_sequence(sequence)
    # Optionally strip modification annotations before further normalization
    if strip_mods:
        # Remove bracketed tokens such as UNIMOD labels from the sequence
        text = strip_modifications(text)
    # Convert the working sequence to uppercase letters
    text = text.upper()
    # Optionally collapse isoleucine into leucine for isobaric matching
    if il_equivalent:
        # Replace isoleucine with leucine when the caller requests I/L equivalence
        text = unify_il(text)
    # Return the fully normalized peptide sequence
    return text


def is_valid_peptide_sequence(sequence: str) -> bool:
    """Return whether a normalized sequence contains only standard amino acids."""
    # Reject empty or missing sequences immediately
    if not sequence:
        # Return false because an empty string cannot be a valid peptide
        return False
    # Check that every character is a standard uppercase amino acid letter
    return bool(AMINO_ACID_PATTERN.match(sequence))
