"""Tests for is_decoy normalization — XPECTRA and QUANTIX must agree."""

from __future__ import annotations

import pandas as pd
import pytest

from msgflex.quantix.filter_psms import ensure_isdecoy


class TestIsDecoyNormalization:
    """ensure_isdecoy must produce int 0/1 regardless of input format."""

    def test_from_xxx_prefix(self):
        df = pd.DataFrame({"Protein": ["P1", "XXX_P2", "sp|P3|GENE3"]})
        out = ensure_isdecoy(df)
        assert out["is_decoy"].tolist() == [0, 1, 0]
        assert out["is_decoy"].dtype.kind == "i"

    def test_from_bool_column(self):
        df = pd.DataFrame({
            "Protein": ["P1", "P2"],
            "is_decoy": [False, True],
        })
        out = ensure_isdecoy(df)
        assert out["is_decoy"].tolist() == [0, 1]
        assert out["is_decoy"].dtype.kind == "i"

    def test_from_string_column(self):
        df = pd.DataFrame({
            "Protein": ["P1", "P2", "P3", "P4"],
            "is_decoy": ["true", "false", "1", "0"],
        })
        out = ensure_isdecoy(df)
        assert out["is_decoy"].tolist() == [1, 0, 1, 0]
        assert out["is_decoy"].dtype.kind == "i"

    def test_from_int_column_passthrough(self):
        df = pd.DataFrame({
            "Protein": ["P1", "P2"],
            "is_decoy": [0, 1],
        })
        out = ensure_isdecoy(df)
        assert out["is_decoy"].tolist() == [0, 1]
        assert out["is_decoy"].dtype.kind == "i"

    def test_handles_empty_dataframe(self):
        df = pd.DataFrame({"Protein": []})
        out = ensure_isdecoy(df)
        assert "is_decoy" in out.columns
        assert len(out) == 0

    def test_preserves_other_columns(self):
        df = pd.DataFrame({
            "Protein": ["P1", "XXX_P2"],
            "Peptide": ["K.PEP1.R", "K.PEP2.R"],
            "QValue": [0.001, 0.5],
        })
        out = ensure_isdecoy(df)
        assert "Peptide" in out.columns
        assert "QValue" in out.columns
        assert out["is_decoy"].tolist() == [0, 1]
