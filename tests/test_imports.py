"""
Import-chain smoke tests.
"""

from __future__ import annotations


class TestCoreImports:

    def test_project_paths(self):
        from msgflex.core import ProjectPaths
        assert ProjectPaths is not None

    def test_exception_hierarchy(self):
        from msgflex.core import (
            ConfigurationError,
            MsgflexError,
            QuantixError,
            ToolsNotFoundError,
            SparxError,
            XpectraError,
        )
        assert issubclass(ConfigurationError, MsgflexError)
        assert issubclass(SparxError, MsgflexError)
        assert issubclass(XpectraError, MsgflexError)
        assert issubclass(QuantixError, MsgflexError)
        assert issubclass(ToolsNotFoundError, ConfigurationError)

    def test_logging(self):
        from msgflex.core import get_logger, setup_logging
        assert callable(setup_logging)
        assert callable(get_logger)


class TestXpectraImports:

    def test_package_run(self):
        from msgflex.xpectra import run
        assert callable(run)

    def test_rxflow(self):
        from msgflex.xpectra.rxflow import rxflow
        assert callable(rxflow)

    def test_model_apis(self):
        from msgflex.xpectra.model import (
            FeatureEngineer,
            cross_validated_rescoring,
            select_model_features,
        )
        assert callable(cross_validated_rescoring)

    def test_features(self):
        from msgflex.xpectra.features import (
            IonCalculator,
            PeptideParser,
            SpectralFeatureExtractor,
            extract_spectral_features,
        )
        assert callable(extract_spectral_features)

    def test_metrics(self):
        from msgflex.xpectra.metrics import (
            compare_gains,
            estimate_fdr,
            estimate_pep_fdr,
        )
        assert callable(estimate_fdr)

    def test_rt_calibrator(self):
        from msgflex.xpectra.xrtpred import RTCalibrator
        assert RTCalibrator is not None

    def test_utils(self):
        from msgflex.xpectra.utils import (
            _compute_sample_weights,
            _infer_is_decoy,
        )


class TestQuantixImports:

    def test_package_run(self):
        from msgflex.quantix import run
        assert callable(run)

    def test_pipeline(self):
        from msgflex.quantix.pipeline import QuantixPipeline
        assert QuantixPipeline is not None

    def test_prepare_sics(self):
        from msgflex.quantix.prepare_sics import prepare_SICs
        assert callable(prepare_SICs)

    def test_merge_shared(self):
        from msgflex.quantix.merge_shared import append_syn_to_fdr
        assert callable(append_syn_to_fdr)

    def test_filter_psms(self):
        from msgflex.quantix.filter_psms import (
            ensure_isdecoy,
            process_fdrdir,
        )

    def test_peptide_crosstab(self):
        from msgflex.quantix.peptide_crosstab import merge_peptide_ctab
        assert callable(merge_peptide_ctab)

    def test_annotator(self):
        from msgflex.quantix.annotator import (
            build_protein_info_map,
            get_protein_function,
        )

    def test_coverage(self):
        from msgflex.quantix.coverage import (
            compute_coverage_for_sample,
            detect_sample_columns,
            ensure_required_columns,
            load_fasta_for_samples,
            prepare_subtable_for_sample,
        )

    def test_rollup(self):
        from msgflex.quantix.rollup import rollup_from_annotated, run_rollup
        assert callable(run_rollup)


class TestPipelineImports:
    """These depend on user-provided modules (masic_merger, dbcurator,
    MASIC_wrapper, inputfile_generator) which live in msgflex/modules/
    in the real installation. If those aren't present we skip —
    the tests are about our own integration surface, not the user's files.
    """

    def test_conventional(self):
        try:
            from msgflex.pipelines.conventional import run as conv_run
            assert callable(conv_run)
        except ImportError as e:
            import pytest
            pytest.skip(f"Sparx modules not bundled: {e}")

    def test_binning(self):
        try:
            from msgflex.pipelines.binning import run as bin_run
            assert callable(bin_run)
        except ImportError as e:
            import pytest
            pytest.skip(f"Sparx modules not bundled: {e}")


class TestOrchestratorImports:

    def test_run_full(self):
        from msgflex.orchestrator import run_full
        assert callable(run_full)


class TestCliImports:

    def test_app(self):
        from msgflex.cli import app, main
        assert callable(main)
