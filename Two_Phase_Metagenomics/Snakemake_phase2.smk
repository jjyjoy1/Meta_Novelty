#!/usr/bin/env python3
"""
Cross-Sample Metagenomics Analysis Pipeline (Phase 2)
Snakefile for multi-sample comparative and temporal analysis with ML enhancement
"""

import os
import pandas as pd
from pathlib import Path

# =============================================================================
# CONFIGURATION
# =============================================================================

# Load configuration
configfile: "config/cross_sample_config.yaml"

# Load sample metadata
SAMPLES_DF = pd.read_csv(config["sample_metadata"], sep='\t')
SAMPLES = SAMPLES_DF['sample_id'].tolist()

# Output directory
OUTDIR = config["output_directory"]
INDIR = config["input_directory"]

# =============================================================================
# TARGET RULES
# =============================================================================

rule all:
    input:
        # Integrated abundance data
        f"{OUTDIR}/abundance/integrated_abundance_matrix.tsv",
        f"{OUTDIR}/abundance/integrated_mag_abundance_matrix.tsv",
        
        # MAG dereplication
        f"{OUTDIR}/mag_dereplication/dereplicated_mags.tsv",
        f"{OUTDIR}/mag_dereplication/mag_clusters.tsv",
        
        # ML analysis
        f"{OUTDIR}/ml_analysis/vae_embeddings.npy",
        f"{OUTDIR}/ml_analysis/vae_model_metrics.json",
        f"{OUTDIR}/ml_analysis/isolation_forest_results.tsv",
        
        # Novelty analysis
        f"{OUTDIR}/novelty_analysis/pan_novelty_results.tsv",
        f"{OUTDIR}/novelty_analysis/novelty_emergence_patterns.tsv",
        
        # Temporal analysis
        f"{OUTDIR}/temporal_analysis/temporal_trends.tsv",
        f"{OUTDIR}/temporal_analysis/changepoint_detection.tsv",
        
        # Pathogen detection
        f"{OUTDIR}/pathogen_detection/pathogen_risk_assessment.tsv",
        f"{OUTDIR}/pathogen_detection/high_risk_alerts.tsv",
        
        # Comprehensive report
        f"{OUTDIR}/reports/comprehensive_analysis_report.html",
        f"{OUTDIR}/reports/executive_summary.pdf"

# Quick analysis targets
rule abundance_only:
    input:
        f"{OUTDIR}/abundance/integrated_abundance_matrix.tsv",
        f"{OUTDIR}/abundance/integrated_mag_abundance_matrix.tsv"

rule ml_only:
    input:
        f"{OUTDIR}/ml_analysis/vae_embeddings.npy",
        f"{OUTDIR}/ml_analysis/isolation_forest_results.tsv"

rule novelty_only:
    input:
        f"{OUTDIR}/novelty_analysis/pan_novelty_results.tsv",
        f"{OUTDIR}/novelty_analysis/novelty_emergence_patterns.tsv"

rule temporal_only:
    input:
        f"{OUTDIR}/temporal_analysis/temporal_trends.tsv",
        f"{OUTDIR}/temporal_analysis/changepoint_detection.tsv"

rule pathogen_only:
    input:
        f"{OUTDIR}/pathogen_detection/pathogen_risk_assessment.tsv",
        f"{OUTDIR}/pathogen_detection/high_risk_alerts.tsv"

# =============================================================================
# ABUNDANCE DATA INTEGRATION
# =============================================================================

rule integrate_abundance_data:
    input:
        abundance_files = expand(f"{INDIR}/{{sample}}/abundance/contig_abundance.tsv", sample=SAMPLES),
        mag_abundance_files = expand(f"{INDIR}/{{sample}}/abundance/mag_abundance.tsv", sample=SAMPLES),
        metadata = config["sample_metadata"]
    output:
        abundance_matrix = f"{OUTDIR}/abundance/integrated_abundance_matrix.tsv",
        mag_abundance_matrix = f"{OUTDIR}/abundance/integrated_mag_abundance_matrix.tsv",
        abundance_summary = f"{OUTDIR}/abundance/abundance_integration_summary.json"
    threads: config["resources"]["default_threads"]
    resources:
        mem_mb = config["resources"]["memory_limits"]["abundance_calculation"] * 1024
    log:
        f"{OUTDIR}/logs/integrate_abundance_data.log"
    script:
        "scripts/integrate_abundance_data.py"

# =============================================================================
# MAG DEREPLICATION
# =============================================================================

rule mag_dereplication:
    input:
        mag_files = expand(f"{INDIR}/{{sample}}/binning/final_bins", sample=SAMPLES),
        mag_quality = expand(f"{INDIR}/{{sample}}/binning/checkm2_results.tsv", sample=SAMPLES),
        abundance_matrix = f"{OUTDIR}/abundance/integrated_mag_abundance_matrix.tsv"
    output:
        dereplicated_mags = f"{OUTDIR}/mag_dereplication/dereplicated_mags.tsv",
        mag_clusters = f"{OUTDIR}/mag_dereplication/mag_clusters.tsv",
        representative_mags = directory(f"{OUTDIR}/mag_dereplication/representative_mags"),
        ani_matrix = f"{OUTDIR}/mag_dereplication/ani_matrix.tsv"
    params:
        ani_threshold = config["mag_dereplication"]["ani_threshold"],
        min_alignment_fraction = config["mag_dereplication"]["min_alignment_fraction"],
        quality_weight = config["mag_dereplication"]["quality_weight"],
        size_weight = config["mag_dereplication"]["size_weight"]
    threads: config["mag_dereplication"]["threads"]
    resources:
        mem_mb = config["resources"]["memory_limits"]["mag_dereplication"] * 1024
    log:
        f"{OUTDIR}/logs/mag_dereplication.log"
    script:
        "scripts/mag_dereplication.py"

# =============================================================================
# MACHINE LEARNING ANALYSIS
# =============================================================================

rule vae_community_embedding:
    input:
        abundance_matrix = f"{OUTDIR}/abundance/integrated_abundance_matrix.tsv",
        metadata = config["sample_metadata"]
    output:
        embeddings = f"{OUTDIR}/ml_analysis/vae_embeddings.npy",
        model_file = f"{OUTDIR}/ml_analysis/vae_model.pkl",
        training_history = f"{OUTDIR}/ml_analysis/vae_training_history.json",
        model_metrics = f"{OUTDIR}/ml_analysis/vae_model_metrics.json",
        reconstruction_error = f"{OUTDIR}/ml_analysis/reconstruction_errors.tsv"
    params:
        latent_dim = config["vae_settings"]["latent_dim"],
        hidden_dims = config["vae_settings"]["hidden_dims"],
        batch_size = config["vae_settings"]["batch_size"],
        learning_rate = config["vae_settings"]["learning_rate"],
        num_epochs = config["vae_settings"]["num_epochs"],
        early_stopping_patience = config["vae_settings"]["early_stopping_patience"],
        beta = config["vae_settings"]["beta"],
        dropout_rate = config["vae_settings"]["dropout_rate"]
    threads: config["resources"]["default_threads"]
    resources:
        mem_mb = config["resources"]["memory_limits"]["vae_training"] * 1024
    log:
        f"{OUTDIR}/logs/vae_community_embedding.log"
    script:
        "scripts/vae_community_embedding.py"

rule isolation_forest_cross_sample:
    input:
        abundance_matrix = f"{OUTDIR}/abundance/integrated_abundance_matrix.tsv",
        vae_embeddings = f"{OUTDIR}/ml_analysis/vae_embeddings.npy",
        novelty_data = f"{OUTDIR}/novelty_analysis/pan_novelty_results.tsv",
        metadata = config["sample_metadata"]
    output:
        anomaly_scores = f"{OUTDIR}/ml_analysis/isolation_forest_results.tsv",
        sample_anomalies = f"{OUTDIR}/ml_analysis/sample_anomaly_scores.tsv",
        feature_importance = f"{OUTDIR}/ml_analysis/feature_importance.tsv",
        outlier_analysis = f"{OUTDIR}/ml_analysis/outlier_analysis.json"
    params:
        n_estimators = config["isolation_forest"]["n_estimators"],
        contamination = config["isolation_forest"]["contamination"],
        features = config["isolation_forest"]["features"],
        scaling_method = config["isolation_forest"]["scaling_method"]
    threads: config["resources"]["default_threads"]
    log:
        f"{OUTDIR}/logs/isolation_forest_cross_sample.log"
    script:
        "scripts/isolation_forest_cross_sample.py"

# =============================================================================
# NOVELTY ANALYSIS
# =============================================================================

rule pan_novelty_analysis:
    input:
        novelty_files = expand(f"{INDIR}/{{sample}}/novelty/combined_novelty_results.tsv", sample=SAMPLES),
        abundance_matrix = f"{OUTDIR}/abundance/integrated_abundance_matrix.tsv",
        metadata = config["sample_metadata"]
    output:
        pan_novelty_results = f"{OUTDIR}/novelty_analysis/pan_novelty_results.tsv",
        novelty_emergence_patterns = f"{OUTDIR}/novelty_analysis/novelty_emergence_patterns.tsv",
        core_novelty = f"{OUTDIR}/novelty_analysis/core_novelty_taxa.tsv",
        rare_novelty = f"{OUTDIR}/novelty_analysis/rare_novelty_taxa.tsv",
        novelty_summary = f"{OUTDIR}/novelty_analysis/novelty_analysis_summary.json"
    params:
        min_samples = config["pan_novelty"]["min_samples"],
        thresholds = config["pan_novelty"]["thresholds"],
        emergence_detection = config["pan_novelty"]["emergence_detection"],
        core_threshold = config["pan_novelty"]["core_novelty_threshold"],
        rare_threshold = config["pan_novelty"]["rare_novelty_threshold"]
    threads: config["resources"]["default_threads"]
    log:
        f"{OUTDIR}/logs/pan_novelty_analysis.log"
    script:
        "scripts/pan_novelty_analysis.py"

# =============================================================================
# TEMPORAL ANALYSIS
# =============================================================================

rule temporal_analysis:
    input:
        abundance_matrix = f"{OUTDIR}/abundance/integrated_abundance_matrix.tsv",
        novelty_data = f"{OUTDIR}/novelty_analysis/pan_novelty_results.tsv",
        metadata = config["sample_metadata"]
    output:
        temporal_trends = f"{OUTDIR}/temporal_analysis/temporal_trends.tsv",
        changepoint_detection = f"{OUTDIR}/temporal_analysis/changepoint_detection.tsv",
        seasonal_decomposition = f"{OUTDIR}/temporal_analysis/seasonal_decomposition.tsv",
        trend_significance = f"{OUTDIR}/temporal_analysis/trend_significance_tests.tsv",
        temporal_plots = f"{OUTDIR}/temporal_analysis/temporal_visualization.html"
    params:
        enabled = config["temporal_analysis"]["enabled"],
        time_unit = config["temporal_analysis"]["time_unit"],
        smoothing = config["temporal_analysis"]["smoothing"],
        trend_detection = config["temporal_analysis"]["trend_detection"],
        changepoint_detection = config["temporal_analysis"]["changepoint_detection"],
        seasonal_decomposition = config["temporal_analysis"]["seasonal_decomposition"]
    threads: config["resources"]["default_threads"]
    log:
        f"{OUTDIR}/logs/temporal_analysis.log"
    script:
        "scripts/temporal_analysis.py"

# =============================================================================
# PATHOGEN DETECTION
# =============================================================================

rule pathogen_detection:
    input:
        abundance_matrix = f"{OUTDIR}/abundance/integrated_abundance_matrix.tsv",
        novelty_data = f"{OUTDIR}/novelty_analysis/pan_novelty_results.tsv",
        contigs = expand(f"{INDIR}/{{sample}}/assembly/contigs.fasta", sample=SAMPLES),
        metadata = config["sample_metadata"]
    output:
        pathogen_risk_assessment = f"{OUTDIR}/pathogen_detection/pathogen_risk_assessment.tsv",
        high_risk_alerts = f"{OUTDIR}/pathogen_detection/high_risk_alerts.tsv",
        virulence_analysis = f"{OUTDIR}/pathogen_detection/virulence_factor_analysis.tsv",
        resistance_analysis = f"{OUTDIR}/pathogen_detection/antibiotic_resistance_analysis.tsv",
        pathogen_timeline = f"{OUTDIR}/pathogen_detection/pathogen_detection_timeline.tsv"
    params:
        databases = config["pathogen_detection"]["databases"],
        blast_params = config["pathogen_detection"]["blast_params"],
        risk_scoring = config["pathogen_detection"]["risk_scoring"],
        alert_thresholds = config["pathogen_detection"]["alert_thresholds"],
        known_pathogens = config["pathogen_detection"]["known_pathogens"]
    threads: config["resources"]["default_threads"]
    log:
        f"{OUTDIR}/logs/pathogen_detection.log"
    script:
        "scripts/pathogen_detection.py"

# =============================================================================
# COMPREHENSIVE REPORTING
# =============================================================================

rule generate_comprehensive_report:
    input:
        abundance_matrix = f"{OUTDIR}/abundance/integrated_abundance_matrix.tsv",
        mag_abundance_matrix = f"{OUTDIR}/abundance/integrated_mag_abundance_matrix.tsv",
        dereplicated_mags = f"{OUTDIR}/mag_dereplication/dereplicated_mags.tsv",
        vae_embeddings = f"{OUTDIR}/ml_analysis/vae_embeddings.npy",
        vae_metrics = f"{OUTDIR}/ml_analysis/vae_model_metrics.json",
        isolation_forest_results = f"{OUTDIR}/ml_analysis/isolation_forest_results.tsv",
        pan_novelty_results = f"{OUTDIR}/novelty_analysis/pan_novelty_results.tsv",
        novelty_emergence = f"{OUTDIR}/novelty_analysis/novelty_emergence_patterns.tsv",
        temporal_trends = f"{OUTDIR}/temporal_analysis/temporal_trends.tsv",
        pathogen_risk = f"{OUTDIR}/pathogen_detection/pathogen_risk_assessment.tsv",
        high_risk_alerts = f"{OUTDIR}/pathogen_detection/high_risk_alerts.tsv",
        metadata = config["sample_metadata"]
    output:
        html_report = f"{OUTDIR}/reports/comprehensive_analysis_report.html",
        pdf_summary = f"{OUTDIR}/reports/executive_summary.pdf",
        interactive_dashboard = f"{OUTDIR}/reports/interactive_dashboard.html",
        analysis_summary = f"{OUTDIR}/reports/analysis_summary.json"
    params:
        formats = config["reporting"]["formats"],
        sections = config["reporting"]["sections"],
        visualization = config["reporting"]["visualization"],
        interactive_plots = config["reporting"]["interactive_plots"]
    threads: config["resources"]["default_threads"]
    resources:
        mem_mb = config["resources"]["memory_limits"]["report_generation"] * 1024
    log:
        f"{OUTDIR}/logs/generate_comprehensive_report.log"
    script:
        "scripts/generate_comprehensive_report.py"

# =============================================================================
# QUALITY CONTROL AND VALIDATION
# =============================================================================

rule quality_control_check:
    input:
        metadata = config["sample_metadata"],
        abundance_matrix = f"{OUTDIR}/abundance/integrated_abundance_matrix.tsv"
    output:
        qc_report = f"{OUTDIR}/quality_control/cross_sample_qc_report.json",
        outlier_detection = f"{OUTDIR}/quality_control/outlier_samples.tsv",
        data_completeness = f"{OUTDIR}/quality_control/data_completeness_summary.tsv"
    params:
        min_samples = config["quality_control"]["min_samples"],
        min_reads_per_sample = config["quality_control"]["min_reads_per_sample"],
        max_missing_data = config["quality_control"]["max_missing_data"],
        outlier_detection = config["quality_control"]["outlier_detection"]
    threads: config["resources"]["default_threads"]
    log:
        f"{OUTDIR}/logs/quality_control_check.log"
    script:
        "scripts/quality_control_cross_sample.py"

# =============================================================================
# UTILITY RULES
# =============================================================================

rule clean_intermediate:
    shell:
        """
        if [ "{config[resources][keep_intermediate]}" = "false" ]; then
            find {OUTDIR} -name "*.tmp" -delete
            find {OUTDIR} -name "temp_*" -delete
            echo "Intermediate files cleaned"
        else
            echo "Keeping intermediate files as configured"
        fi
        """

rule create_output_directories:
    output:
        abundance_dir = directory(f"{OUTDIR}/abundance"),
        mag_dereplic_dir = directory(f"{OUTDIR}/mag_dereplication"),
        ml_analysis_dir = directory(f"{OUTDIR}/ml_analysis"),
        novelty_dir = directory(f"{OUTDIR}/novelty_analysis"),
        temporal_dir = directory(f"{OUTDIR}/temporal_analysis"),
        pathogen_dir = directory(f"{OUTDIR}/pathogen_detection"),
        reports_dir = directory(f"{OUTDIR}/reports"),
        qc_dir = directory(f"{OUTDIR}/quality_control"),
        logs_dir = directory(f"{OUTDIR}/logs")
    shell:
        """
        mkdir -p {output}
        """

# =============================================================================
# ERROR HANDLING AND CHECKPOINTS
# =============================================================================

checkpoint validate_inputs:
    input:
        metadata = config["sample_metadata"]
    output:
        validation_report = f"{OUTDIR}/validation/input_validation.json"
    run:
        import json
        import os
        
        validation_results = {
            "sample_count": len(SAMPLES),
            "metadata_valid": os.path.exists(config["sample_metadata"]),
            "single_sample_results_available": [],
            "missing_files": []
        }
        
        # Check if single-sample results exist
        for sample in SAMPLES:
            sample_dir = f"{INDIR}/{sample}"
            required_files = [
                f"{sample_dir}/abundance/contig_abundance.tsv",
                f"{sample_dir}/novelty/combined_novelty_results.tsv"
            ]
            
            missing_files = [f for f in required_files if not os.path.exists(f)]
            if missing_files:
                validation_results["missing_files"].extend(missing_files)
            else:
                validation_results["single_sample_results_available"].append(sample)
        
        # Ensure output directory exists
        os.makedirs(os.path.dirname(output.validation_report), exist_ok=True)
        
        # Write validation report
        with open(output.validation_report, 'w') as f:
            json.dump(validation_results, f, indent=2)

# =============================================================================
# PERFORMANCE MONITORING
# =============================================================================

if config["advanced"]["performance_monitoring"]["enabled"]:
    rule monitor_performance:
        input:
            f"{OUTDIR}/reports/comprehensive_analysis_report.html"
        output:
            performance_report = f"{OUTDIR}/performance/pipeline_performance.json",
            resource_usage = f"{OUTDIR}/performance/resource_usage.tsv"
        shell:
            """
            python scripts/monitor_pipeline_performance.py \
                --log-dir {OUTDIR}/logs \
                --output-performance {output.performance_report} \
                --output-resources {output.resource_usage}
            """

# =============================================================================
# WORKFLOW HOOKS
# =============================================================================

onsuccess:
    print("="*80)
    print("🎉 Cross-sample analysis pipeline completed successfully!")
    print(f"📊 Results available in: {OUTDIR}")
    print(f"📋 Comprehensive report: {OUTDIR}/reports/comprehensive_analysis_report.html")
    print(f"🚨 High-risk alerts: {OUTDIR}/pathogen_detection/high_risk_alerts.tsv")
    print("="*80)

onerror:
    print("="*80)
    print("❌ Cross-sample analysis pipeline failed!")
    print(f"📋 Check logs in: {OUTDIR}/logs/")
    print("🔧 Common issues and solutions:")
    print("   - Ensure all single-sample results are complete")
    print("   - Check memory allocation for large datasets")
    print("   - Verify sample metadata format")
    print("="*80)

