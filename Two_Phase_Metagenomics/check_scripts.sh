#!/bin/bash

# =============================================================================
# 🔍 Script Verification Tool
# Checks that all required scripts exist and are functional
# =============================================================================

echo "🧬 Checking Pipeline Scripts"
echo "=" * 30

SCRIPTS_DIR="scripts"

# Create scripts directory if it doesn't exist
mkdir -p "$SCRIPTS_DIR"

# List of required scripts
REQUIRED_SCRIPTS=(
    "assembly_stats.py"
    "ml_novelty_detection.py"
    "analyze_homology_novelty.py"
    "combine_novelty_results.py"
    "calculate_contig_abundance.py"
    "calculate_mag_abundance.py"
    "filter_contigs.py"
    "prepare_dastool_input.py"
    "generate_phase1_report.py"
)

echo "📜 Required Scripts Status:"
echo ""

for script in "${REQUIRED_SCRIPTS[@]}"; do
    script_path="$SCRIPTS_DIR/$script"
    
    if [ -f "$script_path" ]; then
        # Check if executable
        if [ -x "$script_path" ]; then
            echo "  ✅ $script (executable)"
        else
            echo "  ⚠️  $script (not executable - fixing...)"
            chmod +x "$script_path"
            echo "     ✅ Made executable"
        fi
        
        # Quick syntax check for Python scripts
        if [[ "$script" == *.py ]]; then
            if python3 -m py_compile "$script_path" 2>/dev/null; then
                echo "     ✅ Syntax valid"
            else
                echo "     ❌ Syntax error!"
            fi
        fi
    else
        echo "  ❌ $script (MISSING)"
    fi
done

echo ""
echo "📊 Summary:"
existing_count=$(find "$SCRIPTS_DIR" -name "*.py" -type f | wc -l)
total_scripts=${#REQUIRED_SCRIPTS[@]}

echo "  Scripts found: $existing_count/$total_scripts"

if [ "$existing_count" -eq "$total_scripts" ]; then
    echo "  Status: ✅ All scripts present"
else
    echo "  Status: ⚠️  Some scripts missing"
    echo ""
    echo "🔧 To create missing scripts:"
    echo "  1. Copy the scripts from the artifacts above"
    echo "  2. Save them in the scripts/ directory"
    echo "  3. Make them executable: chmod +x scripts/*.py"
fi

echo ""
echo "🧪 Quick Test (filter_contigs.py):"
if [ -f "$SCRIPTS_DIR/filter_contigs.py" ]; then
    if python3 "$SCRIPTS_DIR/filter_contigs.py" --help > /dev/null 2>&1; then
        echo "  ✅ filter_contigs.py help works"
    else
        echo "  ❌ filter_contigs.py has issues"
    fi
else
    echo "  ❌ filter_contigs.py not found"
fi

echo ""
echo "🔗 Quick Test (prepare_dastool_input.py):"
if [ -f "$SCRIPTS_DIR/prepare_dastool_input.py" ]; then
    if python3 "$SCRIPTS_DIR/prepare_dastool_input.py" --help > /dev/null 2>&1; then
        echo "  ✅ prepare_dastool_input.py help works"
    else
        echo "  ❌ prepare_dastool_input.py has issues"
    fi
else
    echo "  ❌ prepare_dastool_input.py not found"
fi

echo ""
echo "✨ Script verification complete!"

