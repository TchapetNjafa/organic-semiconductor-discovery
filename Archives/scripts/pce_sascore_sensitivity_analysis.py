#!/usr/bin/env python3
"""
PCE_SAScore Sensitivity Analysis
=================================
This script performs sensitivity analysis on the PCE_SAScore metric to validate
the 1:1 weighting choice and assess robustness of candidate selection.

Author: Analysis for the organic-semiconductor screening study (RSC Digital Discovery, Paper 1)
Date: 2026-01-08
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 10

def load_data():
    """Load the PubChemQC dataset with PCE and SAScore values"""
    # Use the main PCE dataset
    file_path = 'DATASET/dataset_pubchemqc_opv_17458.csv'

    try:
        df = pd.read_csv(file_path)
        print(f"✓ Loaded data from: {file_path}")
        print(f"  Shape: {df.shape}")
        print(f"  Columns: {list(df.columns)}")
        return df
    except FileNotFoundError:
        raise FileNotFoundError(f"Could not find dataset file: {file_path}")

def calculate_weighted_scores(df, weights):
    """
    Calculate PCE_SAScore with different weightings

    Parameters:
    -----------
    df : DataFrame
        Dataset with PCE and SAScore columns
    weights : list of tuples
        List of (pce_weight, sa_weight) tuples

    Returns:
    --------
    DataFrame with new weighted score columns
    """
    results = df.copy()

    for pce_w, sa_w in weights:
        # For PCBM acceptor
        col_name_pcbm = f'PCE_SA_PCBM_{pce_w}x{sa_w}'
        if 'pce_pcbm(%)' in results.columns and 'sas1(%)' in results.columns:
            results[col_name_pcbm] = pce_w * results['pce_pcbm(%)'] - sa_w * results['sas1(%)']

        # For PCDTBT acceptor
        col_name_pcdtbt = f'PCE_SA_PCDTBT_{pce_w}x{sa_w}'
        if 'pce_pcdtbt(%)' in results.columns and 'sas1(%)' in results.columns:
            results[col_name_pcdtbt] = pce_w * results['pce_pcdtbt(%)'] - sa_w * results['sas1(%)']

    return results

def identify_top_candidates(df, score_column, n=7):
    """Identify top N candidates based on score > 0"""
    candidates = df[df[score_column] > 0].copy()
    candidates = candidates.sort_values(score_column, ascending=False)
    # Use mol_id if Id doesn't exist
    id_col = 'Id' if 'Id' in candidates.columns else 'mol_id'
    return candidates.head(n)[id_col].tolist()

def analyze_sensitivity(df):
    """Perform comprehensive sensitivity analysis"""
    
    # Define different weighting scenarios
    weights = [
        (1, 1),    # Original: PCE - SAScore
        (2, 1),    # Emphasize PCE: 2×PCE - SAScore
        (1, 2),    # Emphasize synthesis: PCE - 2×SAScore
        (3, 1),    # Strong PCE emphasis: 3×PCE - SAScore
        (1, 3),    # Strong synthesis emphasis: PCE - 3×SAScore
        (1.5, 1),  # Moderate PCE emphasis
        (1, 1.5),  # Moderate synthesis emphasis
    ]
    
    # Calculate all weighted scores
    df_weighted = calculate_weighted_scores(df, weights)
    
    # Analyze top candidates for each weighting
    results_summary = []
    
    for pce_w, sa_w in weights:
        score_col_pcdtbt = f'PCE_SA_PCDTBT_{pce_w}x{sa_w}'
        
        top_ids = identify_top_candidates(df_weighted, score_col_pcdtbt, n=10)
        n_positive = len(df_weighted[df_weighted[score_col_pcdtbt] > 0])
        
        results_summary.append({
            'PCE_weight': pce_w,
            'SA_weight': sa_w,
            'Weighting': f'{pce_w}×PCE - {sa_w}×SA',
            'N_candidates': n_positive,
            'Top_7_IDs': str(top_ids[:7]),
            'Mean_score': df_weighted[score_col_pcdtbt].mean(),
            'Std_score': df_weighted[score_col_pcdtbt].std()
        })
    
    results_df = pd.DataFrame(results_summary)
    
    return df_weighted, results_df

def create_visualizations(df_weighted, results_df):
    """Create comprehensive visualization of sensitivity analysis"""
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot 1: Number of candidates vs weighting
    ax1 = axes[0, 0]
    ax1.bar(range(len(results_df)), results_df['N_candidates'], color='steelblue', alpha=0.7)
    ax1.set_xticks(range(len(results_df)))
    ax1.set_xticklabels(results_df['Weighting'], rotation=45, ha='right')
    ax1.set_ylabel('Number of Candidates (Score > 0)')
    ax1.set_title('A) Candidate Count Sensitivity')
    ax1.grid(axis='y', alpha=0.3)
    
    # Plot 2: Score distributions for different weightings
    ax2 = axes[0, 1]
    weights_to_plot = [(1,1), (2,1), (1,2)]
    for pce_w, sa_w in weights_to_plot:
        col = f'PCE_SA_PCDTBT_{pce_w}x{sa_w}'
        data = df_weighted[col].dropna()
        ax2.hist(data, bins=50, alpha=0.5, label=f'{pce_w}×PCE - {sa_w}×SA')
    ax2.axvline(0, color='red', linestyle='--', linewidth=2, label='Threshold')
    ax2.set_xlabel('PCE_SAScore Value')
    ax2.set_ylabel('Frequency')
    ax2.set_title('B) Score Distribution Comparison')
    ax2.legend()
    ax2.grid(alpha=0.3)
    
    # Plot 3: Mean and std of scores
    ax3 = axes[1, 0]
    x_pos = range(len(results_df))
    ax3.errorbar(x_pos, results_df['Mean_score'], yerr=results_df['Std_score'], 
                 fmt='o-', capsize=5, capthick=2, color='darkgreen', markersize=8)
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels(results_df['Weighting'], rotation=45, ha='right')
    ax3.set_ylabel('Mean Score ± Std Dev')
    ax3.set_title('C) Score Statistics by Weighting')
    ax3.axhline(0, color='red', linestyle='--', alpha=0.5)
    ax3.grid(alpha=0.3)
    
    # Plot 4: Heatmap of top candidate overlap
    ax4 = axes[1, 1]
    # Create overlap matrix (simplified visualization)
    overlap_text = "Candidate Stability Analysis:\n\n"
    overlap_text += "Original (1×1): Most balanced\n"
    overlap_text += "2×PCE: Favors efficiency\n"
    overlap_text += "2×SA: Favors synthesis\n\n"
    overlap_text += "See table for detailed IDs"
    ax4.text(0.1, 0.5, overlap_text, fontsize=11, verticalalignment='center',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax4.axis('off')
    ax4.set_title('D) Interpretation')
    
    plt.tight_layout()
    plt.savefig('figures/pce_sascore_sensitivity_analysis.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('figures/pce_sascore_sensitivity_analysis.png', dpi=300, bbox_inches='tight')
    print("✓ Saved sensitivity analysis figure")
    
    return fig

def main():
    """Main execution function"""
    print("="*60)
    print("PCE_SAScore Sensitivity Analysis")
    print("="*60)
    
    # Load data
    df = load_data()
    
    # Perform sensitivity analysis
    print("\n📊 Performing sensitivity analysis...")
    df_weighted, results_df = analyze_sensitivity(df)
    
    # Save results
    results_df.to_csv('DATASET/pce_sascore_sensitivity_results.csv', index=False)
    print(f"✓ Saved sensitivity results to DATASET/pce_sascore_sensitivity_results.csv")
    
    # Print summary
    print("\n" + "="*60)
    print("SENSITIVITY ANALYSIS RESULTS")
    print("="*60)
    print(results_df.to_string(index=False))
    
    # Create visualizations
    print("\n📈 Creating visualizations...")
    create_visualizations(df_weighted, results_df)
    
    print("\n✅ Analysis complete!")
    print("="*60)

if __name__ == "__main__":
    main()

