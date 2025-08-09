"""
I/O utilities for saving and loading analysis results.
"""

import pandas as pd
from typing import Dict


def save_statistics(protein_stats: Dict, output_file: str = 'domain_statistics.csv'):
    """
    Save protein domain distribution and lipid statistics.
    
    Parameters
    ----------
    protein_stats : dict
        Dictionary containing protein statistics
    output_file : str
        Path to output CSV file
    """
    stats_data = []
    lipid_types = ['CHOL', 'DIPC', 'DPSM', 'DPG3']
    
    for protein_name, domains in protein_stats.items():
        total_frames = domains['total']
        if total_frames > 0:
            for domain in ['core_cs', 'cs', 'd']:
                if domain in domains:
                    domain_stats = domains[domain]
                    frames_in_domain = domain_stats['frames']
                    
                    if frames_in_domain > 0:
                        # Basic statistics
                        stats_row = {
                            'Protein': protein_name,
                            'Domain': domain.upper(),
                            'Frames_in_Domain': frames_in_domain,
                            'Domain_Percentage': (frames_in_domain / total_frames) * 100,
                            'Total_Frames': total_frames
                        }
                        
                        # Lipid statistics
                        for lipid_type in lipid_types:
                            count_key = f'{lipid_type}_sum'
                            density_key = f'{lipid_type}_density_sum'
                            coloc_key = f'{lipid_type}_gm3_coloc_sum'
                            
                            if count_key in domain_stats and density_key in domain_stats:
                                avg_count = domain_stats[count_key] / frames_in_domain
                                avg_density = domain_stats[density_key] / frames_in_domain
                                
                                stats_row[f'Average_{lipid_type}_Count'] = avg_count
                                stats_row[f'Average_{lipid_type}_Density'] = avg_density
                                
                                # GM3 colocalization data
                                if coloc_key in domain_stats:
                                    avg_coloc = domain_stats[coloc_key] / frames_in_domain
                                    stats_row[f'Average_{lipid_type}_GM3_Colocalization'] = avg_coloc
                        
                        stats_data.append(stats_row)
    
    if stats_data:
        # Save to CSV
        df = pd.DataFrame(stats_data)
        df.to_csv(output_file, index=False)
        print(f"\nStatistics saved to {output_file}")
        
        # Display summary
        print_statistics_summary(df, protein_stats)
    else:
        print("No statistics to save.")


def print_statistics_summary(df: pd.DataFrame, protein_stats: Dict):
    """
    Print a summary of the statistics.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing statistics
    protein_stats : dict
        Dictionary containing protein statistics
    """
    print("\nDomain Distribution Summary:")
    
    for protein_name in protein_stats.keys():
        print(f"\n{protein_name}:")
        protein_data = df[df['Protein'] == protein_name]
        
        for domain in ['CORE_CS', 'CS', 'D']:
            domain_data = protein_data[protein_data['Domain'] == domain]
            if not domain_data.empty:
                if domain == 'CORE_CS':
                    domain_label = "Core-CS rich"
                elif domain == 'CS':
                    domain_label = "CS rich region"
                else:
                    domain_label = "D rich region"
                    
                print(f"\n  {domain_label}:")
                print(f"    Percentage: {domain_data['Domain_Percentage'].iloc[0]:.1f}%")
                print(f"    Frames: {int(domain_data['Frames_in_Domain'].iloc[0])}")
                
                # Lipid statistics
                for lipid_type in ['CHOL', 'DIPC', 'DPSM', 'DPG3']:
                    count_col = f'Average_{lipid_type}_Count'
                    density_col = f'Average_{lipid_type}_Density'
                    coloc_col = f'Average_{lipid_type}_GM3_Colocalization'
                    
                    if count_col in domain_data.columns and density_col in domain_data.columns:
                        print(f"    Average {lipid_type} count: {domain_data[count_col].iloc[0]:.1f}")
                        print(f"    Average {lipid_type} density: {domain_data[density_col].iloc[0]:.4f}/nm²")
                        
                        if coloc_col in domain_data.columns:
                            coloc_percentage = domain_data[coloc_col].iloc[0] * 100
                            print(f"    {lipid_type} GM3 colocalization: {coloc_percentage:.1f}%")
    
    # Overall statistics
    print("\nOverall Domain Distribution:")
    for domain in ['CORE_CS', 'CS', 'D']:
        domain_data = df[df['Domain'] == domain]
        if not domain_data.empty:
            avg_percentage = domain_data['Domain_Percentage'].mean()
            total_frames = domain_data['Frames_in_Domain'].sum()
            
            if domain == 'CORE_CS':
                domain_label = "Core-CS rich regions"
            elif domain == 'CS':
                domain_label = "CS rich regions"
            else:
                domain_label = "D rich regions"
                
            print(f"  {domain_label}: {avg_percentage:.1f}% ({total_frames} frames)")
            
            # Average lipid densities
            for lipid_type in ['CHOL', 'DIPC', 'DPSM', 'DPG3']:
                density_col = f'Average_{lipid_type}_Density'
                if density_col in domain_data.columns:
                    avg_density = domain_data[density_col].mean()
                    print(f"    Average {lipid_type} density: {avg_density:.4f}/nm²")