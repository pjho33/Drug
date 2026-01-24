#!/usr/bin/env python3
"""
Tripod 1ns 분석 결과 시각화
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['axes.unicode_minus'] = False

RESULTS_DIR = "/home/pjho3tr/projects/Drug/2026-01-22_Glycogate_Tripod/results/md_tripod_1ns"
ARM_CSV = os.path.join(RESULTS_DIR, "tripod_arm_analysis.csv")
GLUCOSE_CSV = os.path.join(RESULTS_DIR, "glucose_distances.csv")
OUTPUT_PNG = os.path.join(RESULTS_DIR, "tripod_1ns_plots.png")

def main():
    arm_data = pd.read_csv(ARM_CSV)
    glucose_data = pd.read_csv(GLUCOSE_CSV)
    
    print(f"Loaded {len(arm_data)} frames")
    
    # Time series plot
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle('Tripod MD Simulation Analysis (1 ns)', fontsize=16, fontweight='bold')
    
    # 1. Individual arm lengths
    ax1 = axes[0, 0]
    ax1.plot(arm_data['time_ns'], arm_data['arm1_A'], 'b-', label='Arm 1', linewidth=2, alpha=0.8)
    ax1.plot(arm_data['time_ns'], arm_data['arm2_A'], 'r-', label='Arm 2', linewidth=2, alpha=0.8)
    ax1.plot(arm_data['time_ns'], arm_data['arm3_A'], 'g-', label='Arm 3', linewidth=2, alpha=0.8)
    ax1.axhline(y=40, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax1.axhline(y=50, color='gray', linestyle=':', linewidth=1, alpha=0.5)
    ax1.set_xlabel('Time (ns)', fontsize=12)
    ax1.set_ylabel('Arm Length (Å)', fontsize=12)
    ax1.set_title('Individual Arm Lengths', fontsize=13, fontweight='bold')
    ax1.legend(loc='best', fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # 2. Max arm length
    ax2 = axes[0, 1]
    ax2.plot(arm_data['time_ns'], arm_data['max_arm_A'], 'purple', linewidth=2, alpha=0.8)
    ax2.axhline(y=arm_data['max_arm_A'].mean(), color='darkviolet', linestyle='--', 
                linewidth=2, label=f'Mean: {arm_data["max_arm_A"].mean():.1f} Å')
    ax2.set_xlabel('Time (ns)', fontsize=12)
    ax2.set_ylabel('Max Arm Length (Å)', fontsize=12)
    ax2.set_title('Maximum Arm Length', fontsize=13, fontweight='bold')
    ax2.legend(loc='best', fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    # 3. Radius of gyration
    ax3 = axes[0, 2]
    ax3.plot(arm_data['time_ns'], arm_data['Rg_A'], 'orange', linewidth=2, alpha=0.8)
    ax3.axhline(y=arm_data['Rg_A'].mean(), color='darkorange', linestyle='--', 
                linewidth=2, label=f'Mean: {arm_data["Rg_A"].mean():.1f} Å')
    ax3.set_xlabel('Time (ns)', fontsize=12)
    ax3.set_ylabel('Rg (Å)', fontsize=12)
    ax3.set_title('Radius of Gyration', fontsize=13, fontweight='bold')
    ax3.legend(loc='best', fontsize=10)
    ax3.grid(True, alpha=0.3)
    
    # 4. Glucose distances
    ax4 = axes[1, 0]
    ax4.plot(glucose_data['time_ns'], glucose_data['glucose12_A'], 'b-', 
             label='G1-G2', linewidth=2, alpha=0.8)
    ax4.plot(glucose_data['time_ns'], glucose_data['glucose13_A'], 'r-', 
             label='G1-G3', linewidth=2, alpha=0.8)
    ax4.plot(glucose_data['time_ns'], glucose_data['glucose23_A'], 'g-', 
             label='G2-G3', linewidth=2, alpha=0.8)
    ax4.set_xlabel('Time (ns)', fontsize=12)
    ax4.set_ylabel('Distance (Å)', fontsize=12)
    ax4.set_title('Glucose-Glucose Distances', fontsize=13, fontweight='bold')
    ax4.legend(loc='best', fontsize=10)
    ax4.grid(True, alpha=0.3)
    
    # 5. Arm length distributions
    ax5 = axes[1, 1]
    ax5.hist(arm_data['arm1_A'], bins=10, alpha=0.6, label='Arm 1', color='blue', edgecolor='black')
    ax5.hist(arm_data['arm2_A'], bins=10, alpha=0.6, label='Arm 2', color='red', edgecolor='black')
    ax5.hist(arm_data['arm3_A'], bins=10, alpha=0.6, label='Arm 3', color='green', edgecolor='black')
    ax5.set_xlabel('Arm Length (Å)', fontsize=12)
    ax5.set_ylabel('Frequency', fontsize=12)
    ax5.set_title('Arm Length Distributions', fontsize=13, fontweight='bold')
    ax5.legend(loc='best', fontsize=10)
    ax5.grid(True, alpha=0.3, axis='y')
    
    # 6. Statistics table
    ax6 = axes[1, 2]
    ax6.axis('off')
    
    table_data = [
        ['Metric', 'Mean ± Std', 'Min - Max'],
        ['Arm 1 (Å)', 
         f'{arm_data["arm1_A"].mean():.1f} ± {arm_data["arm1_A"].std():.1f}',
         f'{arm_data["arm1_A"].min():.1f} - {arm_data["arm1_A"].max():.1f}'],
        ['Arm 2 (Å)', 
         f'{arm_data["arm2_A"].mean():.1f} ± {arm_data["arm2_A"].std():.1f}',
         f'{arm_data["arm2_A"].min():.1f} - {arm_data["arm2_A"].max():.1f}'],
        ['Arm 3 (Å)', 
         f'{arm_data["arm3_A"].mean():.1f} ± {arm_data["arm3_A"].std():.1f}',
         f'{arm_data["arm3_A"].min():.1f} - {arm_data["arm3_A"].max():.1f}'],
        ['Max Arm (Å)', 
         f'{arm_data["max_arm_A"].mean():.1f} ± {arm_data["max_arm_A"].std():.1f}',
         f'{arm_data["max_arm_A"].min():.1f} - {arm_data["max_arm_A"].max():.1f}'],
        ['Rg (Å)', 
         f'{arm_data["Rg_A"].mean():.1f} ± {arm_data["Rg_A"].std():.1f}',
         f'{arm_data["Rg_A"].min():.1f} - {arm_data["Rg_A"].max():.1f}'],
        ['G1-G2 (Å)', 
         f'{glucose_data["glucose12_A"].mean():.1f} ± {glucose_data["glucose12_A"].std():.1f}',
         f'{glucose_data["glucose12_A"].min():.1f} - {glucose_data["glucose12_A"].max():.1f}'],
        ['G1-G3 (Å)', 
         f'{glucose_data["glucose13_A"].mean():.1f} ± {glucose_data["glucose13_A"].std():.1f}',
         f'{glucose_data["glucose13_A"].min():.1f} - {glucose_data["glucose13_A"].max():.1f}'],
        ['G2-G3 (Å)', 
         f'{glucose_data["glucose23_A"].mean():.1f} ± {glucose_data["glucose23_A"].std():.1f}',
         f'{glucose_data["glucose23_A"].min():.1f} - {glucose_data["glucose23_A"].max():.1f}']
    ]
    
    table = ax6.table(cellText=table_data, cellLoc='left', loc='center',
                      colWidths=[0.3, 0.35, 0.35])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)
    
    # 헤더 스타일
    for i in range(3):
        table[(0, i)].set_facecolor('#4CAF50')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    ax6.set_title('Statistical Summary', fontsize=13, fontweight='bold', pad=20)
    
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(OUTPUT_PNG, dpi=300, bbox_inches='tight')
    print(f"\n✅ Plot saved: {OUTPUT_PNG}")

if __name__ == "__main__":
    main()
