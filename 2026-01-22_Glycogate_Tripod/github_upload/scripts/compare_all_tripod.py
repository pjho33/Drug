#!/usr/bin/env python3
"""
Tripod 1ns, 10ns, 100ns 비교 분석 및 시각화
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['axes.unicode_minus'] = False

BASE_DIR = "/home/pjho3tr/projects/Drug/2026-01-22_Glycogate_Tripod/results"

# 데이터 로드
data_1ns = pd.read_csv(f"{BASE_DIR}/md_tripod_1ns/tripod_arm_analysis.csv")
data_10ns = pd.read_csv(f"{BASE_DIR}/md_tripod_10ns/tripod_arm_analysis.csv")
data_100ns = pd.read_csv(f"{BASE_DIR}/md_tripod_100ns/tripod_arm_analysis.csv")

OUTPUT_PNG = f"{BASE_DIR}/tripod_comparison_all.png"

def main():
    print("Comparing 1ns, 10ns, 100ns simulations...")
    
    # 통계 계산
    stats = {
        '1ns': {
            'max_arm_mean': data_1ns['max_arm_A'].mean(),
            'max_arm_std': data_1ns['max_arm_A'].std(),
            'rg_mean': data_1ns['Rg_A'].mean(),
            'rg_std': data_1ns['Rg_A'].std(),
            'p40': (data_1ns['max_arm_A'] >= 40).mean() * 100,
            'p50': (data_1ns['max_arm_A'] >= 50).mean() * 100,
            'frames': len(data_1ns)
        },
        '10ns': {
            'max_arm_mean': data_10ns['max_arm_A'].mean(),
            'max_arm_std': data_10ns['max_arm_A'].std(),
            'rg_mean': data_10ns['Rg_A'].mean(),
            'rg_std': data_10ns['Rg_A'].std(),
            'p40': (data_10ns['max_arm_A'] >= 40).mean() * 100,
            'p50': (data_10ns['max_arm_A'] >= 50).mean() * 100,
            'frames': len(data_10ns)
        },
        '100ns': {
            'max_arm_mean': data_100ns['max_arm_A'].mean(),
            'max_arm_std': data_100ns['max_arm_A'].std(),
            'rg_mean': data_100ns['Rg_A'].mean(),
            'rg_std': data_100ns['Rg_A'].std(),
            'p40': (data_100ns['max_arm_A'] >= 40).mean() * 100,
            'p50': (data_100ns['max_arm_A'] >= 50).mean() * 100,
            'frames': len(data_100ns)
        }
    }
    
    # 시각화
    fig = plt.figure(figsize=(20, 12))
    gs = fig.add_gridspec(3, 4, hspace=0.3, wspace=0.3)
    
    fig.suptitle('Tripod MD Simulation: 1ns vs 10ns vs 100ns Comparison', 
                 fontsize=18, fontweight='bold', y=0.98)
    
    # 1. Max Arm Length Time Series (전체)
    ax1 = fig.add_subplot(gs[0, :2])
    ax1.plot(data_1ns['time_ns'], data_1ns['max_arm_A'], 'b-', 
             label='1ns', linewidth=2, alpha=0.8)
    ax1.axhline(y=stats['1ns']['max_arm_mean'], color='blue', 
                linestyle='--', linewidth=1.5, alpha=0.6)
    ax1.set_xlabel('Time (ns)', fontsize=12)
    ax1.set_ylabel('Max Arm Length (Å)', fontsize=12)
    ax1.set_title('1ns Simulation', fontsize=13, fontweight='bold')
    ax1.legend(loc='best', fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 1)
    
    ax2 = fig.add_subplot(gs[0, 2:])
    ax2.plot(data_10ns['time_ns'], data_10ns['max_arm_A'], 'r-', 
             label='10ns', linewidth=1.5, alpha=0.8)
    ax2.axhline(y=stats['10ns']['max_arm_mean'], color='red', 
                linestyle='--', linewidth=1.5, alpha=0.6)
    ax2.set_xlabel('Time (ns)', fontsize=12)
    ax2.set_ylabel('Max Arm Length (Å)', fontsize=12)
    ax2.set_title('10ns Simulation', fontsize=13, fontweight='bold')
    ax2.legend(loc='best', fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 10)
    
    # 2. 100ns Time Series
    ax3 = fig.add_subplot(gs[1, :])
    ax3.plot(data_100ns['time_ns'], data_100ns['max_arm_A'], 'g-', 
             label='100ns', linewidth=1, alpha=0.8)
    ax3.axhline(y=stats['100ns']['max_arm_mean'], color='green', 
                linestyle='--', linewidth=2, alpha=0.6, 
                label=f'Mean: {stats["100ns"]["max_arm_mean"]:.1f} Å')
    ax3.axhline(y=40, color='orange', linestyle='--', linewidth=1, alpha=0.5)
    ax3.axhline(y=50, color='red', linestyle='--', linewidth=1, alpha=0.5)
    ax3.set_xlabel('Time (ns)', fontsize=12)
    ax3.set_ylabel('Max Arm Length (Å)', fontsize=12)
    ax3.set_title('100ns Simulation', fontsize=13, fontweight='bold')
    ax3.legend(loc='best', fontsize=10)
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(0, 100)
    
    # 3. Distribution Comparison
    ax4 = fig.add_subplot(gs[2, 0])
    ax4.hist(data_1ns['max_arm_A'], bins=10, alpha=0.7, label='1ns', 
             color='blue', edgecolor='black')
    ax4.axvline(x=40, color='orange', linestyle='--', linewidth=2, alpha=0.7)
    ax4.axvline(x=50, color='red', linestyle='--', linewidth=2, alpha=0.7)
    ax4.set_xlabel('Max Arm Length (Å)', fontsize=11)
    ax4.set_ylabel('Frequency', fontsize=11)
    ax4.set_title('1ns Distribution', fontsize=12, fontweight='bold')
    ax4.legend(loc='best', fontsize=9)
    ax4.grid(True, alpha=0.3, axis='y')
    
    ax5 = fig.add_subplot(gs[2, 1])
    ax5.hist(data_10ns['max_arm_A'], bins=20, alpha=0.7, label='10ns', 
             color='red', edgecolor='black')
    ax5.axvline(x=40, color='orange', linestyle='--', linewidth=2, alpha=0.7)
    ax5.axvline(x=50, color='red', linestyle='--', linewidth=2, alpha=0.7)
    ax5.set_xlabel('Max Arm Length (Å)', fontsize=11)
    ax5.set_ylabel('Frequency', fontsize=11)
    ax5.set_title('10ns Distribution', fontsize=12, fontweight='bold')
    ax5.legend(loc='best', fontsize=9)
    ax5.grid(True, alpha=0.3, axis='y')
    
    ax6 = fig.add_subplot(gs[2, 2])
    ax6.hist(data_100ns['max_arm_A'], bins=30, alpha=0.7, label='100ns', 
             color='green', edgecolor='black')
    ax6.axvline(x=40, color='orange', linestyle='--', linewidth=2, alpha=0.7)
    ax6.axvline(x=50, color='red', linestyle='--', linewidth=2, alpha=0.7)
    ax6.set_xlabel('Max Arm Length (Å)', fontsize=11)
    ax6.set_ylabel('Frequency', fontsize=11)
    ax6.set_title('100ns Distribution', fontsize=12, fontweight='bold')
    ax6.legend(loc='best', fontsize=9)
    ax6.grid(True, alpha=0.3, axis='y')
    
    # 4. Comparison Table
    ax7 = fig.add_subplot(gs[2, 3])
    ax7.axis('off')
    
    table_data = [
        ['Metric', '1ns', '10ns', '100ns'],
        ['Frames', 
         f'{stats["1ns"]["frames"]}',
         f'{stats["10ns"]["frames"]}',
         f'{stats["100ns"]["frames"]}'],
        ['Max Arm (Å)', 
         f'{stats["1ns"]["max_arm_mean"]:.1f}±{stats["1ns"]["max_arm_std"]:.1f}',
         f'{stats["10ns"]["max_arm_mean"]:.1f}±{stats["10ns"]["max_arm_std"]:.1f}',
         f'{stats["100ns"]["max_arm_mean"]:.1f}±{stats["100ns"]["max_arm_std"]:.1f}'],
        ['Rg (Å)', 
         f'{stats["1ns"]["rg_mean"]:.1f}±{stats["1ns"]["rg_std"]:.1f}',
         f'{stats["10ns"]["rg_mean"]:.1f}±{stats["10ns"]["rg_std"]:.1f}',
         f'{stats["100ns"]["rg_mean"]:.1f}±{stats["100ns"]["rg_std"]:.1f}'],
        ['P(≥40Å) %', 
         f'{stats["1ns"]["p40"]:.1f}',
         f'{stats["10ns"]["p40"]:.1f}',
         f'{stats["100ns"]["p40"]:.1f}'],
        ['P(≥50Å) %', 
         f'{stats["1ns"]["p50"]:.1f}',
         f'{stats["10ns"]["p50"]:.1f}',
         f'{stats["100ns"]["p50"]:.1f}']
    ]
    
    table = ax7.table(cellText=table_data, cellLoc='center', loc='center',
                      colWidths=[0.25, 0.25, 0.25, 0.25])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2.5)
    
    # 헤더 스타일
    for i in range(4):
        table[(0, i)].set_facecolor('#2196F3')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    ax7.set_title('Statistical Comparison', fontsize=12, fontweight='bold', pad=20)
    
    plt.savefig(OUTPUT_PNG, dpi=300, bbox_inches='tight')
    print(f"\n✅ Comparison plot saved: {OUTPUT_PNG}")
    
    # 통계 출력
    print("\n" + "=" * 80)
    print("COMPARISON SUMMARY")
    print("=" * 80)
    for sim in ['1ns', '10ns', '100ns']:
        print(f"\n{sim.upper()}:")
        print(f"  Frames: {stats[sim]['frames']}")
        print(f"  Max Arm: {stats[sim]['max_arm_mean']:.2f} ± {stats[sim]['max_arm_std']:.2f} Å")
        print(f"  Rg: {stats[sim]['rg_mean']:.2f} ± {stats[sim]['rg_std']:.2f} Å")
        print(f"  P(≥40Å): {stats[sim]['p40']:.1f}%")
        print(f"  P(≥50Å): {stats[sim]['p50']:.1f}%")
    print("=" * 80)

if __name__ == "__main__":
    main()
