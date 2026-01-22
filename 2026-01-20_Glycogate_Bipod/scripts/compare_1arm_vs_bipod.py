#!/usr/bin/env python3
"""
1-arm (Glycogate) vs Bipod 비교 분석
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['axes.unicode_minus'] = False

# 파일 경로
# 1-arm (Glycogate) - 200ns
ARM1_LENGTH = "/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/results/md_200ns_rep1/lig04_one_arm_length_last200ns.csv"
ARM1_RG = "/home/pjho3tr/projects/Drug/2026-01-18_Glycogate/results/md_200ns_rep1/lig02_rg_last200ns.csv"

# Bipod - 10ns
BIPOD_10NS_ARM = "/home/pjho3tr/projects/Drug/2026-01-20_Glycogate_Bipod/results/md_bipod_10ns/bipod_arm_analysis.csv"
BIPOD_10NS_GLUCOSE = "/home/pjho3tr/projects/Drug/2026-01-20_Glycogate_Bipod/results/md_bipod_10ns/glucose_distance.csv"

# Bipod - 100ns
BIPOD_100NS_ARM = "/home/pjho3tr/projects/Drug/2026-01-20_Glycogate_Bipod/results/md_bipod_100ns/bipod_arm_analysis.csv"
BIPOD_100NS_GLUCOSE = "/home/pjho3tr/projects/Drug/2026-01-20_Glycogate_Bipod/results/md_bipod_100ns/glucose_distance.csv"

# 출력
OUTDIR = "/home/pjho3tr/projects/Drug/2026-01-20_Glycogate_Bipod/results/comparison"
os.makedirs(OUTDIR, exist_ok=True)

def main():
    # 데이터 로드
    print("Loading data...")
    
    # 1-arm
    arm1_data = pd.read_csv(ARM1_LENGTH)
    arm1_rg = pd.read_csv(ARM1_RG)
    
    # Bipod 10ns
    bipod_10ns_arm = pd.read_csv(BIPOD_10NS_ARM)
    bipod_10ns_glucose = pd.read_csv(BIPOD_10NS_GLUCOSE)
    
    # Bipod 100ns
    bipod_100ns_arm = pd.read_csv(BIPOD_100NS_ARM)
    bipod_100ns_glucose = pd.read_csv(BIPOD_100NS_GLUCOSE)
    
    print(f"1-arm: {len(arm1_data)} frames")
    print(f"Bipod 10ns: {len(bipod_10ns_arm)} frames")
    print(f"Bipod 100ns: {len(bipod_100ns_arm)} frames")
    
    # 통계 계산
    stats = {
        '1-arm (200ns)': {
            'arm_length': arm1_data['arm_A'].values,
            'rg': arm1_rg['Rg_A'].values,
            'n_frames': len(arm1_data)
        },
        'Bipod 10ns': {
            'arm1': bipod_10ns_arm['arm1_A'].values,
            'arm2': bipod_10ns_arm['arm2_A'].values,
            'max_arm': bipod_10ns_arm['max_arm_A'].values,
            'glucose_dist': bipod_10ns_glucose['glucose_distance_A'].values,
            'rg': bipod_10ns_arm['Rg_A'].values,
            'n_frames': len(bipod_10ns_arm)
        },
        'Bipod 100ns': {
            'arm1': bipod_100ns_arm['arm1_A'].values,
            'arm2': bipod_100ns_arm['arm2_A'].values,
            'max_arm': bipod_100ns_arm['max_arm_A'].values,
            'glucose_dist': bipod_100ns_glucose['glucose_distance_A'].values,
            'rg': bipod_100ns_arm['Rg_A'].values,
            'n_frames': len(bipod_100ns_arm)
        }
    }
    
    # 통계 요약 출력
    print("\n" + "=" * 80)
    print("COMPARISON STATISTICS")
    print("=" * 80)
    
    print("\n1-ARM (Glycogate, 200ns):")
    print(f"  Arm length: {stats['1-arm (200ns)']['arm_length'].mean():.2f} ± {stats['1-arm (200ns)']['arm_length'].std():.2f} Å")
    print(f"  Range: {stats['1-arm (200ns)']['arm_length'].min():.2f} - {stats['1-arm (200ns)']['arm_length'].max():.2f} Å")
    print(f"  Rg: {stats['1-arm (200ns)']['rg'].mean():.2f} ± {stats['1-arm (200ns)']['rg'].std():.2f} Å")
    print(f"  P(>= 40Å): {(stats['1-arm (200ns)']['arm_length'] >= 40).mean() * 100:.1f}%")
    print(f"  P(>= 50Å): {(stats['1-arm (200ns)']['arm_length'] >= 50).mean() * 100:.1f}%")
    
    print("\nBIPOD 10ns:")
    print(f"  Arm 1: {stats['Bipod 10ns']['arm1'].mean():.2f} ± {stats['Bipod 10ns']['arm1'].std():.2f} Å")
    print(f"  Arm 2: {stats['Bipod 10ns']['arm2'].mean():.2f} ± {stats['Bipod 10ns']['arm2'].std():.2f} Å")
    print(f"  Max arm: {stats['Bipod 10ns']['max_arm'].mean():.2f} ± {stats['Bipod 10ns']['max_arm'].std():.2f} Å")
    print(f"  Glucose distance: {stats['Bipod 10ns']['glucose_dist'].mean():.2f} ± {stats['Bipod 10ns']['glucose_dist'].std():.2f} Å")
    print(f"  Rg: {stats['Bipod 10ns']['rg'].mean():.2f} ± {stats['Bipod 10ns']['rg'].std():.2f} Å")
    print(f"  P(>= 40Å): {(stats['Bipod 10ns']['max_arm'] >= 40).mean() * 100:.1f}%")
    print(f"  P(>= 50Å): {(stats['Bipod 10ns']['max_arm'] >= 50).mean() * 100:.1f}%")
    
    print("\nBIPOD 100ns:")
    print(f"  Arm 1: {stats['Bipod 100ns']['arm1'].mean():.2f} ± {stats['Bipod 100ns']['arm1'].std():.2f} Å")
    print(f"  Arm 2: {stats['Bipod 100ns']['arm2'].mean():.2f} ± {stats['Bipod 100ns']['arm2'].std():.2f} Å")
    print(f"  Max arm: {stats['Bipod 100ns']['max_arm'].mean():.2f} ± {stats['Bipod 100ns']['max_arm'].std():.2f} Å")
    print(f"  Glucose distance: {stats['Bipod 100ns']['glucose_dist'].mean():.2f} ± {stats['Bipod 100ns']['glucose_dist'].std():.2f} Å")
    print(f"  Rg: {stats['Bipod 100ns']['rg'].mean():.2f} ± {stats['Bipod 100ns']['rg'].std():.2f} Å")
    print(f"  P(>= 40Å): {(stats['Bipod 100ns']['max_arm'] >= 40).mean() * 100:.1f}%")
    print(f"  P(>= 50Å): {(stats['Bipod 100ns']['max_arm'] >= 50).mean() * 100:.1f}%")
    
    print("=" * 80)
    
    # 비교 플롯 생성
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle('1-arm vs Bipod Comparison', fontsize=16, fontweight='bold')
    
    # 1. Arm length 분포 비교
    ax1 = axes[0, 0]
    ax1.hist(stats['1-arm (200ns)']['arm_length'], bins=50, alpha=0.6, label='1-arm (200ns)', 
             color='blue', edgecolor='black', density=True)
    ax1.hist(stats['Bipod 10ns']['max_arm'], bins=50, alpha=0.6, label='Bipod 10ns', 
             color='green', edgecolor='black', density=True)
    ax1.hist(stats['Bipod 100ns']['max_arm'], bins=50, alpha=0.6, label='Bipod 100ns', 
             color='red', edgecolor='black', density=True)
    ax1.set_xlabel('Arm Length (Å)', fontsize=11)
    ax1.set_ylabel('Density', fontsize=11)
    ax1.set_title('Arm Length Distribution', fontsize=12, fontweight='bold')
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    
    # 2. Rg 분포 비교
    ax2 = axes[0, 1]
    ax2.hist(stats['1-arm (200ns)']['rg'], bins=50, alpha=0.6, label='1-arm (200ns)', 
             color='blue', edgecolor='black', density=True)
    ax2.hist(stats['Bipod 10ns']['rg'], bins=50, alpha=0.6, label='Bipod 10ns', 
             color='green', edgecolor='black', density=True)
    ax2.hist(stats['Bipod 100ns']['rg'], bins=50, alpha=0.6, label='Bipod 100ns', 
             color='red', edgecolor='black', density=True)
    ax2.set_xlabel('Radius of Gyration (Å)', fontsize=11)
    ax2.set_ylabel('Density', fontsize=11)
    ax2.set_title('Rg Distribution', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    
    # 3. Box plot 비교 - Arm length
    ax3 = axes[0, 2]
    data_to_plot = [
        stats['1-arm (200ns)']['arm_length'],
        stats['Bipod 10ns']['max_arm'],
        stats['Bipod 100ns']['max_arm']
    ]
    bp = ax3.boxplot(data_to_plot, labels=['1-arm\n(200ns)', 'Bipod\n10ns', 'Bipod\n100ns'],
                     patch_artist=True)
    colors = ['lightblue', 'lightgreen', 'lightcoral']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    ax3.set_ylabel('Arm Length (Å)', fontsize=11)
    ax3.set_title('Arm Length Comparison', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3, axis='y')
    
    # 4. Box plot 비교 - Rg
    ax4 = axes[1, 0]
    data_to_plot_rg = [
        stats['1-arm (200ns)']['rg'],
        stats['Bipod 10ns']['rg'],
        stats['Bipod 100ns']['rg']
    ]
    bp2 = ax4.boxplot(data_to_plot_rg, labels=['1-arm\n(200ns)', 'Bipod\n10ns', 'Bipod\n100ns'],
                      patch_artist=True)
    for patch, color in zip(bp2['boxes'], colors):
        patch.set_facecolor(color)
    ax4.set_ylabel('Radius of Gyration (Å)', fontsize=11)
    ax4.set_title('Rg Comparison', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3, axis='y')
    
    # 5. 통계 비교 테이블
    ax5 = axes[1, 1]
    ax5.axis('off')
    
    table_data = [
        ['Metric', '1-arm\n(200ns)', 'Bipod\n10ns', 'Bipod\n100ns'],
        ['Arm (Å)', 
         f'{stats["1-arm (200ns)"]["arm_length"].mean():.1f}±{stats["1-arm (200ns)"]["arm_length"].std():.1f}',
         f'{stats["Bipod 10ns"]["max_arm"].mean():.1f}±{stats["Bipod 10ns"]["max_arm"].std():.1f}',
         f'{stats["Bipod 100ns"]["max_arm"].mean():.1f}±{stats["Bipod 100ns"]["max_arm"].std():.1f}'],
        ['Rg (Å)',
         f'{stats["1-arm (200ns)"]["rg"].mean():.1f}±{stats["1-arm (200ns)"]["rg"].std():.1f}',
         f'{stats["Bipod 10ns"]["rg"].mean():.1f}±{stats["Bipod 10ns"]["rg"].std():.1f}',
         f'{stats["Bipod 100ns"]["rg"].mean():.1f}±{stats["Bipod 100ns"]["rg"].std():.1f}'],
        ['P(≥40Å) %',
         f'{(stats["1-arm (200ns)"]["arm_length"] >= 40).mean() * 100:.1f}',
         f'{(stats["Bipod 10ns"]["max_arm"] >= 40).mean() * 100:.1f}',
         f'{(stats["Bipod 100ns"]["max_arm"] >= 40).mean() * 100:.1f}'],
        ['P(≥50Å) %',
         f'{(stats["1-arm (200ns)"]["arm_length"] >= 50).mean() * 100:.1f}',
         f'{(stats["Bipod 10ns"]["max_arm"] >= 50).mean() * 100:.1f}',
         f'{(stats["Bipod 100ns"]["max_arm"] >= 50).mean() * 100:.1f}']
    ]
    
    table = ax5.table(cellText=table_data, cellLoc='center', loc='center',
                      colWidths=[0.25, 0.25, 0.25, 0.25])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)
    
    # 헤더 스타일
    for i in range(4):
        table[(0, i)].set_facecolor('#4CAF50')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    ax5.set_title('Statistical Comparison', fontsize=12, fontweight='bold', pad=20)
    
    # 6. Bipod glucose distance (Bipod만 해당)
    ax6 = axes[1, 2]
    ax6.hist(stats['Bipod 10ns']['glucose_dist'], bins=30, alpha=0.6, 
             label='Bipod 10ns', color='green', edgecolor='black', density=True)
    ax6.hist(stats['Bipod 100ns']['glucose_dist'], bins=30, alpha=0.6, 
             label='Bipod 100ns', color='red', edgecolor='black', density=True)
    ax6.set_xlabel('Glucose-Glucose Distance (Å)', fontsize=11)
    ax6.set_ylabel('Density', fontsize=11)
    ax6.set_title('Bipod: Glucose Distance', fontsize=12, fontweight='bold')
    ax6.legend(fontsize=9)
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    
    output_png = os.path.join(OUTDIR, "1arm_vs_bipod_comparison.png")
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    print(f"\n✅ Comparison plot saved: {output_png}")
    
    # CSV 요약 저장
    summary_csv = os.path.join(OUTDIR, "comparison_summary.csv")
    summary_data = {
        'System': ['1-arm (200ns)', 'Bipod 10ns', 'Bipod 100ns'],
        'Arm_mean': [
            stats['1-arm (200ns)']['arm_length'].mean(),
            stats['Bipod 10ns']['max_arm'].mean(),
            stats['Bipod 100ns']['max_arm'].mean()
        ],
        'Arm_std': [
            stats['1-arm (200ns)']['arm_length'].std(),
            stats['Bipod 10ns']['max_arm'].std(),
            stats['Bipod 100ns']['max_arm'].std()
        ],
        'Rg_mean': [
            stats['1-arm (200ns)']['rg'].mean(),
            stats['Bipod 10ns']['rg'].mean(),
            stats['Bipod 100ns']['rg'].mean()
        ],
        'Rg_std': [
            stats['1-arm (200ns)']['rg'].std(),
            stats['Bipod 10ns']['rg'].std(),
            stats['Bipod 100ns']['rg'].std()
        ],
        'P_40A': [
            (stats['1-arm (200ns)']['arm_length'] >= 40).mean() * 100,
            (stats['Bipod 10ns']['max_arm'] >= 40).mean() * 100,
            (stats['Bipod 100ns']['max_arm'] >= 40).mean() * 100
        ],
        'P_50A': [
            (stats['1-arm (200ns)']['arm_length'] >= 50).mean() * 100,
            (stats['Bipod 10ns']['max_arm'] >= 50).mean() * 100,
            (stats['Bipod 100ns']['max_arm'] >= 50).mean() * 100
        ]
    }
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(summary_csv, index=False)
    print(f"✅ Summary CSV saved: {summary_csv}")

if __name__ == "__main__":
    main()
