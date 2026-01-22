#!/usr/bin/env python3
"""
Bipod 분석 결과 시각화
- Arm 길이 time series
- Glucose 간 거리 time series
- Radius of gyration
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# 한글 폰트 설정
mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['axes.unicode_minus'] = False

# 파일 경로
RESULTS_DIR = "/home/pjho3tr/projects/Drug/2026-01-20_Glycogate_Bipod/results/md_bipod_10ns"
ARM_CSV = os.path.join(RESULTS_DIR, "bipod_arm_analysis.csv")
GLUCOSE_CSV = os.path.join(RESULTS_DIR, "glucose_distance.csv")
OUTPUT_PNG = os.path.join(RESULTS_DIR, "bipod_analysis_plots.png")

def main():
    # 데이터 로드
    arm_data = pd.read_csv(ARM_CSV)
    glucose_data = pd.read_csv(GLUCOSE_CSV)
    
    print(f"Loaded {len(arm_data)} frames from arm analysis")
    print(f"Loaded {len(glucose_data)} frames from glucose distance")
    
    # Figure 생성 (2x2 subplot)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Bipod MD Simulation Analysis (10 ns)', fontsize=16, fontweight='bold')
    
    # 1. Arm 길이 time series
    ax1 = axes[0, 0]
    ax1.plot(arm_data['time_ns'], arm_data['arm1_A'], 'b-', label='Arm 1', linewidth=1.5, alpha=0.8)
    ax1.plot(arm_data['time_ns'], arm_data['arm2_A'], 'r-', label='Arm 2', linewidth=1.5, alpha=0.8)
    ax1.axhline(y=40, color='gray', linestyle='--', linewidth=1, alpha=0.5, label='40 A (4 nm)')
    ax1.axhline(y=50, color='gray', linestyle=':', linewidth=1, alpha=0.5, label='50 A (5 nm)')
    ax1.set_xlabel('Time (ns)', fontsize=12)
    ax1.set_ylabel('Arm Length (A)', fontsize=12)
    ax1.set_title('Individual Arm Lengths', fontsize=13, fontweight='bold')
    ax1.legend(loc='best', fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # 2. Max arm 길이
    ax2 = axes[0, 1]
    ax2.plot(arm_data['time_ns'], arm_data['max_arm_A'], 'g-', linewidth=2, alpha=0.8)
    ax2.fill_between(arm_data['time_ns'], arm_data['max_arm_A'], alpha=0.3, color='green')
    ax2.axhline(y=arm_data['max_arm_A'].mean(), color='darkgreen', linestyle='--', 
                linewidth=2, label=f'Mean: {arm_data["max_arm_A"].mean():.1f} A')
    ax2.set_xlabel('Time (ns)', fontsize=12)
    ax2.set_ylabel('Max Arm Length (A)', fontsize=12)
    ax2.set_title('Maximum Arm Length', fontsize=13, fontweight='bold')
    ax2.legend(loc='best', fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    # 3. Glucose 간 거리
    ax3 = axes[1, 0]
    ax3.plot(glucose_data['time_ns'], glucose_data['glucose_distance_A'], 
             'purple', linewidth=2, alpha=0.8)
    ax3.fill_between(glucose_data['time_ns'], glucose_data['glucose_distance_A'], 
                      alpha=0.3, color='purple')
    ax3.axhline(y=glucose_data['glucose_distance_A'].mean(), color='darkviolet', 
                linestyle='--', linewidth=2, 
                label=f'Mean: {glucose_data["glucose_distance_A"].mean():.1f} A')
    ax3.set_xlabel('Time (ns)', fontsize=12)
    ax3.set_ylabel('Distance (A)', fontsize=12)
    ax3.set_title('Distance between Two L-glucose Groups', fontsize=13, fontweight='bold')
    ax3.legend(loc='best', fontsize=10)
    ax3.grid(True, alpha=0.3)
    
    # 4. Radius of Gyration
    ax4 = axes[1, 1]
    ax4.plot(arm_data['time_ns'], arm_data['Rg_A'], 'orange', linewidth=2, alpha=0.8)
    ax4.fill_between(arm_data['time_ns'], arm_data['Rg_A'], alpha=0.3, color='orange')
    ax4.axhline(y=arm_data['Rg_A'].mean(), color='darkorange', linestyle='--', 
                linewidth=2, label=f'Mean: {arm_data["Rg_A"].mean():.1f} A')
    ax4.set_xlabel('Time (ns)', fontsize=12)
    ax4.set_ylabel('Rg (A)', fontsize=12)
    ax4.set_title('Radius of Gyration', fontsize=13, fontweight='bold')
    ax4.legend(loc='best', fontsize=10)
    ax4.grid(True, alpha=0.3)
    
    # 레이아웃 조정
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    
    # 저장
    plt.savefig(OUTPUT_PNG, dpi=300, bbox_inches='tight')
    print(f"\n✅ Plot saved: {OUTPUT_PNG}")
    
    # 통계 출력
    print("\n" + "=" * 80)
    print("STATISTICS SUMMARY")
    print("=" * 80)
    print("\nArm 1:")
    print(f"  Mean ± Std: {arm_data['arm1_A'].mean():.2f} ± {arm_data['arm1_A'].std():.2f} A")
    print(f"  Min - Max: {arm_data['arm1_A'].min():.2f} - {arm_data['arm1_A'].max():.2f} A")
    
    print("\nArm 2:")
    print(f"  Mean ± Std: {arm_data['arm2_A'].mean():.2f} ± {arm_data['arm2_A'].std():.2f} A")
    print(f"  Min - Max: {arm_data['arm2_A'].min():.2f} - {arm_data['arm2_A'].max():.2f} A")
    
    print("\nMax Arm:")
    print(f"  Mean ± Std: {arm_data['max_arm_A'].mean():.2f} ± {arm_data['max_arm_A'].std():.2f} A")
    print(f"  Min - Max: {arm_data['max_arm_A'].min():.2f} - {arm_data['max_arm_A'].max():.2f} A")
    
    print("\nGlucose Distance:")
    print(f"  Mean ± Std: {glucose_data['glucose_distance_A'].mean():.2f} ± {glucose_data['glucose_distance_A'].std():.2f} A")
    print(f"  Min - Max: {glucose_data['glucose_distance_A'].min():.2f} - {glucose_data['glucose_distance_A'].max():.2f} A")
    
    print("\nRadius of Gyration:")
    print(f"  Mean ± Std: {arm_data['Rg_A'].mean():.2f} ± {arm_data['Rg_A'].std():.2f} A")
    print(f"  Min - Max: {arm_data['Rg_A'].min():.2f} - {arm_data['Rg_A'].max():.2f} A")
    print("=" * 80)
    
    # 추가 플롯: 분포 히스토그램
    fig2, axes2 = plt.subplots(2, 2, figsize=(14, 10))
    fig2.suptitle('Bipod Analysis - Distributions', fontsize=16, fontweight='bold')
    
    # Arm 1 분포
    ax1 = axes2[0, 0]
    ax1.hist(arm_data['arm1_A'], bins=30, color='blue', alpha=0.7, edgecolor='black')
    ax1.axvline(arm_data['arm1_A'].mean(), color='red', linestyle='--', linewidth=2, 
                label=f'Mean: {arm_data["arm1_A"].mean():.1f} A')
    ax1.set_xlabel('Arm 1 Length (A)', fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    ax1.set_title('Arm 1 Length Distribution', fontsize=13, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Arm 2 분포
    ax2 = axes2[0, 1]
    ax2.hist(arm_data['arm2_A'], bins=30, color='red', alpha=0.7, edgecolor='black')
    ax2.axvline(arm_data['arm2_A'].mean(), color='darkred', linestyle='--', linewidth=2,
                label=f'Mean: {arm_data["arm2_A"].mean():.1f} A')
    ax2.set_xlabel('Arm 2 Length (A)', fontsize=12)
    ax2.set_ylabel('Frequency', fontsize=12)
    ax2.set_title('Arm 2 Length Distribution', fontsize=13, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Glucose 거리 분포
    ax3 = axes2[1, 0]
    ax3.hist(glucose_data['glucose_distance_A'], bins=30, color='purple', alpha=0.7, edgecolor='black')
    ax3.axvline(glucose_data['glucose_distance_A'].mean(), color='darkviolet', linestyle='--', 
                linewidth=2, label=f'Mean: {glucose_data["glucose_distance_A"].mean():.1f} A')
    ax3.set_xlabel('Glucose Distance (A)', fontsize=12)
    ax3.set_ylabel('Frequency', fontsize=12)
    ax3.set_title('Glucose Distance Distribution', fontsize=13, fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Rg 분포
    ax4 = axes2[1, 1]
    ax4.hist(arm_data['Rg_A'], bins=30, color='orange', alpha=0.7, edgecolor='black')
    ax4.axvline(arm_data['Rg_A'].mean(), color='darkorange', linestyle='--', linewidth=2,
                label=f'Mean: {arm_data["Rg_A"].mean():.1f} A')
    ax4.set_xlabel('Radius of Gyration (A)', fontsize=12)
    ax4.set_ylabel('Frequency', fontsize=12)
    ax4.set_title('Rg Distribution', fontsize=13, fontweight='bold')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    
    output_hist = os.path.join(RESULTS_DIR, "bipod_analysis_distributions.png")
    plt.savefig(output_hist, dpi=300, bbox_inches='tight')
    print(f"✅ Distribution plot saved: {output_hist}")

if __name__ == "__main__":
    main()
