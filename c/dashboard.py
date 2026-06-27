# dashboard.py
import os
import pandas as pd
import matplotlib.pyplot as plt

def generate_mock_metrics():
    """Generates a metrics file if your C loop hasn't written one yet."""
    generations = list(range(1, 101))
    # Simulating a genetic algorithm optimizing over 100 generations
    fitness = [-3.5 + (i ** 0.35) * 0.8 for i in generations]
    charge_penalty = [max(0.0, 4.0 - (i * 0.1)) for i in generations]
    
    df = pd.DataFrame({
        'Generation': generations,
        'FitnessScore': fitness,
        'ChargePenalty': charge_penalty
    })
    df.to_csv('metrics.csv', index=False)
    print("[PYTHON] Mock metrics.csv created for dashboard preview.")

def plot_dashboard():
    if not os.path.exists('metrics.csv'):
        generate_mock_metrics()
        
    # Read the data from the evolutionary run
    df = pd.read_csv('metrics.csv')
    
    # Set style for a clean presentation look
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Graph 1: Fitness Score Optimization Convergence
    ax1.plot(df['Generation'], df['FitnessScore'], color='#007AFF', linewidth=2.5, label='Target Optimization')
    ax1.set_title('Genomic Fitness Convergence (1,000 Variants)', fontsize=12, fontweight='bold', pad=10)
    ax1.set_xlabel('Generation Count', fontsize=10)
    ax1.set_ylabel('Raw Fitness Value', fontsize=10)
    ax1.legend(loc='lower right')
    
    # Graph 2: Negative Design Penalty Decay
    ax2.plot(df['Generation'], df['ChargePenalty'], color='#FF3B30', linewidth=2.5, label='Biophysical Penalties')
    ax2.set_title('Decoy Cross-Reactivity Penalty Decay', fontsize=12, fontweight='bold', pad=10)
    ax2.set_xlabel('Generation Count', fontsize=10)
    ax2.set_ylabel('Penalty Score', fontsize=10)
    ax2.legend(loc='upper right')
    
    plt.tight_layout()
    # Save directly to an image for your slides
    plt.savefig('virovore_dashboard.png', dpi=300)
    print("[PYTHON] Dashboard rendering successful! Saved as 'virovore_dashboard.png'.")

if __name__ == "__main__":
    plot_dashboard()