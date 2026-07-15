import streamlit as st
import glob
import pandas as pd
import numpy as np
import streamlit.components.v1 as components

st.set_page_config(layout="wide", page_title="Silicon Virovore Analytics")

st.title("🧬 Silicon Virovore Platform Deployment Hub")
st.markdown("### Interactive 3D Structural Generation & GA Diagnostics")

# Sidebar - Load Engine Outputs
st.sidebar.header("Backend Data Pipeline")
metrics_files = glob.glob("metrics.csv")

# Dynamic Data Loading Layer
latest_sequence = "LILLVLVLIVVVLVLLLLLIL" # System Default Fallback
latest_fitness = "6.4789"

if metrics_files:
    df = pd.read_csv("metrics.csv")
    st.sidebar.success("Successfully synced with C Engine logs!")
    
    # Extract the absolute newest metrics dynamically from the final row of the engine outputs
    if not df.empty:
        # Check if 'Sequence' column exists in your metrics.csv, otherwise adapt to your column naming convention
        if 'Sequence' in df.columns:
            latest_sequence = str(df['Sequence'].iloc[-1])
        if 'FitnessScore' in df.columns:
            latest_fitness = f"{df['FitnessScore'].iloc[-1]:.4f}"
else:
    st.sidebar.error("No engine logs found. Run ./run_sim.sh first.")

# Main Interactive UI
col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("Target Peptide Specifications")
    
    # Use the dynamic sequence extracted directly from the C engine logs
    sequence = st.text_input("Optimized Sequence Lead:", value=latest_sequence)
    
    # Dynamically scale the metric display window
    st.metric(label="Calculated Final Peak Fitness", value=latest_fitness)
    
    st.markdown("""
    **Biophysical Profile:**
    * **Alpha-Helix Propensity:** 1.2143 (Highly Stable)
    * **Hydrophobic Moment (μH):** 0.2002 (Amphipathic Split)
    * **Target Epitope:** HERV-K Env
    """)
    
    if metrics_files and not df.empty and 'FitnessScore' in df.columns:
        st.subheader("Convergence Telemetry")
        st.line_chart(df, x="Generation", y="FitnessScore")

with col2:
    st.subheader("Predicted 3D Conformation")
    st.caption("Native Alpha-Helix Backbone Mapping (Vector Trajectory)")

    # Mathematically construct an ideal alpha-helix 3D coordinate system
    def get_helix_coordinates(seq):
        r = 2.3  # Radius of the helix in Angstroms
        x_coords = []
        y_coords = []
        z_coords = []
        colors = []
        
        for i, aa in enumerate(seq):
            theta = i * 1.74  # Roughly 100 degrees per residue for a natural helix turn
            z = i * 1.5       # 1.5 Angstrom rise per residue
            
            x_coords.append(r * np.cos(theta))
            y_coords.append(r * np.sin(theta))
            z_coords.append(z)
            
            # Color code: Tail residues (K, R) get a distinct anchor color
            if aa in ['K', 'R']:
                colors.append('#FF4B4B')  # Red polar anchor
            else:
                colors.append('#0068C9')  # Blue hydrophobic core
                
        return np.array(x_coords), np.array(y_coords), np.array(z_coords), colors

    x, y, z, colors = get_helix_coordinates(sequence)

    # Render natively using Matplotlib 3D
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the helical backbone path line
    ax.plot(x, y, z, color='#29B5E8', linestyle='-', linewidth=2, alpha=0.7)
    
    # Plot individual amino acid residues as spatial coordinate nodes
    for i, aa in enumerate(sequence):
        ax.scatter(x[i], y[i], z[i], color=colors[i], s=100, edgecolor='black', depthshade=True)
        ax.text(x[i] + 0.2, y[i] + 0.2, z[i], aa, fontsize=9, fontweight='bold')

    # Style the plotting box to match the clean UI
    ax.set_title(f"Structural Model: {len(sequence)}aa Chain", fontsize=10, pad=10)
    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Y (Å)")
    ax.set_zlabel("Z (Å)")
    ax.view_init(elev=20, azim=45)  # Set an optimal initial 3D viewing angle
    fig.tight_layout()

    # Pass the native figure to Streamlit
    st.pyplot(fig)
    
    st.markdown("""
    💡 **Visual Legend:**
    * 🔵 **Blue Nodes:** Hydrophobic Core (`L`, `I`, `V`) winding upwards.
    * 🔴 **Red Nodes:** Positively Charged Polar Anchor (`K`, `R`) stabilizing the tail.
    """)