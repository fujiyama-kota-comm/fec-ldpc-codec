import pandas as pd
import matplotlib.pyplot as plt
import os

# =============================================================================
#  Plot BER graphs for LDPC (SPA Decoder)
# =============================================================================

# 出力フォルダ
os.makedirs("images", exist_ok=True)

# =============================================================================
#  Matplotlib settings (論文品質)
# =============================================================================
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.size"] = 14

# =============================================================================
#  Load CSV data
# =============================================================================
df = pd.read_csv("results/ldpc_ber_data.csv")

EbN0 = df["EbN0_dB"]
ber_info = df["BER_info"]
ber_bpsk = df["BER_bpsk"]

# =============================================================================
#  Create BER figure
# =============================================================================
plt.figure(figsize=(7.5, 6))

# --- LDPC BPSK ---
plt.semilogy(
    EbN0,
    ber_info,
    marker="o",
    markersize=8,
    markerfacecolor="none",
    markeredgewidth=1.8,
    linewidth=2.5,
    color="green",
    label="LDPC BPSK",
)

# --- Uncoded BPSK (theory) ---
plt.semilogy(
    EbN0,
    ber_bpsk,
    linewidth=2.5,
    color="red",
    label="Uncoded BPSK (theory)",
)


# Axis scale
plt.ylim(1e-5, 1)

# Labels
plt.xlabel("Eb/N0 [dB]", fontsize=18)
plt.ylabel("Bit Error Rate (BER)", fontsize=18)

# Grid
plt.grid(True, which="both", linestyle="--", linewidth=0.6, alpha=0.6)

# Legend
plt.legend(fontsize=14, loc="upper right", frameon=True, edgecolor="black")

plt.tight_layout()

# =============================================================================
#  Save PNG + SVG
# =============================================================================
plt.savefig("images/ldpc_ber_graph.png", dpi=300, bbox_inches="tight")
plt.savefig("images/ldpc_ber_graph.svg", bbox_inches="tight")  # 高品質ベクタ形式

plt.show()
