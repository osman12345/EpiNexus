"""
Motif Logo Visualization

Provides:
- Sequence logo generation from PWM/PFM
- Interactive motif visualization
- Motif comparison
- Information content calculation
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass


@dataclass
class MotifMatrix:
    """Position Weight Matrix for a TF motif."""

    name: str
    matrix: np.ndarray  # 4 x length (A, C, G, T)
    alphabet: str = "ACGT"
    source: str = ""

    @property
    def length(self) -> int:
        return self.matrix.shape[1]

    @property
    def consensus(self) -> str:
        """Get consensus sequence."""
        indices = np.argmax(self.matrix, axis=0)
        return "".join([self.alphabet[i] for i in indices])


class MotifLogoGenerator:
    """Generate sequence logos from motif matrices."""

    # Standard nucleotide colors (WebLogo style)
    COLORS = {
        "A": "#2ECC71",  # Green
        "C": "#3498DB",  # Blue
        "G": "#F39C12",  # Orange/Yellow
        "T": "#E74C3c",  # Red
    }

    # Alternative color schemes
    COLOR_SCHEMES = {
        "weblogo": {"A": "#2ECC71", "C": "#3498DB", "G": "#F39C12", "T": "#E74C3C"},
        "classic": {"A": "#00CC00", "C": "#0000CC", "G": "#FFB300", "T": "#CC0000"},
        "colorblind": {"A": "#009E73", "C": "#0072B2", "G": "#D55E00", "T": "#CC79A7"},
        "grayscale": {"A": "#404040", "C": "#808080", "G": "#606060", "T": "#A0A0A0"},
    }

    def __init__(self, color_scheme: str = "weblogo"):
        self.colors = self.COLOR_SCHEMES.get(color_scheme, self.COLOR_SCHEMES["weblogo"])

    def pfm_to_pwm(self, pfm: np.ndarray, pseudocount: float = 0.01, background: np.ndarray = None) -> np.ndarray:
        """Convert Position Frequency Matrix to Position Weight Matrix."""

        if background is None:
            background = np.array([0.25, 0.25, 0.25, 0.25])

        # Normalize PFM to frequencies
        pfm = pfm + pseudocount
        ppm = pfm / pfm.sum(axis=0, keepdims=True)

        # Calculate log-odds scores
        pwm = np.log2(ppm / background.reshape(-1, 1))

        return pwm

    def calculate_information_content(self, ppm: np.ndarray, background: np.ndarray = None) -> np.ndarray:
        """
        Calculate information content (bits) at each position.

        IC = 2 + sum(p * log2(p)) for DNA
        """

        if background is None:
            background = np.array([0.25, 0.25, 0.25, 0.25])

        # Avoid log(0)
        ppm = np.clip(ppm, 1e-10, 1.0)

        # Calculate entropy at each position
        entropy = -np.sum(ppm * np.log2(ppm), axis=0)

        # Maximum entropy for DNA = 2 bits
        max_entropy = 2.0

        # Information content = max_entropy - entropy
        ic = max_entropy - entropy

        return ic

    def calculate_height_matrix(self, ppm: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate letter heights for logo visualization.

        Returns:
            heights: Height of each letter at each position
            ic: Information content at each position
        """

        ic = self.calculate_information_content(ppm)

        # Scale probabilities by IC
        heights = ppm * ic

        return heights, ic

    def generate_logo_data(self, motif: MotifMatrix) -> Dict:
        """
        Generate data for logo visualization.

        Returns dict with position, letter, height, color for each element.
        """

        # Normalize matrix to PPM
        ppm = motif.matrix / motif.matrix.sum(axis=0, keepdims=True)
        heights, ic = self.calculate_height_matrix(ppm)

        logo_data = []

        for pos in range(motif.length):
            # Get heights and sort by height
            pos_heights = heights[:, pos]
            sorted_indices = np.argsort(pos_heights)

            # Stack letters from bottom to top
            y_offset = 0
            for idx in sorted_indices:
                letter = motif.alphabet[idx]
                height = pos_heights[idx]

                if height > 0.01:  # Skip very small letters
                    logo_data.append(
                        {
                            "position": pos + 1,
                            "letter": letter,
                            "height": height,
                            "y_start": y_offset,
                            "y_end": y_offset + height,
                            "color": self.colors[letter],
                            "ic": ic[pos],
                        }
                    )
                    y_offset += height

        return {
            "elements": logo_data,
            "length": motif.length,
            "max_ic": ic.max(),
            "total_ic": ic.sum(),
            "consensus": motif.consensus,
        }

    def generate_plotly_logo(self, motif: MotifMatrix, width: int = 600, height: int = 200):
        """Generate a Plotly figure for the motif logo."""
        import plotly.graph_objects as go

        logo_data = self.generate_logo_data(motif)

        fig = go.Figure()

        # Create letter shapes
        for elem in logo_data["elements"]:
            fig.add_trace(
                go.Scatter(
                    x=[elem["position"] - 0.4, elem["position"] - 0.4, elem["position"] + 0.4, elem["position"] + 0.4],
                    y=[elem["y_start"], elem["y_end"], elem["y_end"], elem["y_start"]],
                    fill="toself",
                    fillcolor=elem["color"],
                    line=dict(width=0),
                    mode="lines",
                    name=elem["letter"],
                    showlegend=False,
                    hoverinfo="text",
                    hovertext=f"Position: {elem['position']}<br>Letter: {elem['letter']}<br>Height: {elem['height']:.3f}",
                )
            )

            # Add letter text
            if elem["height"] > 0.1:
                fig.add_annotation(
                    x=elem["position"],
                    y=(elem["y_start"] + elem["y_end"]) / 2,
                    text=f"<b>{elem['letter']}</b>",
                    showarrow=False,
                    font=dict(
                        size=min(24, int(elem["height"] * 30)),
                        color="white" if elem["letter"] in ["G", "T"] else "white",
                        family="Arial Black",
                    ),
                )

        fig.update_layout(
            title=f"Motif: {motif.name}",
            xaxis_title="Position",
            yaxis_title="Information Content (bits)",
            width=width,
            height=height,
            xaxis=dict(tickmode="linear", tick0=1, dtick=1, range=[0.5, motif.length + 0.5]),
            yaxis=dict(range=[0, 2.2]),
            showlegend=False,
            plot_bgcolor="white",
        )

        return fig

    def compare_motifs(self, motif1: MotifMatrix, motif2: MotifMatrix) -> Dict[str, float]:
        """
        Compare two motifs using various metrics.

        Returns similarity scores.
        """

        # Normalize to PPM
        ppm1 = motif1.matrix / motif1.matrix.sum(axis=0, keepdims=True)
        ppm2 = motif2.matrix / motif2.matrix.sum(axis=0, keepdims=True)

        # Pad shorter motif
        len1, len2 = motif1.length, motif2.length
        max_len = max(len1, len2)

        if len1 < max_len:
            pad = np.full((4, max_len - len1), 0.25)
            ppm1 = np.hstack([ppm1, pad])
        if len2 < max_len:
            pad = np.full((4, max_len - len2), 0.25)
            ppm2 = np.hstack([ppm2, pad])

        # Pearson correlation
        corr = np.corrcoef(ppm1.flatten(), ppm2.flatten())[0, 1]

        # KL divergence (average)
        kl_div = np.mean(np.sum(ppm1 * np.log2(ppm1 / (ppm2 + 1e-10) + 1e-10), axis=0))

        # Euclidean distance
        euclidean = np.sqrt(np.sum((ppm1 - ppm2) ** 2))

        return {
            "pearson_correlation": corr,
            "kl_divergence": kl_div,
            "euclidean_distance": euclidean,
            "length_diff": abs(len1 - len2),
        }


# Pre-defined motifs for common TFs
KNOWN_MOTIFS = {
    "CTCF": MotifMatrix(
        name="CTCF",
        matrix=np.array(
            [
                [0.1, 0.1, 0.7, 0.1, 0.1, 0.1, 0.6, 0.1, 0.1, 0.7, 0.1, 0.1],  # A
                [0.2, 0.7, 0.1, 0.1, 0.7, 0.6, 0.1, 0.1, 0.7, 0.1, 0.6, 0.2],  # C
                [0.6, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.7, 0.1, 0.1, 0.2, 0.1],  # G
                [0.1, 0.1, 0.1, 0.7, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.6],  # T
            ]
        ),
        source="JASPAR",
    ),
    "SP1": MotifMatrix(
        name="SP1",
        matrix=np.array(
            [
                [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],  # A
                [0.1, 0.1, 0.8, 0.1, 0.8, 0.8, 0.1, 0.8, 0.8, 0.1],  # C
                [0.7, 0.7, 0.1, 0.8, 0.1, 0.1, 0.8, 0.1, 0.1, 0.8],  # G
                [0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1],  # T
            ]
        ),
        source="JASPAR",
    ),
    "MYC": MotifMatrix(
        name="MYC (E-box)",
        matrix=np.array(
            [
                [0.1, 0.7, 0.1, 0.1, 0.1, 0.7],  # A
                [0.7, 0.1, 0.8, 0.8, 0.1, 0.1],  # C
                [0.1, 0.1, 0.1, 0.1, 0.8, 0.1],  # G
                [0.1, 0.1, 0.0, 0.0, 0.0, 0.1],  # T
            ]
        ),
        source="JASPAR",
    ),
    "GATA1": MotifMatrix(
        name="GATA1",
        matrix=np.array(
            [
                [0.1, 0.1, 0.8, 0.1, 0.9, 0.1, 0.1],  # A
                [0.1, 0.1, 0.1, 0.1, 0.0, 0.1, 0.1],  # C
                [0.7, 0.8, 0.1, 0.8, 0.1, 0.1, 0.7],  # G
                [0.1, 0.0, 0.0, 0.0, 0.0, 0.7, 0.1],  # T
            ]
        ),
        source="JASPAR",
    ),
}


def get_motif(name: str) -> Optional[MotifMatrix]:
    """Get a known motif by name."""
    return KNOWN_MOTIFS.get(name.upper())


def list_available_motifs() -> List[str]:
    """List available pre-defined motifs."""
    return list(KNOWN_MOTIFS.keys())
