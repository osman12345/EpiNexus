"""
Interactive visualization components using Plotly.

Provides publication-quality plots for:
- Volcano plots
- MA plots
- PCA plots
- Heatmaps
- Chromatin state distributions
"""

import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from typing import Optional, List, Tuple


def create_volcano_plot(
    df: pd.DataFrame,
    fold_change_col: str = "Fold",
    fdr_col: str = "FDR",
    gene_col: str = "gene_name",
    fdr_threshold: float = 0.1,
    lfc_threshold: float = 0.5,
    title: str = "Volcano Plot",
    highlight_genes: Optional[List[str]] = None,
) -> go.Figure:
    """
    Create an interactive volcano plot.

    Args:
        df: DataFrame with differential peak results
        fold_change_col: Column name for log2 fold change
        fdr_col: Column name for FDR/adjusted p-value
        gene_col: Column name for gene names
        fdr_threshold: FDR significance threshold
        lfc_threshold: Log2 fold change threshold
        title: Plot title
        highlight_genes: List of gene names to highlight

    Returns:
        Plotly figure
    """
    df = df.copy()

    # Calculate -log10(FDR)
    df["neg_log10_fdr"] = -np.log10(df[fdr_col].clip(lower=1e-300))

    # Classify points
    def classify_point(row):
        if row[fdr_col] >= fdr_threshold:
            return "Not Significant"
        elif row[fold_change_col] > lfc_threshold:
            return "Gained"
        elif row[fold_change_col] < -lfc_threshold:
            return "Lost"
        else:
            return "Not Significant"

    df["category"] = df.apply(classify_point, axis=1)

    # Color mapping
    color_map = {
        "Gained": "#E41A1C",  # Red
        "Lost": "#377EB8",  # Blue
        "Not Significant": "#999999",
    }

    # Create figure
    fig = px.scatter(
        df,
        x=fold_change_col,
        y="neg_log10_fdr",
        color="category",
        color_discrete_map=color_map,
        hover_data=[gene_col, fdr_col, fold_change_col],
        title=title,
        labels={fold_change_col: "Log2 Fold Change", "neg_log10_fdr": "-Log10(FDR)"},
    )

    # Add threshold lines
    fig.add_hline(
        y=-np.log10(fdr_threshold), line_dash="dash", line_color="gray", annotation_text=f"FDR = {fdr_threshold}"
    )
    fig.add_vline(x=lfc_threshold, line_dash="dash", line_color="gray")
    fig.add_vline(x=-lfc_threshold, line_dash="dash", line_color="gray")

    # Highlight specific genes
    if highlight_genes:
        highlight_df = df[df[gene_col].isin(highlight_genes)]
        fig.add_trace(
            go.Scatter(
                x=highlight_df[fold_change_col],
                y=highlight_df["neg_log10_fdr"],
                mode="markers+text",
                marker=dict(size=12, color="gold", symbol="star"),
                text=highlight_df[gene_col],
                textposition="top center",
                name="Highlighted",
                hoverinfo="text",
            )
        )

    # Layout
    fig.update_layout(template="plotly_white", legend_title="Category", hovermode="closest", width=800, height=600)

    return fig


def create_ma_plot(
    df: pd.DataFrame,
    fold_change_col: str = "Fold",
    mean_col: str = "Conc",
    fdr_col: str = "FDR",
    gene_col: str = "gene_name",
    fdr_threshold: float = 0.1,
    title: str = "MA Plot",
) -> go.Figure:
    """
    Create an interactive MA plot (log2FC vs mean expression).

    Args:
        df: DataFrame with differential results
        fold_change_col: Column for log2 fold change
        mean_col: Column for mean concentration/expression
        fdr_col: Column for FDR
        gene_col: Column for gene names
        fdr_threshold: Significance threshold
        title: Plot title

    Returns:
        Plotly figure
    """
    df = df.copy()

    # Classify significance
    df["significant"] = df[fdr_col] < fdr_threshold
    df["category"] = np.where(
        ~df["significant"], "Not Significant", np.where(df[fold_change_col] > 0, "Gained", "Lost")
    )

    color_map = {"Gained": "#E41A1C", "Lost": "#377EB8", "Not Significant": "#999999"}

    fig = px.scatter(
        df,
        x=mean_col,
        y=fold_change_col,
        color="category",
        color_discrete_map=color_map,
        hover_data=[gene_col, fdr_col],
        title=title,
        labels={mean_col: "Mean Concentration (log2)", fold_change_col: "Log2 Fold Change"},
    )

    # Add zero line
    fig.add_hline(y=0, line_dash="solid", line_color="black", line_width=1)

    fig.update_layout(template="plotly_white", width=800, height=600)

    return fig


def create_pca_plot(
    df: pd.DataFrame,
    pc1_col: str = "PC1",
    pc2_col: str = "PC2",
    color_col: str = "condition",
    label_col: str = "sample",
    variance_explained: Optional[Tuple[float, float]] = None,
    title: str = "PCA Plot",
) -> go.Figure:
    """
    Create an interactive PCA plot.

    Args:
        df: DataFrame with PCA coordinates
        pc1_col: Column for PC1
        pc2_col: Column for PC2
        color_col: Column for coloring points
        label_col: Column for point labels
        variance_explained: Tuple of (PC1_var, PC2_var) percentages
        title: Plot title

    Returns:
        Plotly figure
    """
    # Labels
    xlabel = f"PC1 ({variance_explained[0]:.1f}%)" if variance_explained else "PC1"
    ylabel = f"PC2 ({variance_explained[1]:.1f}%)" if variance_explained else "PC2"

    fig = px.scatter(
        df,
        x=pc1_col,
        y=pc2_col,
        color=color_col,
        text=label_col,
        title=title,
        labels={pc1_col: xlabel, pc2_col: ylabel},
    )

    fig.update_traces(textposition="top center", marker=dict(size=12))

    fig.update_layout(template="plotly_white", width=700, height=600)

    return fig


def create_heatmap(
    data: pd.DataFrame, title: str = "Correlation Heatmap", color_scale: str = "RdBu_r", show_values: bool = True
) -> go.Figure:
    """
    Create a correlation/expression heatmap.

    Args:
        data: DataFrame (samples x features) or correlation matrix
        title: Plot title
        color_scale: Plotly color scale
        show_values: Show values on heatmap

    Returns:
        Plotly figure
    """
    # If not a correlation matrix, compute it
    if data.shape[0] != data.shape[1]:
        data = data.corr()

    fig = go.Figure(
        data=go.Heatmap(
            z=data.values,
            x=data.columns,
            y=data.index,
            colorscale=color_scale,
            text=data.values.round(2) if show_values else None,
            texttemplate="%{text}" if show_values else None,
            textfont={"size": 10},
            hoverongaps=False,
        )
    )

    fig.update_layout(title=title, template="plotly_white", width=700, height=600)

    return fig


def create_chromatin_state_plot(
    df: pd.DataFrame, state_col: str = "chromatin_state", title: str = "Chromatin State Distribution"
) -> go.Figure:
    """
    Create a bar plot of chromatin state distribution.

    Args:
        df: DataFrame with chromatin state assignments
        state_col: Column name for chromatin states
        title: Plot title

    Returns:
        Plotly figure
    """
    # Count states
    state_counts = df[state_col].value_counts().reset_index()
    state_counts.columns = ["State", "Count"]

    # Define colors for each state
    state_colors = {
        "Active_enhancer": "#4DAF4A",
        "Active_promoter": "#377EB8",
        "Poised_enhancer": "#FF7F00",
        "Bivalent": "#984EA3",
        "Repressed": "#E41A1C",
        "Complex": "#A65628",
        "Other": "#999999",
    }

    # Map colors
    state_counts["Color"] = state_counts["State"].map(lambda x: state_colors.get(x, "#999999"))

    fig = go.Figure(
        data=[
            go.Bar(
                x=state_counts["State"],
                y=state_counts["Count"],
                marker_color=state_counts["Color"],
                text=state_counts["Count"],
                textposition="outside",
            )
        ]
    )

    fig.update_layout(
        title=title,
        xaxis_title="Chromatin State",
        yaxis_title="Number of Genes",
        template="plotly_white",
        showlegend=False,
        width=800,
        height=500,
    )

    return fig


def create_integration_scatter(
    df: pd.DataFrame,
    x_col: str = "ac_fold_change",
    y_col: str = "me3_fold_change",
    color_col: str = "integration_category",
    gene_col: str = "gene_symbol",
    title: str = "H3K27ac vs H3K27me3 Changes",
) -> go.Figure:
    """
    Create scatter plot for two-mark integration.

    Args:
        df: Integrated results DataFrame
        x_col: Column for x-axis (e.g., H3K27ac fold change)
        y_col: Column for y-axis (e.g., H3K27me3 fold change)
        color_col: Column for coloring points
        gene_col: Column for gene names
        title: Plot title

    Returns:
        Plotly figure
    """
    # Color mapping for categories
    color_map = {
        "Strong_activation": "#1b9e77",
        "Strong_repression": "#d95f02",
        "Activation_ac_only": "#7fc97f",
        "Activation_me3_only": "#beaed4",
        "Repression_ac_only": "#fdc086",
        "Repression_me3_only": "#ffff99",
        "Discordant_both_gained": "#984ea3",
        "Discordant_both_lost": "#e7298a",
        "No_significant_change": "#999999",
    }

    fig = px.scatter(
        df,
        x=x_col,
        y=y_col,
        color=color_col,
        color_discrete_map=color_map,
        hover_data=[gene_col],
        title=title,
        labels={x_col: "H3K27ac Log2 Fold Change", y_col: "H3K27me3 Log2 Fold Change"},
    )

    # Add quadrant lines
    fig.add_hline(y=0, line_dash="dash", line_color="gray")
    fig.add_vline(x=0, line_dash="dash", line_color="gray")

    # Add quadrant labels
    fig.add_annotation(x=2, y=2, text="Both Gained", showarrow=False, font=dict(color="gray"))
    fig.add_annotation(x=-2, y=-2, text="Both Lost", showarrow=False, font=dict(color="gray"))
    fig.add_annotation(x=2, y=-2, text="Activation", showarrow=False, font=dict(color="green"))
    fig.add_annotation(x=-2, y=2, text="Repression", showarrow=False, font=dict(color="red"))

    fig.update_layout(template="plotly_white", width=800, height=700)

    return fig


def create_upset_plot(df: pd.DataFrame, set_columns: List[str], title: str = "Mark Overlap") -> go.Figure:
    """
    Create an UpSet-style plot for set intersections.

    Args:
        df: DataFrame with boolean columns for set membership
        set_columns: List of column names representing sets
        title: Plot title

    Returns:
        Plotly figure (simplified bar chart version)
    """
    # Calculate intersections
    intersections = []
    for i in range(1, 2 ** len(set_columns)):
        # Convert to binary to get combination
        combo = []
        for j, col in enumerate(set_columns):
            if i & (1 << j):
                combo.append(col)

        # Filter rows matching this exact combination
        mask = pd.Series([True] * len(df))
        for col in set_columns:
            if col in combo:
                mask &= df[col]
            else:
                mask &= ~df[col]

        count = mask.sum()
        if count > 0:
            intersections.append({"combination": " + ".join(combo), "count": count, "n_sets": len(combo)})

    int_df = pd.DataFrame(intersections)
    int_df = int_df.sort_values("count", ascending=False)

    fig = px.bar(
        int_df,
        x="combination",
        y="count",
        color="n_sets",
        color_continuous_scale="Viridis",
        title=title,
        labels={"combination": "Mark Combination", "count": "Gene Count", "n_sets": "# Marks"},
    )

    fig.update_layout(template="plotly_white", xaxis_tickangle=45, width=900, height=500)

    return fig
