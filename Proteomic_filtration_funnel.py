import matplotlib.pyplot as plt
import matplotlib.patches as patches

def generate_refined_funnel():
    # 1. Define Updated Data
    stages = [
        {"label": "Brain Lysosome Atlas", "n": 5353, "color": "#90afc5"},
        {"label": "LSD-Associated", "n": 647, "color": "#afbc9f"},
        {"label": "Palmitoylated Candidates\n(Verified + Predicted)", "n": 626, "color": "#d3abb1", 
         "sub": "Verified: 432 | Predicted: 194"}, # Combined to one line to prevent escaping
        {"label": "Top 10 PPI Clusters", "n": 151, "color": "#a84646"},
        {"label": "Allen Brain Atlas\nSpatiotemporal mapping", "n": 109, "color": "#7b6094"}
    ]
    
    excluded_counts = [4706, 21, 475, 42]
    excluded_reasons = [
        "Non-LSD Associated",
        "Non-Palmitoylated",
        "Other PPI Clusters",
        "Unmapped Proteins"
    ]

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Dimensions for a "proper" funnel shape
    top_width = 100
    bottom_width = 35
    height_per_tier = 20
    gap = 1.5
    current_y = 100
    
    # Calculate width step for a linear taper
    width_step = (top_width - bottom_width) / len(stages)

    for i, stage in enumerate(stages):
        # Calculate trapezoid dimensions
        w_top = top_width - (i * width_step)
        w_bot = top_width - ((i + 1) * width_step)
        
        # Coordinates (centered at 0)
        coords = [
            [-w_top/2, current_y], 
            [w_top/2, current_y], 
            [w_bot/2, current_y - height_per_tier], 
            [-w_bot/2, current_y - height_per_tier]
        ]
        
        # Draw the Tier
        poly = patches.Polygon(coords, facecolor=stage['color'], edgecolor='white', linewidth=1.5)
        ax.add_patch(poly)
        
        # Text Positioning
        mid_y = current_y - (height_per_tier / 2)
        text_color = 'black' if i < 3 else 'white'
        
        # Main Label (slightly smaller font to ensure fit)
        ax.text(0, mid_y + 2.5, stage['label'], ha='center', va='center', 
                fontsize=11, fontweight='bold', color=text_color)
        
        # N value
        ax.text(0, mid_y - 2.5, f"n = {stage['n']:,}", ha='center', va='center', 
                fontsize=14, fontweight='bold', color=text_color)
        
        # Sub-text for Palmitoylation (kept inside the box)
        if 'sub' in stage:
            ax.text(0, mid_y - 7, stage['sub'], ha='center', va='center', 
                    fontsize=9, color='black', style='italic', fontweight='medium')

        # Exclusion Arrows and Labels
        if i < len(excluded_counts):
            arrow_y = current_y - height_per_tier
            # Dynamically anchor arrow to the edge of the funnel
            ax.annotate('', xy=(w_top/2 + 15, arrow_y), xytext=(w_bot/2, arrow_y),
                        arrowprops=dict(arrowstyle="->", color="#555555", lw=1))
            
            ex_text = f"Excluded:\n{excluded_reasons[i]}\n(n = {excluded_counts[i]:,})"
            ax.text(w_top/2 + 17, arrow_y, ex_text, va='center', fontsize=10, color='#444444', linespacing=1.2)

        current_y -= (height_per_tier + gap)

    # Styling
    ax.set_xlim(-75, 140)
    ax.set_ylim(0, 110)
    ax.axis('off')
    plt.title("Proteomic Refinement Pipeline", fontsize=22, fontweight='bold', pad=30)
    
    # Save as high-resolution PNG
    plt.savefig('Proteomic_Refinement_Funnel_HighRes.png', dpi=300, bbox_inches='tight', transparent=False, facecolor='white')
    
    print("Funnel generated and saved as 'Proteomic_Refinement_Funnel_HighRes.png'")
    plt.show()

if __name__ == "__main__":
    generate_refined_funnel()