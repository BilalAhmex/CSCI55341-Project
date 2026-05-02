# SpatialDE Analysis: Our Results vs. Github Results

This document explains the expected numerical differences between the results generated in this pipeline and the examples provided in the SpatialDE authors' github.


---

## 1. The Length Scale ($l \approx 1.41$ vs $l \approx 155.73$)

The length scale (`l`) metric in SpatialDE simply reflects the **units of measurement** of the X and Y coordinates. SpatialDE does not inherently know what a millimeter or a pixel is; it only evaluates the numerical range it is fed.

*   **The Authors' Data ($l = 1.41$):** The sample data provided for the SpatialDE tutorials features cells arranged on a perfectly uniform, integer grid. The spots are located at exact integer coordinates such as `(1, 1)`, `(1, 2)`, `(2, 2)`. In a grid system, the diagonal distance between two adjacent spots is exactly $\sqrt{2}$, which equals **1.41**. Their model is finding patterns that span roughly 1.41 grid spots.
*   **Our seqFISH Data ($l \approx 155.73$):** Our data is not mapped to an integer grid. It utilizes raw, continuous physical imaging coordinates that scale up to ~1,000. Therefore, our model identifies spatial patterns that span roughly 156 pixels/micrometers across the actual tissue.

**Takeaway:** Neither scale is "wrong." They merely reflect the different underlying measuring sticks used in the respective coordinate data.

---

## 2. Statistical Strength

Our code produced highly definitive statistical results (e.g., `qval = 0.0`), compared to the authors' sample snippet (e.g., `qval = 0.43`). For Seqfish:

*   **Sample Size ($N$):** The tutorial sample processes an $N$ of **257** cells/spots. In contrast, our dataset processes an $N$ of **3,585** cells/spots.
*   **Statistical Power:** The more data points provided to the model, the higher the mathematical confidence in the output. Because our dataset contains roughly 14 times more data points than the tutorial, the resulting p-values shrink from "marginally significant" to "absolutely, undeniably significant".
*   **Biological Resolution:** seqFISH is a modern, single-cell resolution technology. I am assuming the dataset used in the authors' quick tutorial is an older, lower-resolution array that was heavily downsampled so users could successfully run the tutorial on a standard laptop in a few seconds.

---
