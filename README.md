# Using 
For the Mouse Olfactory Bulb: Use MOB Rep 11  (or rep 12)
For the Breast Cancer: Use BC Layer 1 (or Layer 2). T


# Results for MOB
It is completely normal to be suspicious of a flat 0.0 in statistics, but your results are exactly what they should be. Here is why you are seeing those zeros, and why that weird 1.11e-16 number you spotted is the ultimate proof that your run was successful.

1. The Computer Science Reason: "Machine Epsilon"

The math behind SpatialDE calculates how likely it is that a gene's expression pattern is just random noise. For genes like Apoe or Trak2, the biological pattern is so incredibly strong that the probability of it being random is something like 0.000000000000000000000000000001.

Standard computers use 64-bit floating-point memory to store decimals. There is a hard limit to how small a decimal a computer can track before it just runs out of memory space.

That list of numbers you pasted starting with 1.110223e-16? That is a very famous number in computer science called Machine Epsilon. It is the absolute smallest difference between numbers that Python can calculate.

When SpatialDE calculates a p-value that is mathematically smaller than that 1.11e-16 threshold, Python literally runs out of decimal places, throws its hands up, and rounds it to an absolute 0.0.

So, 0.0 doesn't mean "error." It means "a signal so strong it broke the computer's decimal limit."

2. The Biological Reason: Famous Marker Genes

The genes that hit that 0.0 limit in your results are actually the ultimate biological proof that your algorithm worked.

If you look up the top genes your model spit out (Apoe, Apc, Fabp7), they are not random. Apoe and Fabp7 are famous, heavily researched genes in the mouse brain. They are heavily restricted to specific layers of glial cells and astrocytes. If SpatialDE didn't flag them with the highest possible significance, we would know the algorithm was broken!

3. Q-Values Follow P-Values

The qval is just the False Discovery Rate (FDR) adjusted version of the pval. Because the initial pval was rounded to 0.0, the math to adjust it also resulted in 0.0.

What you should do next:
You have a perfect output file for SpatialDE. You can confidently take this spatialDE_Rep1_MOB_results.csv, look at all the genes with a qval < 0.05, and consider them your "Found Spatially Variable Genes" to compare against SPARK and your SVG ground truth file!