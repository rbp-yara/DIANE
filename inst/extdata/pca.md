# Principal Component Analysis

------------------------------------------------------------------------

Performing PCA on normalized RNA-Seq counts can be really informative about the conditions that impact gene expression the most. During PCA, new variables are computed, as linear combinations of your initial variables (e.g. experimental conditions). Those new variables, also called principal components, are designed to carry the maximum of the data variability.

We can plot the correlations of the initial variables to the new variables, the principal components, to see which ones contribute the most to the overall expression changes.

Each principal component explains a certain amount of the total variability, and those relative percentages are shown in what is called the screeplot.

In the first tab, you can have a quick view of the first 4 components. You can print specific components in the "Specific PCA plot". The "PCA correlation plot" shows correlations between groups defined in the experimental design and the different components.
