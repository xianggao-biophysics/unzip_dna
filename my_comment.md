1. While calculating LUT in CUDA, we can load DNA sequence to host memory
2. We can change the index from [j][ext] to [j][j-ext]. Before we do it, we need to investigate what region of the LUT is used, we can plot a **heatmap**.
3. We might need to use double buffer?
4. Think about the trace alignment using CCF
5. CCF lib to calculate CCF at different [stretch][shift]
