# Blind Demixing of Sparse Signals Diffused on Graphs

MATLAB solver for blind separation of sparse signals diffused on graphs. [1]

Diffusion of signals defined on the nodes of a graph is a generalization of the convolution from classical signal processing. [2] [3]

[1] [F. J. Iglesias](https://github.com/iglesias), S. Segarra, S. Rey-Escudero, A. G. Marques and D. Ramírez, *Demixing and Blind Deconvolution of Graph-Diffused Sparse Signals*. 2018 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), Calgary (2018), pp. 4189-4193 ([check in Google scholar](https://scholar.google.com/citations?view_op=view_citation&hl=en&user=H0okuHUAAAAJ&citation_for_view=H0okuHUAAAAJ:_xSYboBqXhAC)).

[2] D. I. Shuman, et al., *The Emerging Field of Signal Processing on Graphs: Extending High-Dimensional Data Analysis to Networks and Other Irregular Domains*. IEEE signal processing magazine 30.3 (2013): 83-98.

[3] A. Ortega, *[Introduction to Graph Signal Processing](https://www.graph-signal-processing-book.org/)*. Cambridge University Press (2022).


## ICASSP 2018

### Scenario 1. Single random graph (multiple diffusing filters on the same graph).

```Matlab
>> addpath('experiment_singlegraph_bss_logdet')
>> singlegraph_bss_logdet_N_S_numFilters
>> % Long for 1000 simulations!
>> Alternative: use pre-generated mat files in experiment_singlegraph_bss_logdet.

>> % Crunch data from mat files and produce mean and median RMSE plots.
>> cd experiment_singlegraph_bss_logdet
>> analysis_singlegraph_bss_logdet
```

### Scenario 2. Two coupled random graphs.
TODO

### Scenario 3. Multiple brain graphs.

```Matlab
>> addpath('experiment_brain_bss_logdet')
>> brain_bss_logdet_S_numGraphs
>> % Very long for 1000 simulations!
>> Alternative: use pre-generated mat files in experiment_brain_bss_logdet.

>> % Produce plot.
>> cd experiment_brain_bss_logdet
>> brain_bss_logdet_S_numGraphs_figure
```
