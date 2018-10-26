### icassp 2018

#### Scenario 1. Two coupled random graphs.
```Matlab
>> addpath('experiment_singlegraph_bss_logdet')
>> singlegraph_bss_logdet_N_S_numFilters
>> % Long for 1000 simulations!
>> Alternative: use pre-generated mat files in experiment_singlegraph_bss_logdet.

>> % Crunch data from mat files and produce mean and median RMSE plots.
>> cd experiment_singlegraph_bss_logdet
>> analysis_singlegraph_bss_logdet
```

#### Scenario 2. Two coupled random graphs.
TODO

#### Scenario 3. Multiple brain graphs.

```Matlab
>> addpath('experiment_brain_bss_logdet')
>> brain_bss_logdet_S_numGraphs
>> % Very long for 1000 simulations!
>> Alternative: use pre-generated mat files in experiment_brain_bss_logdet.

>> % Produce plot.
>> cd experiment_brain_bss_logdet
>> brain_bss_logdet_S_numGraphs_figure
```
