# vision_pipeline_hg19_0526

## Dependence:
#### Python/2.7
###### numpy, matplotlib, scipy
#### R



## Input data
###### cell_list.txt in input folder
```
cell_list.txt
>>> head cell_list.txt 
H1
H2
CD34

```

###### mark_list.txt in input folder
```
mark_list.txt
>>> head mark_list.txt 
H3K4me1
H3K4me3
H3K27ac

```

###### replicate_list.txt in input folder
```
replicate_list.txt
>>> head replicate_list.txt 
R1
R2

```

## run pipeline
###### The 1st parameter is the script folder of the pipeline; 2nd parameter is the input folder of the data
###### The 3rd parameter is the reference sample for across mark normalization (e.g. H3K4me3_H1_R1)
###### The 4th parameter is the reference cell type within mark normalization (e.g. H1)
###### The 5th parameter is the reference replicate id within mark normalization (e.g. R1)
```
time sh /storage/home/software/shared_pipeline/get_vision_human_rep.sh /storage/home/software/shared_pipeline/ /storage/home/data2run/test_data/ H3K4me3_H1_R1 H1 R1
```

