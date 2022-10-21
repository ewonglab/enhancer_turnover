# Nerual network model architecture and datasets for model

All models related code and datasets can be found here.

## Code requirement
Model requires Python 3.6+. Before you start, it might be easier to create your python evironment by using "requirement.txt" file.

```
python -m venv -r requirement.txt
```

## Modles & Datasets
The trained models we used for prediction can be found in the /models file folder, and datasets in /data folder. All files are organized by the below structure.

```
- /model

└── models
    ├── hg38
    │   ├── CEBPA
    │   └── Hnf4a
    └── mm10
        ├── CEBPA
        └── Hnf4a
- /data	

└── data
    ├── hg38
    │   ├── CEBPA
    │   └── Hnf4a
    └── mm10
        ├── CEBPA
        └── Hnf4a
```

For each TF and source speices:

- we trained 15 models and the best performance one is selected for further prediction.

- we seperate how genome window into three parts for training, held-out validation and held-out testing.
  - all data excluded chromosome 1 and 2 for model training(chr3toY_shuf.bed (chr3toY_pos_shuf.bed/chr3toY_neg_shuf.bed))
  - chromosome 1 for model validition (chr1_random_1m.bed)
  - chromosome 2 for model testing (chr2.bed)
  
 
If you want to train your own model, please refer to the oirginal domain adaptive model work, more detailed steps can be found in their [github](https://github.com/seqcode/cross-species-domain-adaptation).


