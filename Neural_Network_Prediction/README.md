# Nerual network model architecture and datasets for model

All models related code and datasets can be found here.

## Code requirement
Script requires Python 3.6+. Before you start, it might be easier to create your python evironment by using "requirement.txt" file. ```python -m venv -r requirement.txt```

Before you repeat all analysis, please run ```git clone https://github.com/seqcode/cross-species-domain-adaptation.git``` first to make sure you have all needed files under your working directory.

## Models & Datasets

For each TF and source speices:

- we trained 15 models and the best performance one is selected for further prediction.

- we seperate how genome window into three parts for training, held-out validation and held-out testing.
  - all data excluded chromosome 1 and 2 for model training(chr3toY_shuf.bed (including chr3toY_pos_shuf.bed and chr3toY_neg_shuf.bed))
  - chromosome 1 for model validition (chr1_random_1m.bed)
  - chromosome 2 for model testing (chr2.bed)
  

The **best trained models** we used for prediction can be found in the models.tar.gz file. Due to file size limitation, you can access it from the [link](https://drive.google.com/file/d/1h3egck0zs-d7TsbJpkNQUrtMWGiI33HO/view?usp=sharing). The model can directely use for TF binding prediction.

- CEBPA_hg38trained.model
- CEBPA_mm10trained.model
- Hnf4a_hg38_trained.model
- Hnf4a_mm10trained.model

**All datasets files** are organized by the below structure, you can easily re-create from the github script according to our paper. 

```
└── data
    ├── hg38
    │   ├── CEBPA
    │   └── Hnf4a
    └── mm10
        ├── CEBPA
        └── Hnf4a
```
 
If you want to train your own model, please refer to the oirginal domain adaptive model work, more detailed steps can be found in their [github](https://github.com/seqcode/cross-species-domain-adaptation).

## Datasets for prediction

Datasets are used in "Evolutionary sequence divergence predicts functional changes to TF binding in enhancer turnover" section.

There are two sets of data.
- human-specific enhancers and their matched non-enhancer counterpart in mouse
- mouse-specific enhancers and their matched non-enhancer counterpart in human

You can find the files in the data folder.

## Analysis Scripts

A few scripts we created to do model prediction and further analysis can be found in the analysis folder.

- model_predict: do model prediction and save predicted values in *.npy.
- piechart.ipynb: draw piechart for bind prediction for human/mouse datasets respectively.
- roc_plot.ipynb: ROC curve of our trained models.

Please refer to our paper "DNA replication timing predicts transcription factor binding and the rate of enhancer turnover" for more details.

