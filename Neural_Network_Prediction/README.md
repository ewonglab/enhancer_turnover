You can find all details about Nerual network model architecture and datasets for model training, validation and testing.

the trained models we used for prediction can be found in the /models file folder, and datasets in /data folder.

For each TF and source speices:

- we trained 15 models and the best performance one is selected for further prediction.

- we seperate how genome window into three parts for training, held-out validation and held-out testing.
  - all data excluded chromosome 1 and 2 for model training()
  - chromosome 1 for model validition ()
  - chromosome 2 for model testing ()

Detailed datasets info are listed in the below table.
Table I. Training data
Training data
Cebpa
Hnf4a
	Total windows	Bind windows 	Total windows	Bind windows 
Human windows	159,304
79,652
260,372
130,186

Mouse windows	332,866
166,433
199,792
99,896

![Uploading image.png…]()

 
If you want to train your own model, please refer to the oirginal domain adaptive model work, more detailed steps can be found in their [github].(https://github.com/seqcode/cross-species-domain-adaptation)
