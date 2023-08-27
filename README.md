# DeepEpiTEX
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; DeepEpiTEX is a novel deep learning framework to characterize and infer the developmental hierarchy and functional states of exhausted T cells by integrating single- or multiple-layer epigenetic data modalities across 30 solid tumor types. This framework sheds light on potential epigenetic mechanisms governing differentiation and functional control, offering a reference atlas for the hierarchical categorization of TEX in various disease contexts. Furthermore, our model reveals potential associations between TEX subsets and immunotherapy responses, thus contributing to the advancement of personalized and precision immunotherapeutic strategies.

## Download
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 
An online platform, DeepEpiTEX (http://bio-data.cn/DeepEpiTEX/), offers a comprehensive resource for accessing models and exemplary datasets derived from different epigenetic data modalities. Due to the size limitations imposed by GitHub for file uploads, the model was unable to be uploaded here.
<img src="img/DeepEpiTEX_home.jpg" width="60%">


## Load model

``` python
## DeepEpiTEX constructed using the TensorFlow framework
from tensorflow.keras.models import load_model
classifier=load_model("LncRNA_MiRNA_Methylation")
```

## Example data

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;We incorporate distinct modality of features originating from diverse epigenetic data. Specifically, gene features prefixed with 'RNA_' are derived from RNA sequencing (RNA-seq) data, while those prefixed with 'Meth_' stem from DNA methylation data. Ensuring the accuracy and consistency of the model mandates a strong emphasis on maintaining feature order consistency. This practice is pivotal to correctly aligning features from disparate data sources with their respective data types, thereby ensuring precision in predictive outcomes and streamlining model maintenance and utilization procedures.

``` python
import pandas as pd
test_data=pd.read_csv("MiRNA_LncRNA_Meth_example.txt",sep=" ")
``` 
``` python
## Rows represent samples, and columns represent genes.
                Meth_ZYX	Meth_ZZEF1	Meth_ZZZ3  RNA_EHD4.AS1	RNA_RP11.166P13
TCGA.OR.A5JG.01	0.859171	0.895151	0.100498	10.584045	0.000000
TCGA.OR.A5L9.01	0.753927	0.804374	0.091260	0.000000	0.000000
```
## Predict TEX subsets

``` python
## Data Standardization
from sklearn.preprocessing import StandardScaler

sc = StandardScaler()
test_data = sc.fit_transform(test_data)

pre_label=pd.DataFrame(model.predict(test_data))
label=pre_label.idxmax(axis=1)

``` 
Within the labels, the values 0 through 4 correspond to TEX-S1 through TEX-S5.
