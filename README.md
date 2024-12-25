# PDATC-NCPMKL-updated
This repository provides the data and codes for the prediction of drug-ATC code associations. An inference method is first designed to identify highly related target proteins, fingerprints, and side effects of ATC codes, which adopts random walk with restart algorithm and permutation test. The procedures are illustrated in the following figure.
![Figure1](https://github.com/Lywhere/PDATC-NCPMKL-updated/blob/main/model/Figure1.png)

Based on the highly related target proteins, fingerprints, and side effects of ATC codes, three new ATC code kernels are constructed, which are fused with the ATC code kernels in one previous model. The updated model provides high performance with AUROC and AUPR exceeding 0.96. 
![Figure3](https://github.com/Lywhere/PDATC-NCPMKL-updated/blob/main/model/Figure3.png)
# Requirements
+ python = 3.10
+ cvxopt = 1.3.2
+ ecos = 2.0.12
+ numpy+mkl = 1.26.2
+ osqp = 0.6.3
+ pandas = 1.3.4
+ scikit-learn = 1.3.2
+ scipy = 1.11.3
+ scs = 3.2.3
# Files:
 filename | explain
 ------------------------- | -------------------------
 drug_ATC | The adjacency matrix of drugs and ATC codes at the second,third,fourth level.
 drug_fingerprint | The drug representation based on their fingerprints.
 drug_side_effects | The drug representation based on their side effects.
 drug_target_protein | The drug representation based on their target proteins.
 drug_interaction | The drug kernel using the interaction information collected in STITCH.
 ATC_fingerprint | The ATC code representation based on their fingerprints at the second, third, fourth level.
 ATC_side_effects | The ATC code representation based on their side effects at the second, third, fourth level.
 ATC_target_protein | The ATC code representation based on their target proteins at the second, third, fourth level.

Drug SMILES information, drug ATC code information, and drug target protein information were obtained from DrugBank database (https://go.drugbank.com/), drug side effect information were obtained from SIDER database (http://sideeffects.emblde/), and drug interaction information was obtained from STITCH website (http://stitch4.embl.de/).
# Usage
## How to use it?
## 1. Use the data set we provide
### 1.1 Cross verification
If you use our dataset for cross-validation, all you need to do is enter the following command in the termi

    python main.py 

### 1.2 Modify model parameters 
You just need to adjust the following code in the main.py file.

```python
if __name__ == "__main__":
    drug_atc_path = 'data/drug_ATC/second_ATC.csv'
    atc_target_protein_path = 'data/ATC_target_protein/second_uniprot.csv'
    atc_side_effects_path = 'data/ATC_side_effects/second_side_effects.csv'
    atc_fingerprint_path = 'data/ATC_fingerprint/second_fingerprint.csv'
    op = Options(drug_atc_path=drug_atc_path, atc_target_protein_path=atc_target_protein_path, atc_side_effects_path=atc_side_effects_path, atc_fingerprint_path=atc_fingerprint_path, level=4, omega=0.9)
    op.train(k=10)
```

+ **drug_atc_path** is the file path storing the adjacency matrix of drug and ATC codes.
+ **atc_target_protein_path** is the file path storing the adjacency matrix of ATC codes and target protein.
+ **atc_side_effects_path** is the file path storing the adjacency matrix of ATC codes and side effects.
+ **atc_fingerprint_path** is the file path storing the adjacency matrix of ATC codes and fingerprint.
+ The parameter level represents the level of ATC codes. **It can be 2, 3 and 4**.
+ The parameter omega represents parameter in WKNKN when reformulating the adjacency matrix. It can be any numbers **between 0.0 and 1.0**.
+ The parameter k represents the number of folds in cross-validation. **k was set to 10 in our study**.
## 2. Use your own data set
### 2.1 Preprocessed data set
You need to prepare some files, which are all in CSV format. The detailed format is displayed as below:
#### 1. The adjacency matrix of drug-ATC code associations
DrugBankID | ATCcode*1* | ATCcode*2* | ATCcode*3* | ATCcode*4* | ... | ATCcode*m*
------ | ------ | ------ | ------ | ------ | ------ | -------
drugID*1* | 0 | 1 | 1 | 0 | ... | 0
drugID*2* | 1 | 0 | 0 | 1 | ... | 0
drugID*3* | 1 | 1 | 0 | 0 | ... | 1
drugID*4* | 1 | 1 | 0 | 0 | ... | 1
... | ... | ... | ... | ... | ... | ...
drugID*n* | 0 | 0 | 1 | 1 | ... | 0

#### 2. Drug fingerprints matrix
DrugBankID | F*1* | F*2* | F*3* | F*4* | ... 
------ | ------ | ------ | ------ | ------ | ------ 
drugID*1* | 0 | 1 | 1 | 0 | ... 
drugID*2* | 1 | 0 | 0 | 1 | ... 
drugID*3* | 1 | 1 | 0 | 0 | ... 
drugID*4* | 1 | 1 | 0 | 0 | ... 
... | ... | ... | ... | ... | ... 
drugID*n* | 0 | 0 | 1 | 1 | ..

#### 3. Drug side effects matrix
DrugBankID | side*1* | side*2* | side*3* | side*4* | ... 
------ | ------ | ------ | ------ | ------ | ------ 
drugID*1* | 1 | 0 | 0 | 1 | ... 
drugID*2* | 1 | 1 | 0 | 0 | ... 
drugID*3* | 0 | 0 | 0 | 1 | ... 
drugID*4* | 1 | 1 | 0 | 0 | ... 
... | ... | ... | ... | ... | ... 
drugID*n* | 0 | 0 | 1 | 1 | ... 

#### 4. Drug target proteins matrix
DrugBankID | side*1* | side*2* | side*3* | side*4* | ... 
------ | ------ | ------ | ------ | ------ | ------ 
drugID*1* | 0 | 1 | 0 | 1 | ... 
drugID*2* | 1 | 0 | 0 | 0 | ... 
drugID*3* | 0 | 0 | 0 | 1 | ... 
drugID*4* | 0 | 0 | 1 | 1 | ... 
... | ... | ... | ... | ... | ... 
drugID*n* | 0 | 0 | 1 | 0 | ... 

#### 5. Drug interaction kernel
DrugBankID | drugID*1* | drugID*2* | drugID*3* | drugID*4* | ... | drugID*n*
------ | ------ | ------ | ------ | ------ | ------ | -------
drugID*1* | 1 | 0.3 | 0 | 0.75 | ... | 0.33
drugID*2* | 0.3 | 1 | 0.9 | 0.22 | ... | 0.68
drugID*3* | 0 | 0.9 | 1 | 0 | ... | 0.47
drugID*4* | 0.75 | 0.22 | 0 | 1 | ... | 0.92
... | ... | ... | ... | ... | ... | ...
drugID*n* | 0.33 | 0.68 | 0.47 | 0.92 | ... | 1

#### 6. ATC code fingerprints matrix
ATCcodeID | F*1* | F*2* | F*3* | F*4* | ...
------ | ------ | ------ | ------ | ------ | ------
ATCcode*1* | 0 | 1 | 0 | 1 | ... 
ATCcode*2* | 1 | 0 | 0 | 0 | ... 
ATCcode*3* | 0 | 0 | 0 | 1 | ... 
ATCcode*4* | 0 | 0 | 1 | 1 | ... 
... | ... | ... | ... | ... | ... 
ATCcode*m* | 0 | 0 | 1 | 0 | ... 

#### 7. ATC code side effects matrix
ATCcodeID | side*1* | side*2* | side*3* | side*4* | ... 
------ | ------ | ------ | ------ | ------ | ------ 
ATCcode*1* | 0 | 1 | 1 | 0 | ... 
ATCcode*2* | 1 | 0 | 0 | 1 | ... 
ATCcode*3* | 1 | 1 | 0 | 0 | ... 
ATCcode*4* | 1 | 1 | 0 | 0 | ... 
... | ... | ... | ... | ... | ... 
ATCcode*m* | 0 | 0 | 1 | 1 | ... 

#### 8. ATC code target proteins matrix
ATCcodeID | F*1* | F*2* | F*3* | F*4* | ... 
------ | ------ | ------ | ------ | ------ | ------ 
ATCcode*1* | 1 | 0 | 0 | 1 | ... 
ATCcode*2* | 1 | 1 | 0 | 0 | ... 
ATCcode*3* | 0 | 0 | 0 | 1 | ... 
ATCcode*4* | 1 | 1 | 0 | 0 | ... 
... | ... | ... | ... | ... | ... 
ATCcode*m* | 0 | 0 | 1 | 1 | ... 

#### 9. Because it involves ATC code tree structure to find the shortest path, different data sets involve different ATC, so you should prepare the shortest path file for ATC codes at different levels. It is also in CSV format, as shown below
ATCcodeID | ATCcode*1* | ATCcode*2* | ATCcode*3* | ATCcode*4* | ... | ATCcode*m*
------ | ------ | ------ | ------ | ------ | ------ | -------
ATCcode*1* | 0 | 2 | 2 | 4 | ... | 4
ATCcode*2* | 2 | 0 | 2 | 8 | ... | 6
ATCcode*3* | 2 | 2 | 0 | 4 | ... | 8
ATCcode*4* | 4 | 8 | 4 | 0 | ... | 2
... | ... | ... | ... | ... | ... | ...
ATCcode*m* | 4 | 6 | 8 | 2 | ... | 0

+ You should put this file in the **NAME/shortest_path/** folder, and it should have **the same file name as mine. (For example, the second level ATC code file is named new_2ATC_shortest_path_length_matrix.csv)**
+ In addition, in order to prevent the accuracy of SPro kernel matrix calculation, ensure that the order of ATCcode here is consistent with that in **the adjacency matrix of drug-ATC code**.
### 2.2 Cross verification
You just need to modify the following code in the main.py file to run it:

```python
def file_path(self):
    drug_fingerprint_path = 'your own drug fingerprint file path'
    drug_side_effects_path = 'your own drug side effect file path'
    drug_target_protein_path = 'your own drug target protein file path'
    drug_interaction_path = 'your own drug interaction file path'
    return drug_fingerprint_path, drug_side_effects_path, drug_target_protein_path, drug_interaction_path
```

```python
if __name__ == "__main__":
    drug_atc_path = 'your drug-ATC code adjacency matrix file path'
    atc_target_protein_path = 'your ATC-target protein adjacency matrix file path'
    atc_side_effects_path = 'your ATC-side effect adjacency matrix file path'
    atc_fingerprint_path = 'your ATC-fingerprint adjacency matrix file path'
    op = Options(drug_atc_path=drug_atc_path, atc_target_protein_path=atc_target_protein_path, atc_side_effects_path=atc_side_effects_path, atc_fingerprint_path=atc_fingerprint_path, level=4, omega=0.9)
    op.train(k=10)
```
    
### The results predicted by the model
After running our model, the **NAME_predict.csv** file and **NAME_actual.csv** file will be generated, where the **NAME_predict.csv** file will store the **predicted score**, the **actual value** is saved in the **NAME_actual.csv** file.

# Additional experiment
There were six ATC code kernels in this model. It was necessary to analyze their importance for the model. In view of this, each ATC code kernel was removed one by one, producing six models. These models were still evaluated by ten-fold cross-validation. 
```python
atc_global_kernel = atc_matrix_combination([ATC_atc_columns_kernel, ATC_atc_probabilistic_kernel, ATC_atc_SM_kernel, ATC_target_protein_kernel, ATC_side_effects_kernel, ATC_fingerprint_kernel])
```
## The models by removing one ATC code kernel
```python
atc_global_kernel = atc_matrix_combination([ATC_atc_probabilistic_kernel, ATC_atc_SM_kernel, ATC_target_protein_kernel, ATC_side_effects_kernel, ATC_fingerprint_kernel])
atc_global_kernel = atc_matrix_combination([ATC_atc_columns_kernel ATC_atc_SM_kernel, ATC_target_protein_kernel, ATC_side_effects_kernel, ATC_fingerprint_kernel])
atc_global_kernel = atc_matrix_combination([ATC_atc_columns_kernel, ATC_atc_probabilistic_kernelATC_target_protein_kernel, ATC_side_effects_kernel, ATC_fingerprint_kernel])
atc_global_kernel = atc_matrix_combination([ATC_atc_columns_kernel, ATC_atc_probabilistic_kernel, ATC_atc_SM_kernelATC_side_effects_kernel, ATC_fingerprint_kernel])
atc_global_kernel = atc_matrix_combination([ATC_atc_columns_kernel, ATC_atc_probabilistic_kernel, ATC_atc_SM_kernel, ATC_target_protein_kernel, ATC_fingerprint_kernel])
atc_global_kernel = atc_matrix_combination([ATC_atc_columns_kernel, ATC_atc_probabilistic_kernel, ATC_atc_SM_kernel, ATC_target_protein_kernel, ATC_side_effects_kernel])
```
## The models of single ATC code kernel
```python
atc_global_kernel = atc_matrix_combination([ATC_atc_columns_kernel])
atc_global_kernel = atc_matrix_combination([ATC_atc_probabilistic_kernel])
atc_global_kernel = atc_matrix_combination([ATC_atc_SM_kernel])
atc_global_kernel = atc_matrix_combination([ATC_target_protein_kernel])
atc_global_kernel = atc_matrix_combination([ATC_fingerprint_kernel])
```
## The models of single fingerprint kernel
```python
WK = WKNKN(drug_global_kernel, atc_global_kernel, train_data.values, self.omega)
```
```python
WK = WKNKN(drug_fingerprint_kernel, ATC_fingerprint_kernel, train_data.values, self.omega)
```
# Result
The PR curves and ROC curves predicted by our model on the dataset are shown below:
1. The PR curves 
2. The ROC curves



