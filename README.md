# MetaRib V1.0
## Manual 
MeteRib is an iterative tool for ribosomal gene reconstruction from whole RNA meta-transcriptomic data.   
         
__Dependencies:__ .     
* EMIRGE: https://github.com/csmiller/EMIRGE
* BBtools: https://jgi.doe.gov/data-and-tools/bbtools/  
* Python: python2.7, pandas, seaborn.   


__Usage:__     

Assume your current folder is: */project*. Then the following setup would be appropriate to run MetaRib.    

***Step1: Data preparation.***       
MetaRib supports multiple samples. If you data folder is: */project/data*, merge all samples reads to generate all.1.fq and all.2.fq (cat SAMPLE*.1.fq > all.1.fq, cat SAMPLE*.2.fq > all.2.fq), your data folder structure need be like this:    

\________________________________   
Content of *samples.list.txt* files:    
*SAMPLE_1*   
*SAMPLE_2*   
\________________________________    
Data folder files:  
*/project/data/samples.list.txt*   
*/project/data/SAMPLE_1.1.fq*      
*/project/data/SAMPLE_1.2.fq*          
*/project/data/SAMPLE_2.1.fq*         
*/project/data/SAMPLE_2.1.fq*             
*/project/data/all.1.fq*   
*/project/data/all.2.fq*     

***Step2: Setup configure file.***       
Content of *MetaRib.cfg* need to have the folowing information:   
```
[BASE]    
# your project folder       
PROJECT_DIR : /project/   
# your data folder        
DATA_DIR : /project/data/        
# Subsampling reads number in each iteration        
SAMPLING_NUM : 100000     
# Number of cores used in MetaRib     
THREAD : 20   

[EMIRGE]     
# EMIRGE path       
EM_PATH : ~/bin/EMIRGE/emirge_amplicon.py  
# EMIRGE parameters   
EM_PARA : --phred33 -l 125 -i 250 -s 50 -a 20 -n 20  
# EMIRGE references   
EM_REF : ~/bin/EMIRGE/references/silva.v123.fa  
# EMIRGE reference index   
EM_BT : ~/bin/EMIRGE/references/silva.v123.fa.btindex   

[BBTOOL]   
# BBTOOLS path   
BB_PATH: ~/bin/BBTOOL/   
# BBTOOLS mapping parameters   
MAP_PARA : minid=0.96 maxindel=1 minhits=2 idfilter=0.98   
# BBTOOLS cluster parameters   
CLS_PARA : fo=t ow=t c=t mcs=1 e=5 mid=99   

[FILTER]    
# Minimium averge coverage in filter process    
MIN_COV : 2   
# Minimium coverge percent in filter process   
MIN_PER : 80    
```

***Step3: Run MetaRib.***  
MetaRib is developed in Python2.7 as EMIRGE is impletmented in Python2.7. After you prepare the data and configure file, run the code like this:    
```python
python2 run_MetaRib.py MetaRib.cfg
```
***Step4: Output.***       
All output will be stored at */project/MetaRib/*.    
*/project/MetaRib/Iteration/* includes running result of each iteration.   
*/project/MetaRib/Abundance/* includes a final constructed fasta *all.dedup.filtered.fasta*, an estimated abundance table file *all.dedup.filtered.est.ab.txt* and a heatmap distribution of top 20 most abudant contigs across samples *top20.heatmap.pdf*.     




