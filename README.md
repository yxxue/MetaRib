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

MetaRib supports multiple samples. If you data folder is: */project/data*. Your data structure need to be like this:   

\________________________________   
Content of *samples.list.txt* files:    
*SAMPLE_1*   
*SAMPLE_2*   
\________________________________    
Sample files:   
*/project/data/SAMPLE_1.1.fq*      
*/project/data/SAMPLE_1.2.fq*          
*/project/data/SAMPLE_2.1.fq*         
*/project/data/SAMPLE_2.1.fq*             
    
***Step2: Setup configure file.***       
Content of *MetaRib.cfg* need to have the folowing information:   

[BASE]      
PROJECT_DIR : /project/   # your project folder   
DATA_DIR : /project/data/  # your data folder     
SAMPLING_NUM : 100000   # Subsampling reads number in each iteration      
THREAD : 20 # Number of cores used in MetaRib
SEED : 100 # SEED number for random subsampling         

[EMIRGE]         
EM_PATH : ~/bin/EMIRGE/emirge_amplicon.py # EMIRGE path
EM_PARA : --phred33 -l 125 -i 250 -s 50 -a 20 -n 20 # EMIRGE parameters 
EM_REF : ~/bin/EMIRGE/references/silva.v123.fa # EMIRGE references     
EM_BT : ~/bin/EMIRGE/references/silva.v123.fa.btindex # EMIRGE reference index

[BBMAP]   
BM_PATH: ~/bin/BBTOOL/ # BBTOOLS path     
MAP_PARA : minid=0.96 maxindel=1 minhits=2 idfilter=0.98 # BBTOOLS mapping parameters    
CLS_PARA : fo=t ow=t c=t mcs=1 e=5 mid=99 # BBTOOLS cluster parameters          



