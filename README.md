# MetaRib V1.0
## Manual 
MeteRib is an iterative tool for ribosomal gene reconstruction from whole RNA meta-transcriptomic data.   
         
__Dependencies:__ .     
* EMIRGE: https://github.com/csmiller/EMIRGE
* BBMAP tools: https://jgi.doe.gov/data-and-tools/bbtools/  
* Python: python2.7, pandas, seaborn.   


__Usage:__     

Assume your current folder is: */project*. Then the following setup would be appropriate to run MetaRib.    

Step1: Data preparation.    

MetaRib supports multiple samples. If you data folder is: */project/data*. Your data structure need to be like this:   

\________________________________   
Content of samples.list.txt files:    
*SAMPLE_1*   
*SAMPLE_2*   
\________________________________    
Sample files:   
*/project/data/SAMPLE_1.1.fq*      
*/project/data/SAMPLE_1.2.fq*          
*/project/data/SAMPLE_2.1.fq*         
*/project/data/SAMPLE_2.1.fq*             
    
Step2: Set configure file.    

[BASE] .    
PROJECT_DIR :  # Project Dir.   
DATA_DIR : # Raw data storage dir.     
SAMPLING_NUM : # Subsampling reads number in each iteration.     
THREAD : # Number of cores used in mapping and assemble.        

[EMIRGE] .      
EM_PATH : ~/bin/EMIRGE/emirge_amplicon.py
EM_PARA : --phred33 -l 125 -i 250 -s 50 -a 20 -n 20
EM_REF : /export/jonassenfs/yaxinx/meta_rna/final_test/MetaRib/references/s1.ref.fa
EM_BT : /export/jonassenfs/yaxinx/meta_rna/final_test/MetaRib/references/s1.ref.btindex

[BBMAP]
BM_PATH: /export/jonassenfs/yaxinx/local/software/bbmap
MAP_PARA : minid=0.96 maxindel=1 minhits=2 idfilter=0.98
CLS_PARA : fo=t ow=t c=t mcs=1 e=5 mid=99
