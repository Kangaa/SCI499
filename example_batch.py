import re,sys,os,time
models = range(1,9+1) 
cities = range(1,4+1)
for m in models:
 for c in cities:
  outFile = open('getRefTables-'+str(m)+'-'+str(c)+'.sh','w')
  print >> outFile, "#!/bin/bash"
  print >> outFile, "#SBATCH -p batch        	     "   
  print >> outFile, "#SBATCH -N 1               	     "        
  print >> outFile, "#SBATCH -n 1            	     "             
  print >> outFile, "#SBATCH --time=12:00:00    	     "                   
  print >> outFile, "#SBATCH --mem=4GB         	     "                  
  #if month == 12:
  print >> outFile, "#SBATCH --mail-type=BEGIN         "                      
  print >> outFile, "#SBATCH --mail-type=END           "                          
  print >> outFile, "#SBATCH --mail-type=FAIL              "                    
  print >> outFile, "#SBATCH --mail-user=robert.cope@adelaide.edu.au"            
  print >> outFile, "module load R/3.3.0-foss-2016uofa     "  
  print >> outFile, "R --no-save --args %d %d < PH_03_ref_tables.R > ref_tables_%d_%d.out" %(m,c,m,c)	                              
  outFile.close()
  os.system('sbatch getRefTables-'+str(m)+'-'+str(c)+'.sh')
  time.sleep(1)

