for i in OCT_age_prediction_all
do
# define path
dataDir="/data/PheWAS_baselineVisit"
# results File
CatDir="${dataDir}/results/"
mkdir -p $CatDir
codeDir='/software/PHESANT-master'
# expouse file and  trait of interests name
expFile="${dataDir}/${i}_PheWASinput.csv"
trait_name="value"
 
# ---Step1 Running a phenome scan in UK Biobank-----
## pheno data File
phenodataFile="${dataDir}/withbaseline_data_PheWASinput.csv"
confounderFile="${dataDir}/cova_data_more_withPC_all2visits_PheWASinput.csv"

# run PHESANT
cd $codeDir/WAS
Rscript phenomeScan.r \
--phenofile="$phenodataFile" \
--traitofinterestfile="$expFile" \
--confounderfile="$confounderFile" \
--variablelistfile="${codeDir}/variable-info/outcome-info.tsv" \
--datacodingfile="${codeDir}/variable-info/data-coding-ordinal-info.txt" \
--traitofinterest="$trait_name" \
--resDir="${CatDir}"  \
--userId="userId"  \
--sensitivity \
--genetic=FALSE 

# ---Step 2 Post-processing of results----
codeDir2="/software/PHESANT-master/resultsProcessing/"
cd $codeDir2
chmod 755 mainCombineResults.r
Rscript mainCombineResults.r \
--resDir="${CatDir}"  \
--variablelistfile="${codeDir}/variable-info/outcome-info.tsv" \


# ---change the foldername---
cd /data/PheWAS_baselineVisit
mkdir $i
mv results $i
done
