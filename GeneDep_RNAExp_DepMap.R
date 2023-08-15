library("dplyr")
library("tibble")
library("gridExtra")
library("stringr")
library(tidyverse)
library(magrittr)
library(annotables)
library(data.table)

# the root directory containing the script and the proteome and rna measurements file
base = "C:/Users/Documents/DepMap/gene_dependency/"


########################
#####1. Data Download###
########################

#read RNA expression dataframe. I used the most recent DEPMAP RNA Expression file 23Q2 downloaded from https://depmap.org/portal/download/all/
target_expression= read_csv(paste0(base,"TPMLogp1_23Q2.csv"))
gene_dep = read_csv(paste0(base,"CRISPRGeneDependency.csv"))


##########################################################################################
#############RNA Expression & Gene Dependency DATA CLEAN UP#############################
########################################################################################

###2. Convert Gene Dependency data to long format to make it easier to use tidyverse functions
keycol <- "gene"
valuecol <- "GeneDep"
gathercols <- colnames(gene_dep)[2:ncol(gene_dep)]

gene_dep = gather(gene_dep, keycol, valuecol, gathercols)

gene_dep2 <-  gene_dep %>% 
  rename ("gene"= keycol) %>%
  rename ("GeneDep"= valuecol)
gene_dep3 = data.table(gene_dep2)



###3. Convert RNA expression data to long format to make it easier to use tidyverse functions
keycol <- "gene"
valuecol <- "TPMlogp1"
gathercols <- colnames(target_expression)[2:ncol(target_expression)]

target_expression = gather(target_expression, keycol, valuecol, gathercols)

target_expression2 <- target_expression %>% 
  rename("TPMlogp1" = valuecol)%>%
  rename("gene" = keycol)





###4.combine crispr gene dependency data to gene expression
target_expression3 <- data.table(target_expression2)

target_expression3[gene_dep3, GeneDep := i.GeneDep, on = .(ModelID, gene)]

###5. Gene symbol in table is written as Gene (Num)- remove gene number and keep only the gene symbol#################
#target_expression2 <-  separate(target_expression, col=keycol,into=c("gene_symbol", "number"),sep=" \\(")

target_expression4 <- data.table(target_expression3)
# Split the `keycol` into two separate columns `gene_symbol` and `number`
target_expression4[, c("gene_symbol", "number") := tstrsplit(gene, "\\s\\(")]

# Remove the closing parenthesis ")" from the `number` column
target_expression4[, number := gsub("\\)", "", number)]

# delete Keycol and number
target_expression4[, c("gene", "number") := NULL]

############################
#data table includes modelID instead of cell line name. Download Model Table to add cell line name. Download Model table from https://depmap.org/portal/download/all/

###6. add cell line name


cell_line_name= 
  read.csv(paste0(base,"Model.csv")) %>%
  dplyr::select(ModelID,  StrippedCellLineName,OncotreePrimaryDisease, OncotreeLineage)
colnames(cell_line_name) = c("ModelID", "Stripped_cell_line_name", "cancer_type_final", "Tissue")

##7. combine both tables to include actual cell line names 

data =
  left_join(cell_line_name, target_expression4, by = "ModelID", copy = FALSE,
            keep = NULL,
            na_matches = c("na", "never"),
            multiple = "all",
            unmatched = "drop",
            relationship = NULL) 
data2 <- data %>%
  select(-ModelID)



#8.write table

write.table(data2, file="C:/Users/Documents/DepMap/gene_dependency_depmap/downloaded_data_depmap2_genedep.tsv",row.names=FALSE, sep="\t")


###################################################################################################################################
##################################RNA Expression and Gene Dependency Determination#################################################
###################################################################################################################################

# top X cell lines to find for this protein / gene- Currently DepMap has data for 1864 cancer cell lines
top_n = 2000

# most interesting cell lines we have #change to include other cell lines that you may be interested in
top_cells = c("HELA","HAP1", "THP1", "U2OS", "HCT116", "A549", "A431", "U87MG")


## provide a list of gene names here
## made a new file without header info from Mona's E3 ligase list
A_search = 
  read_tsv(paste0(base,"e3_ligase_uniprot.tsv")) %>%
  pull(`Gene Symbol`)


#################
#### outputs ####
#################

# output diretory containing the output files
# one file will be produced for each gene
out_dir = paste0(base,"output_rna_ma/")
dir.create(out_dir,showWarnings = F)

#one output file for all genes
##You HAVE TO reset the table every time you run the script, otherwise it will keep on adding to the old values

selectlines_uniprot_rna = data.frame()
allcell_data_rna = data.frame()



temp =
  dat %>%
  #filter(Study.x %in% studies) %>%
  dplyr::select(gene_symbol, Stripped_cell_line_name, TPMlogp1, GeneDep, cancer_type_final,Tissue) %>%
  rename(RNA_count=TPMlogp1)  %>%
  rename(cell_line=Stripped_cell_line_name)


##################################################
#### 1. FIND TOP RANKED CELL LINES IN RNA-SEQ ####
##################################################
#one output file for all genes A= "ERBB4"
##You HAVE TO reset the table every time you run the script, otherwise it will keep on adding to the old values


for(A in A_search){
  
  # print the gene being processed
  cat(A,"\n")
  
  temp2 =
    temp %>%
    filter(gene_symbol==!!A)
  
  # skip if nothing returned
  if(nrow(temp2)==0){
    next
  }
  
  # top N by RNA count
  top_rna =
    temp2 %>%
    arrange(desc(RNA_count)) %>%
    filter(RNA_count >=3 )
  #slice(1:top_n)
  
  
  # save this table if entries
  # skip if nothing returned
  if(nrow(top_rna)>=1){
    #df_temp_rna=c("RNA_logTPMplus1_",A,"_Nusinow_",top_n)
    allcell_data_rna <-rbind.data.frame(allcell_data_rna, top_rna)
    
  }
  
}


write.table(allcell_data_rna, file="C:/Users/Documents/DepMap/gene_dependency_depmap/allcell_data_rna.tsv",row.names=FALSE, sep="\t")




########################################################
##### 2. Find Expression in Preferred Cell lines###########
#########################################################

## find any cell lines that are in common in the lists
## rank them and then inner join them
## only if both studies returned results

# explore the cell lines we have that are interesting
select_cell_lines =
  temp %>%
  filter(cell_line %in% top_cells) %>%
  pull(cell_line) %>%
  unique()



for(A in A_search){
  
  # print the gene being processed
  cat(A,"\n")
  
  temp2 =
    temp %>%
    filter(gene_symbol==!!A)
  
  # skip if nothing returned
  if(nrow(temp2)==0){
    next
  }
  
  # top N by RNA count
  top_rna =
    temp2 %>%
    arrange(desc(RNA_count)) %>%
    filter(RNA_count >=3 )  %>%
    filter(cell_line %in% top_cells)  
  #slice(1:top_n)
  
  
  # save this table if entries
  # skip if nothing returned
  if(nrow(top_rna)>=1){
    #df_temp_rna=c("RNA_logTPMplus1_",A,"_Nusinow_",top_n)
    selectlines_uniprot_rna <-rbind.data.frame(selectlines_uniprot_rna, top_rna)
    
  }
  
}



write.table(selectlines_uniprot_rna, file="C:/Users/Documents/DepMap/gene_dependency_depmap/output_ma/uniprot_rna_selectlines.tsv",row.names=FALSE, sep="\t")



