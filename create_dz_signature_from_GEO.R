#this code is to create a disease signature
#samr and rankprod are used (or maybe bayesian); SAMR runs faster than rankprod, but usually gives a less number of genes.
#probe annotation can be downloaded from GEO 


set.seed(97)

#packages
library(GEOquery)
library(Biobase)
library(preprocessCore)
library("gplots")


#######################
#functions#####
# Very crude for now
is_logged <- function(sample_frame) {
  mean_sample_var <- mean(apply(sample_frame,2,var,na.rm=T),na.rm=T)
  ifelse(mean_sample_var > 10,FALSE,TRUE)
}

stripNumberPrefix <-function(v, prefix) {
  # Fix this to only work if there is a prefix!  Right now, will error!
  # Definitely need to fix this up!!
  
  #return(as.numeric(as.matrix(data.frame(strsplit(v, paste("^", prefix, sep="") ))[2,] )))
  #use GSUB instead?  
  return(as.matrix(gsub(paste("^", prefix, sep="" ), "", v)))
}

getEntrezMappingsFromGEO <- function(gpl_id){
  gpl <- getGEO(gpl_id, destdir=".")
  if (exists(paste0(gpl_id, "_mapped.csv")))  geneMappings = read.csv(paste0(gpl_id, "_mapped.csv"))
  if (gpl_id == "GPL4133"){
    geneMappings <- Table(gpl)[,c("ID","GENE","GENE_SYMBOL","GENE_SYMBOL")]   # the name of GPL4133 is messed up     
  }else if (gpl_id == "GPL10558"){
    geneMappings <- Table(gpl)[,c("ID","Entrez_Gene_ID","Symbol","Definition")]        
  }else if (gpl_id == "GPL10687"){
    geneMappings <- Table(gpl)[,c("ID","EntrezGeneID","GeneSymbol","GeneSymbol")]   
    # }else if (gpl_id == "GPL571" | gpl_id == "GPL3921"){ #this two gpls' annotation is not right from the website, need to re-processed, so we use from our server
    #geneMappings <- Table(gpl)[,c("ID","ENTREZ_GENE_ID","Gene Symbol","Gene Title")]
    #  geneMappings <- getEntrezMappings(gpl_id, con)
    #   geneMappings$Name <- geneMappings$Symbol # name does not return
  }else if (gpl_id == "GPL10687"){
    geneMappings <- Table(gpl)[,c("ID","EntrezGeneID","GeneSymbol","GeneSymbol")]
  }else if (gpl_id == "GPL7202"){
    geneMappings <- Table(gpl)[,c("ID","GENE","GENE_SYMBOL","GENE_NAME")]  
  }else if (gpl_id == "GPL3738"){
    geneMappings <- Table(gpl)[,c("ID","ENTREZ_GENE_ID","Gene Symbol","Gene Title")]    
  }else if (gpl_id == "GPL96" | gpl_id == "GPL1352"){
    geneMappings <- Table(gpl)[,c("ID","ENTREZ_GENE_ID","Gene Symbol",	"Gene Symbol")]
  }else if (gpl_id == "GPL8177"){
    geneMappings <- Table(gpl)[,c("ID","UGCluster","Symbol","Name")]
    gene2unigene <- getGene2Unigene(con)
    geneMappings <- merge(geneMappings, gene2unigene, by.x = "UGCluster", by.y="Unigene" )
    geneMappings <- subset(geneMappings, select = c("ID", "GeneID", "Symbol", "Name"))
  }else {
    #try to guess column names of entrez gene id and symbol name
    putative_fields <- grep("ENTREZ|GENE|Symbol", names(Table(gpl)), ignore.case=T)
    house_keeping_gene_ids <- c("23521", "7316", "567")
    house_keeping_gene_symbols <- c("RPL13A", "UBC", "B2M")
    gene_field_name <- ""
    symbol_field_name <- ""
    for (field in putative_fields){
      if (sum(toupper(Table(gpl)[,field]) %in% house_keeping_gene_ids) > 0 ){
        gene_field_name <- names(Table(gpl))[field]
      }
    }
    
    for (field in putative_fields){
      if (sum(Table(gpl)[,field] %in% house_keeping_gene_symbols) > 0 ){
        symbol_field_name <- names(Table(gpl))[field]
      }
    }     
    
    if (symbol_field_name != "" & gene_field_name != ""){
      geneMappings <- Table(gpl)[,c("ID", gene_field_name, symbol_field_name,symbol_field_name)]
    }
  }
  
  names(geneMappings) <- c("probe","GeneID","Symbol","NAME")
  geneMappings <- subset(geneMappings, !is.na(GeneID) & GeneID != "")
  dupProbes <- unique(geneMappings$probe[duplicated(geneMappings$probe)])
  geneMappings <- geneMappings[! geneMappings$probe %in% dupProbes,] #remove duplicated probes
  rownames(geneMappings) <- geneMappings$probe
  
  #some mapping files, GENEID composed by muliple gene ids, should split.. this may cause problems for the symbol and name
  if (class(geneMappings$GeneID) != "integer"){
    geneMappings$GeneID <- as.character(geneMappings$GeneID)
    geneMappings$GeneID  = sapply(geneMappings$GeneID, function(x) unlist(strsplit(x, " "))[1]) #not a perfect approach
    geneMappings$Symbol  = sapply(as.character(geneMappings$Symbol), function(x) unlist(strsplit(x, " "))[1]) #not a perfect approach
    
    'newGeneMappings <- data.frame()
    for (i in 1:nrow(geneMappings)){
      newGeneMappings <- rbind(newGeneMappings, data.frame(probe = geneMappings$probe[i], 
                                                          GeneID = unique(na.omit(as.numeric(unlist(strsplit(geneMappings$GeneID[i], "[^0-9]+"))))),
                                                          Symbol = geneMappings$Symbol[i], 
                                                          NAME = geneMappings$NAME[i]))
    }
    geneMappings <- newGeneMappings
    '
  }
  
  return(geneMappings)
}

######################
data <- (getGEO(GSE, destdir=".", GSEMatrix=TRUE)[[1]])
expr_matrix <- exprs(data)
if (length(grep("GSM", colnames(expr_matrix)[1] )) == 0 ) colnames(expr_matrix) = phenoData(data)@data$geo_accession
nn <- phenoData(data)@data

#if it is a regular expression, then run case selection. Otherwise, assuming case is a list of samples
if (length(case_reg) == 1){
case <- as.character(nn$geo_accession[grep(case_reg, nn$title)])
control <- as.character(nn$geo_accession[grep(control_reg, nn$title)])
}else{
  case = case_reg
  control = control_reg
}

#
case <- case[!case %in% control]

platform <- annotation(data)

######the annotations from GEO are twice bigger than those from AILLUN (which might be out of date)
egm <- getEntrezMappingsFromGEO(platform)
##################################

#creat comparison
comparison_frame=subset(expr_matrix,select=c(control,case))
sample_class=c(rep(0,length(control)),rep(1,length(case)))

#remove probes if half of values are missed
complete_probes <- apply(comparison_frame,1,function(row) {
  ctl_vals <- row[sample_class==0]
  dz_vals <- row[sample_class==1]
  missing_ctl <- length(ctl_vals[is.na(ctl_vals)]) / length(ctl_vals)
  missing_dz <- length(dz_vals[is.na(dz_vals)]) / length(dz_vals)
  ifelse(((missing_ctl > 0.5) || (missing_dz > 0.5)),FALSE,TRUE)
})
comparison_frame <- comparison_frame[complete_probes,]

# deal negative values
#colSummarizeMedian(as.matrix(comparison_frame))
for (i in 1:ncol(comparison_frame)){
  comparison_frame[,i]=comparison_frame[,i]+abs(min(comparison_frame[,i],na.rm=TRUE))+1 #transform negatives to postivies, as quantile normailzation will be applied, it would affect 
}

# deal log
values_logged <- is_logged(comparison_frame)
if (!values_logged) {
  comparison_frame <- log2(comparison_frame+0.001) # add 1 to avoid log2(0)
  values_logged <- TRUE
}

#normalization
#boxplot(comparison_frame)
comparison_frame1 <- normalize.quantiles(as.matrix(comparison_frame ))
row.names(comparison_frame1) <- row.names(comparison_frame)
colnames(comparison_frame1) <- colnames(comparison_frame)
comparison_frame <- comparison_frame1

#collapse probe id using entrez id
comparison_frame <- merge(egm[, c("probe", "GeneID")], data.frame(probe = rownames(comparison_frame), comparison_frame), by = "probe")
#remove probe column
comparison_frame <- comparison_frame[, -1]
comparison_frame <- aggregate(. ~ GeneID, comparison_frame, mean)
rownames(comparison_frame) <- comparison_frame$GeneID
comparison_frame <- comparison_frame[, -1]

if (method_id==2) {
  method <- 'RANKPROD_SINGLE'
} else if (method_id==3){
  method <- "SIGGENES_SAMR"
} else if (method_id==4){
  method <- "SAMR"
}

#probe names
genenames<-rownames(comparison_frame)

if (method == 'RANKPROD_SINGLE') {
  library(RankProd)
  
  # Evaluate the impact of na.rm = TRUE. Alex M seems to think it's OK
  RP_result <- RP(comparison_frame, sample_class,gene.names=genenames, num.perm = 100, logged = T, na.rm = TRUE, plot = FALSE, rand = 123)
  # Leave logged=FALSE because topGene() converts the fold-change incorrectly!
  siggenes <- topGene(RP_result,cutoff=0.05,method="pfp",logged=FALSE,gene.names=genenames)
  # Normalize the results across methods
  siggenes.result <- list()
  siggenes.result$UP <- siggenes$Table1[,3:5]
  colnames(siggenes.result$UP) <- c("fold.change","q.value","p.value")
  siggenes.result$DOWN <- siggenes$Table2[,3:5]
  colnames(siggenes.result$DOWN) <- c("fold.change","q.value","p.value")
  # Since RANKPROD does the goofy condition 1 / condition 2, inverse the fold-change and convert to log if not logged. 
  if (values_logged) {
    siggenes.result$UP[,"fold.change"] <- -siggenes.result$UP[,"fold.change"]
    siggenes.result$DOWN[,"fold.change"] <- -siggenes.result$DOWN[,"fold.change"]
  } else {
    siggenes.result$UP[,"fold.change"] <- log2(1/siggenes.result$UP[,"fold.change"])
    siggenes.result$DOWN[,"fold.change"] <- log2(1/siggenes.result$DOWN[,"fold.change"])
  }
  siggenes.result$UP=data.frame(siggenes.result$UP,probe=rownames(siggenes.result$UP))
  siggenes.result$DOWN=data.frame(siggenes.result$DOWN,probe=rownames(siggenes.result$DOWN))
  
  
} else if (method == 'SIGGENES_SAMR') {
  # Using only the SAM implementation in SIGGENES, not EBAM
  library(siggenes)
  
  SAM_result <- sam(comparison_frame,sample_class,rand=123,R.unlog=T,gene.names=genenames) #q.version=1
  delta.table <- rbind(c(0,0,0),findDelta(SAM_result,fdr=0.9)) #fdr=0.9, too loose
  siggenes.table <- summary(SAM_result, delta.table[  dim(delta.table)[1],  1] );
  siggenes <- siggenes.table@mat.sig
  
  if( nrow(siggenes) > 0 ) {
    siggenes.result <- list()
    siggenes.result$UP <- siggenes[siggenes$d.value > 0,c("R.fold","q.value","rawp")]
    colnames(siggenes.result$UP) <- c("fold.change","q.value","p.value")
    siggenes.result$DOWN <- siggenes[siggenes$d.value < 0,c("R.fold","q.value","rawp")]
    colnames(siggenes.result$DOWN) <- c("fold.change","q.value","p.value")
    
    
    siggenes.result$UP=data.frame(siggenes.result$UP,probe=rownames(siggenes.result$UP))
    siggenes.result$DOWN=data.frame(siggenes.result$DOWN,probe=rownames(siggenes.result$DOWN))
    
  }
  
} else if (method == 'SAMR') {
  library(samr)
  # class labels for SAMR are 1 or 2
  input_data <- list(
    x=data.matrix(comparison_frame),
    y=sample_class+1,
    geneid=rownames(comparison_frame),
    genenames=rownames(comparison_frame),
    logged2=T
  )
  
  samr.obj <- samr(input_data,resp.type="Two class unpaired",testStatistic="standard",nperms=100)
  delta.table <- samr.compute.delta.table(samr.obj)
  delta.index <- which.max( delta.table[,5] < 0.4 )
  delta=delta.table[delta.index,1]
  
  #replace data with input_data
  siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, input_data, delta.table)
  
  siggenes.result <- list()
  
  if (!is.null(siggenes.table$genes.up)>0){
    UP <- as.data.frame(subset(siggenes.table$genes.up,select=c("Gene ID","Fold Change","q-value(%)")))
    colnames(UP) <- c("probe","fold.change","q.value")
    UP$p.value  <- NA    
    UP$probe <- as.character(UP$probe)
    UP$fold.change <- as.double(as.character(UP$fold.change))
    UP$q.value <- as.double(as.character(UP$q.value))
    siggenes.result$UP <- UP
  }
  if (!is.null(siggenes.table$genes.lo)){
    DOWN <- as.data.frame(subset(siggenes.table$genes.lo,select=c("Gene ID","Fold Change","q-value(%)")))
    colnames(DOWN) <- c("probe","fold.change","q.value")
    DOWN$p.value  <- NA    
    DOWN$probe <- as.character(DOWN$probe)
    DOWN$fold.change <- as.double(as.character(DOWN$fold.change))
    DOWN$q.value <- as.double(as.character(DOWN$q.value))
    siggenes.result$DOWN <- DOWN
  }
  
  
} 

#annotate probe by merging with GPL
genes.up <- merge(unique(egm[, c("GeneID", "Symbol")]), siggenes.result$UP, by.x ="GeneID", by.y = "probe", sort=F, all.x=T)
genes.down <- merge( unique(egm[, c("GeneID", "Symbol")]), siggenes.result$DOWN, by.x="GeneID",by.y = "probe", sort=F, all.x=T)

#removed missed match probes
genes.up <- genes.up[!is.na(genes.up$GeneID) & genes.up$GeneID != "", ]
genes.down <- genes.down[!is.na(genes.down$GeneID) & genes.down$GeneID != "", ]

genes.up[,"up_down"] <- "up"
genes.down[,"up_down"] <- "down"

genes.down$fold.change <- -1 * 1/genes.down$fold.change
genes.sig.up <- subset(genes.up, q.value < q_thresh & abs(fold.change) > fold_change , select=c("GeneID","Symbol","fold.change","q.value","p.value","up_down"))
genes.sig.down <- subset(genes.down, q.value < q_thresh & abs(fold.change) > fold_change, select=c("GeneID","Symbol","fold.change","q.value","p.value","up_down"))

#write disease signature
genes.sig.up <- genes.sig.up[order(-genes.sig.up$fold.change),]
genes.sig.down <- genes.sig.down[order(genes.sig.down$fold.change),]

dz_sig_all <- rbind(genes.sig.up,genes.sig.down)
colnames(dz_sig_all) <- c("GeneID","Symbol","value","q.value","p.value","up_down")
write.table(dz_sig_all,dz_sig_all.file,sep="\t",quote=F,row.names=F,col.names=T)

#mainly for cmap prediction
genes.sig.up <- genes.sig.up[1:min(nrow(genes.sig.up),150),]
genes.sig.down <- genes.sig.down[1:min(nrow(genes.sig.down),150),]
dz_sig <- rbind(genes.sig.up,genes.sig.down)

colnames(dz_sig) <- c("GeneID","Symbol","value","q.value","p.value","up_down")
write.table(dz_sig,dz_sig.file,sep="\t",quote=F,row.names=F,col.names=T)

#################################
#visualize disease signatures
annotation <- data.frame(type = c(rep("case", length(case)), rep("control", length(control))))
rownames(annotation) <- c(case, control)
annotation$type <- as.factor(annotation$type)
annotation <- subset(annotation, select=c("type"))

Var1        <- c("lightblue", "green")
names(Var1) <- c("case", "control")
anno_colors <- list(type = Var1)

my.cols <- greenred(100) 
comparison_frame_subset <- comparison_frame[rownames(comparison_frame) %in% dz_sig$GeneID, ]
pheatmap(t(scale(t(comparison_frame_subset))), col = my.cols, annotation = annotation,  annotation_colors = anno_colors,
         show_colnames=F, legend=T, show_rownames=F, filename=paste(disease, "/dz_sig_validation.pdf", sep="")
)

dz_expr = comparison_frame
save(dz_expr,annotation, my.cols, anno_colors , file=paste(disease, "/dz_expr.RData", sep=""))



