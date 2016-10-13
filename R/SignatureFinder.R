#' Find gene signature or/and built clasifier.
#'
#' @param table Data frame with gene in columns and patients or sample
#' in rows
#' @param x Name of target column. (With name ofcancers)
#' @param k Number of select genes in first step. (In Kruskal-wallis
#' statistic). In default k=100.
#' @param m Number of important characteristic function in which our
#' group of gene are enriched. In default m=20.
#' @param p Treshold to chose important characteristic function in which our
#' group of gene are enriched. In default p=0.05.
#' @param signature.method Method to calculate the signature. Possible
#' names: all.all.med, all.min.max, all.5.med, m.all.max, p.all.med.
#' In default signature.method="p.all.med".
#' @param out The parameter describing what we would like to receive
#' as a result of algorithm. Possible names: signature, SVM.method,
#' random.forest.method. In default out="signature".
#'
#'
#' @return Signature - vector with genes name. Clasifier of method SVM or
#' Random Forest.
#'
#' @export

signature.algorithm <- function(table, x, k=100,m=20,p=0.05,
                                signature.method="p.all.med",
                                out="signature"){

   if(!is.data.frame(table)){
      stop(paste("Your data is not a data.frame"))
   }
   if(length(which(colnames(table)==x)) == 0 ){
      stop(paste("There is no column called",x))
   }
   if(k>dim(table)[2]){
      stop(paste("Your k is to big"))
   }
   if(p>1 |p<0){
      stop(paste("p mast be a value in [0,1]"))
   }
   if(!(signature.method == "all.all.med" |
        signature.method == "all.min.max" |
        signature.method == "all.5.med" |
        signature.method == "m.all.max" |
        signature.method == "p.all.med" )){
      stop(paste("There is no signature method called",signature.method))
   }
   if(!(out == "signature" |
        out == "SVM.method" |
        out == "random.forest.method" )){
      stop(paste("There is no output called",out))
   }
   convert.data <- function(table, x){
      index_of_target_wariable <- which(colnames(table)==x)
      table <- cbind(table[,-index_of_target_wariable],table[,index_of_target_wariable])
      colnames(table)[dim(table)[2]] <- "target"
      name_of_column <- colnames(table)
      output_data <- data.matrix(table, rownames.force = NA)
      output_data <- output_data[!rowSums(table[,-dim(table)[2]])==0,]
      index_of_non_zeros_column <- which(!(colSums(table[,-dim(table)[2]])==0))
      output_data <- output_data[,c(index_of_non_zeros_column,dim(table)[2])]
      return(output_data)
   }


  
   matrix.gene.GOfunction <- function(table){
      names_of_gene <- colnames(dane)[-dim(dane)[2]]
      gene_function <- select(org.Hs.eg.db,
                              keys = names_of_gene,
                              keytype = "SYMBOL",
                              columns = "GO")
      BP <- which(gene_function[,4]=="BP")
      gene_function <- gene_function[BP,1:2]
      gene_function <- unique(gene_function)
      colnames(gene_function) <- c("id_gene","id_GO")
      return(gene_function)
   }

   KW.ranking <- function(table){
      K_W <- numeric()
      n <- dim(table)[2]-1
      for(i in 1:n){
         K_W <- c(K_W,kruskal.test(target ~ table[,i] ,
                                   data = table)$p.value)
      }
      names_of_gene <- colnames(table)[-dim(table)[2]]
      KW_data <- data.frame(names_of_gene,K_W)
      rownames(KW_data) <- NULL
      colnames(KW_data) <- c("id_gene","KW")
      return(KW_data)

   }

   enriched.GO <- function(table, KW_table, k){
      KW_order <- KW_table[order(KW_table$KW),]
      greatest_genes <- as.vector(KW_order[1:k,1])
      names_of_gene <- colnames(table)[-dim(table)[2]]
      geneList <- factor(as.integer(names_of_gene %in% greatest_genes))
      names(geneList) <- names_of_gene
      BPterms <- ls(GOBPTerm)
      data_GO <- new("topGOdata",
                     ontology="BP",
                     allGenes=geneList,
                     annot=annFUN.org,
                     mapping = "org.Hs.eg.db",
                     ID = "symbol")
      data_topGOresult <- runTest(data_GO
                                  ,algorithm="classic"
                                  ,statistic = "fisher")
      char_function <- GenTable(data_GO,
                                classicFisher = data_topGOresult)[,c(1,6)]
      colnames(char_function) <- c("id_GO","fisher")
      return(char_function)
   }

   merge.all.data <- function(KW_table,gene_function,char_function){
      data_base <- merge(x = KW_table,
                         y = gene_function,
                         by = "id_gene",
                         all.x = TRUE)
      data_base <- merge(x = data_base,
                         y = char_function,
                         by = "id_GO",
                         all.x = TRUE)
      data_base$fisher[is.na(data_base$fisher)] <- 1
      data_base$fisher <- as.numeric(data_base$fisher)
      multiply <- data_base$KW*data_base$fisher
      data_base <- cbind(data_base, multiply)
      data_base <- data_base[order(data_base$id_gene,data_base$fisher),]
      return(data_base)
   }


   merge.par.m <- function(KW_table,gene_function,char_function,m=20){
      char_function <- char_function[order(char_function$fisher),]
      char_fun_m <- char_function[1:m,]
      data_base_m <- merge(x = char_fun_m,
                           y = gene_function,
                           by = "id_GO",
                           all.x = TRUE)
      data_base_m <- data_base_m[-c(which(is.na(data_base_m$id_gene)==1)),]
      data_base_m <- merge(x = KW_table,
                           y = data_base_m,
                           by = "id_gene",
                           all.x = TRUE)
      data_base_m$fisher[is.na(data_base_m$fisher)] <- 1
      data_base_m$fisher <- as.numeric(data_base_m$fisher)
      multiply <- data_base_m$KW*data_base_m$fisher
      data_base_m <- cbind(data_base_m, multiply)
      data_base_m <- data_base_m[order(data_base_m$id_gene,data_base_m$fisher),]
      return(data_base_m)
   }




   merge.par.p <- function(KW_table,gene_function,char_function,p=0.05){
      char_function <- char_function[order(char_function$fisher),]
      char_fun_p <- char_function[which(char_function$fisher<p),]
      data_base_p <- merge(x = char_fun_p,
                           y = gene_function,
                           by = "id_GO",
                           all.x = TRUE)
      data_base_p <- data_base_p[-c(which(is.na(data_base_p$id_gene)==1)),]
      data_base_p <- merge(x = KW_table,
                           y = data_base_p,
                           by = "id_gene",
                           all.x = TRUE)
      data_base_p$fisher[is.na(data_base_p$fisher)] <- 1
      data_base_p$fisher <- as.numeric(data_base_p$fisher)
      multiply <- data_base_p$KW*data_base_p$fisher
      data_base_p <- cbind(data_base_p, multiply)
      data_base_p <- data_base_p[order(data_base_p$id_gene,data_base_p$fisher),]
      return(data_base_p)
   }


   all.all.med <- function(table_all,KW_table,k=100){
      score <- aggregate(. ~ id_gene, table_all[c(2,5)], prod)
      colnames(score) <- c("id_gene","score")
      all_all_med <- merge(x = KW_table,
                           y = score,
                           by = "id_gene",
                           all.x = TRUE)
      rank_score <- 1:dim(all_all_med)[1]
      rank_KW <- 1:dim(all_all_med)[1]
      all_all_med <- cbind(all_all_med[order(all_all_med$score,all_all_med$id_gene),],rank_score)
      all_all_med <- cbind(all_all_med[order(all_all_med$KW,all_all_med$id_gene),],rank_KW)
      all_all_med <- all_all_med[order(all_all_med$id_gene),]
      med <- median(all_all_med$score[which(all_all_med$rank_KW<=k)])
      signature_all_all_med <- all_all_med$id_gene[which(all_all_med$score<=med | all_all_med$rank_KW<=k )]
      return(signature_all_all_med)

   }



   all.min.max <- function(table_all,KW_table,k=100){
      score <- aggregate(. ~ id_gene, table_all[c(2,5)], min)
      colnames(score) <- c("id_gene","score")
      all_min_max <- merge(x = KW_table,
                           y = score,
                           by = "id_gene",
                           all.x = TRUE)
      rank_score <- 1:dim(all_min_max)[1]
      rank_KW <- 1:dim(all_min_max)[1]
      all_min_max <- cbind(all_min_max[order(all_min_max$score,all_min_max$id_gene),],rank_score)
      all_min_max <- cbind(all_min_max[order(all_min_max$KW,all_min_max$id_gene),],rank_KW)
      all_min_max <- all_min_max[order(all_min_max$id_gene),]
      m <- max(all_min_max$score[which(all_min_max$rank_KW<=k)])
      signature_all_min_max <- all_min_max$id_gene[which(all_min_max$score<=m )]
      return(signature_all_min_max)

   }


   all.5.med <- function(table_all,KW_table,k=100){
      first_5 <- by(table_all,table_all["id_gene"],head,n=5)
      data_base_5 <- Reduce(rbind, first_5)
      score <- aggregate(. ~ id_gene, data_base_5[c(2,5)], prod)
      colnames(score) <- c("id_gene","score")
      all_5_med <- merge(x = KW_table,
                         y = score,
                         by = "id_gene",
                         all.x = TRUE)
      rank_score <- 1:dim(all_5_med)[1]
      rank_KW <- 1:dim(all_5_med)[1]
      all_5_med <- cbind(all_5_med[order(all_5_med$score,all_5_med$id_gene),],rank_score)
      all_5_med <- cbind(all_5_med[order(all_5_med$KW,all_5_med$id_gene),],rank_KW)
      all_5_med <- all_5_med[order(all_5_med$id_gene),]
      med <- median(all_5_med$score[which(all_5_med$rank_KW<=k)])
      signature_all_5_med <- all_5_med$id_gene[which(all_5_med$score<=med | all_5_med$rank_KW<=k )]
      return(signature_all_5_med)
   }



   m.all.max <- function(table_m,KW_table,k=100){
      score <- aggregate(. ~ id_gene, table_m[c(1,5)], prod)
      colnames(score) <- c("id_gene","score")
      m_all_max <- merge(x = KW_table,
                         y = score,
                         by = "id_gene",
                         all.x = TRUE)
      rank_score <- 1:dim(m_all_max)[1]
      rank_KW <- 1:dim(m_all_max)[1]
      m_all_max <- cbind(m_all_max[order(m_all_max$score,m_all_max$id_gene),],rank_score)
      m_all_max <- cbind(m_all_max[order(m_all_max$KW,m_all_max$id_gene),],rank_KW)
      m_all_max <- m_all_max[order(m_all_max$id_gene),]
      m <- max(m_all_max$score[which(m_all_max$rank_KW<=k)])
      signature_m_all_max <- m_all_max$id_gene[which(m_all_max$score<=m )]
      return(signature_m_all_max)
   }

   p.all.med <- function(table_p,KW_table,k=100){
      score <- aggregate(. ~ id_gene, table_p[c(1,5)], prod)
      colnames(score) <- c("id_gene","score")
      p_all_med <- merge(x = KW_table,
                         y = score,
                         by = "id_gene",
                         all.x = TRUE)
      rank_score <- 1:dim(p_all_med)[1]
      rank_KW <- 1:dim(p_all_med)[1]
      p_all_med <- cbind(p_all_med[order(p_all_med$score,p_all_med$id_gene),],rank_score)
      p_all_med <- cbind(p_all_med[order(p_all_med$KW,p_all_med$id_gene),],rank_KW)
      p_all_med <- p_all_med[order(p_all_med$id_gene),]
      med <- median(p_all_med$score[which(p_all_med$rank_KW<=k)])
      signature_p_all_med <- p_all_med$id_gene[which(p_all_med$score<=med | p_all_med$rank_KW<=k )]
      return(signature_p_all_med)
   }



   SVM.method <- function(table,signature){
      names_of_gene <- colnames(table)[-dim(table)[2]]
      index_of_gene <- which(names_of_gene %in% signature)
      y <- table[,dim(table)[2]]
      y <- as.factor(y)
      x <- as.matrix(table[,index_of_gene])
      svm_method = svm(x, y, cost = 10, cachesize=500,
                       scale=T, type="C-classification", kernel="linear" )
      return(list(svm_method, signature))

   }

   random.forest.method <- function(table,signature){
      names_of_gene <- colnames(table)[-dim(table)[2]]
      index_of_gene <- which(names_of_gene %in% signature)
      train_data <- table[,c(index_of_gene,dim(table)[2])]
      train_data <- data.frame(train_data)
      train_data$target <- as.factor(train_data$target)
      RF_method <- randomForest(target~.,
                                data=train_data,
                                importance=TRUE,
                                proximity=TRUE)
      return(list(RF_method, signature))

   }




   data <- convert.data(table,x)
   matrix_gene_function <-matrix.gene.GOfunction(data)
   KW_data <- KW.ranking(data)
   charakteristic_function <- enriched.GO(data,KW_data,k)
   if(signature.method=="all.all.med"){
      data_all <- merge.all.data(KW_data,matrix_gene_function,charakteristic_function)
      signature <- all.all.med(data_all,KW_data,k)
   }else if(signature.method=="all.min.max"){
      data_all <- merge.all.data(KW_data,matrix_gene_function,charakteristic_function)
      signature <- all.min.max(data_all,KW_data,k)
   }else if(signature.method=="all.5.med"){
      data_all <- merge.all.data(KW_data,matrix_gene_function,charakteristic_function)
      signature <- all.5.med(data_all,KW_data,k)
   } else if(signature.method=="m.all.max"){
      data_m <- merge.par.m(KW_data,matrix_gene_function,charakteristic_function,m)
      signature <- m.all.max(data_m,KW_data,k)
   } else if(signature.method=="p.all.med"){
      data_p <- merge.par.p(KW_data,matrix_gene_function,charakteristic_function,p)
      signature <- p.all.med(data_p,KW_data,k)
   } else {
      print("You choose wrong method to find signature: all.all.med,all.min.max,all.5.med,m.all.max,p.all.med are correct value")
   }



   if(out=="signature"){
      output <- signature
   } else if(out=="predictive.SVM"){
      predictive_method <- SVM.method(data,signature)
      output <- predictive_method
   } else if(out=="predictive.RandomForest"){
      predictive_method <- random.forest.method(data,signature)
      output <- predictive_method
   } else {
      print("You choose wrong output: signature, predictive.SVM, predictive.RandomForest are only correct value ")
   }

   return(output)

}
