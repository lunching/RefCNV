GB_CNVs_EL<- function(case = coverage_all_cases, ref = coverage_all, CNV_q = .05, avg_cov_rm = 30, TMR = c(143156910, 143150524, 164584547, 319290992, 159630660, 152358720, 333807880, 158395099, 151853258), TMR_new = c(125576640, 117235819,154939360,130158617,115061086,141652146,140649220), gene_map = gene_map){
	 ##	LOOCV	##
	 nr<- nrow(ref[[1]])
	 ns<- length(filenames)
	 count_region<- 1
	 MSRs_loocv<- matrix(NA, nrow = nr, ncol = ns)
	 Reg_model_ref<- list()
	 Genes_all<- c()
		for(i in 1: nr)
		{
		 y_tmp<- sapply(1: ns, function(k) ref[[k]][i,4])
			if(mean(y_tmp) <= avg_cov_rm | sd(y_tmp)==0){}
			else
			{
			 Reg_model_ref[[count_region]]<- lm(y_tmp~TMR)
				for(j in 1: ns)
				{
				 TMR_tmp<- TMR[-j]
	 			 y_tmp_tmp<- y_tmp[-j]
				 fit_tmp<- lm(y_tmp_tmp ~ TMR_tmp)
				 TMR_testing_tmp<- data.frame(TMR_tmp = TMR[j])
				 PC_95CI_PDX_tmp<- predict(fit_tmp, TMR_testing_tmp, interval = "prediction",level = .99, se.fit=T)
				 coverage_tmp<- y_tmp[j]
				 std_residuals_new_tmp<- (coverage_tmp - PC_95CI_PDX_tmp$fit[,1])/(PC_95CI_PDX_tmp$residual.scale*sqrt(1 + (PC_95CI_PDX_tmp$se.fit/PC_95CI_PDX_tmp$residual.scale)^2))
				 MSRs_loocv[i, j]<- std_residuals_new_tmp 
				}
			 count_region<- count_region + 1
			}
		  chr_tmp<- sub("chr", "", ref[[1]][i,1])
		  start_tmp<- ref[[1]][i,2]
		  end_tmp<- ref[[1]][i,3]
		  gene_map_tmp<- gene_map[which(gene_map[,"Chromosome"]==chr_tmp), ]
		  region_tmp<- which(start_tmp<gene_map_tmp[,"Start"]& end_tmp>gene_map_tmp[,"End"]|start_tmp>gene_map_tmp[,"Start"]& end_tmp<gene_map_tmp[,"End"]|start_tmp<gene_map_tmp[,"Start"]& end_tmp>gene_map_tmp[,"Start"]|start_tmp<gene_map_tmp[,"End"]& end_tmp>gene_map_tmp[,"End"])
		  gene_tmp<- gene_map_tmp[region_tmp, ]
		 	if(length(nrow(gene_tmp))>0)
		 	{Genes_all[[i]]<- paste(gene_tmp[,"Gene_Symbol"], collapse=", ")}
		 	else
		 	{Genes_all[[i]]<- "NA"}
		}
	 ##	Summarize into gene-based	##
	 MSRs_loocv2<- MSRs_loocv[!is.na(MSRs_loocv[,1]), ]
	 na_index<- which(is.na(MSRs_loocv[,1]))
	 Genes_all2<- Genes_all[-na_index]
	 MSRs_loocv2<- cbind(MSRs_loocv2, Genes_all2)
	 multiple_gene_index<- grep(", ", MSRs_loocv2[,ns+1])
	 matrix_multiple_gene<- c()
		for(i in 1: length(multiple_gene_index))
		{
		 index_tmp<- multiple_gene_index[i]
		 gene_tmp<- unlist(strsplit(MSRs_loocv2[index_tmp,ns+1], ", "))
		 matrix_tmp<- matrix(NA, nrow=length(gene_tmp), ncol = ns+1)
			for(j in 1: length(gene_tmp))
			{
			 matrix_tmp[j, ]<- MSRs_loocv2[multiple_gene_index[i], ]
			 matrix_tmp[j, ns+1]<- gene_tmp[j]
			}
		 matrix_multiple_gene<- rbind(matrix_multiple_gene, matrix_tmp)
		}
	 new_summary_matrix_PDX<- rbind(MSRs_loocv2[-c(multiple_gene_index),], matrix_multiple_gene)
	 D<- split(as.data.frame(new_summary_matrix_PDX),as.factor(unlist(new_summary_matrix_PDX[,ns+1])))
	 rm_index_tmp<- which(sapply(1: length(D), function(k) nrow(D[[k]]))<=2)
	 D<- D[-c(rm_index_tmp)]
	 MSR_ref<- matrix(0, nrow = length(D), ncol = ns)
	 rownames(MSR_ref)<- names(D)
		for(i in 1: length(D))
		{
		 MSR_ref[i, ]<- apply(matrix(as.numeric(as.matrix(D[[i]][,1:ns])), ncol = ns), 2, median)
		}
	 
	 sMSR_ref<- t(t(MSR_ref) - apply(MSR_ref, 2 , median))

	 Cut_CNVs<- c(quantile(unlist(sMSR_ref), CNV_q/2, na.rm=T), quantile(unlist(sMSR_ref), 1-CNV_q/2, na.rm=T))
	 ##	Predict sMSR of new case	##
	 ns_new<- length(case)
	 nr_new<- nrow(case[[1]])
	 TMR_new<- data.frame(TMR = TMR_new)
	 predict_SR_cases<- matrix(0, nrow = nr_new - length(na_index), ncol = ns_new)
	 for(i in 1: ns_new)
	 {
	  cov_tmp<- sapply(1: nr_new, function(k) case[[i]][k,4])
	  cov_tmp<- cov_tmp[-na_index]
		for(j in 1: length(cov_tmp))
		{
		 PC_95CI_PDX_tmp<- predict(Reg_model_ref[[j]], TMR_new, interval = "prediction", level = .99, se.fit=T)
		 std_residuals_new_tmp<- (cov_tmp[j] - PC_95CI_PDX_tmp$fit[i,1])/(PC_95CI_PDX_tmp$residual.scale*sqrt(1 + (PC_95CI_PDX_tmp$se.fit[i]/PC_95CI_PDX_tmp$residual.scale)^2))
		 predict_SR_cases[j, i]<- std_residuals_new_tmp
		}
	 }
	 Genes_all<- c()
		for(i in 1: nrow(case[[1]]))
		{
		 chr_tmp<- sub("chr", "", ref[[1]][i,1])
		 start_tmp<- ref[[1]][i,2]
		 end_tmp<- ref[[1]][i,3]
		 gene_map_tmp<- gene_map[which(gene_map[,"Chromosome"]==chr_tmp), ]
		 region_tmp<- which(start_tmp<gene_map_tmp[,"Start"]& end_tmp>gene_map_tmp[,"End"]|start_tmp>gene_map_tmp[,"Start"]& end_tmp<gene_map_tmp[,"End"]|start_tmp<gene_map_tmp[,"Start"]& end_tmp>gene_map_tmp[,"Start"]|start_tmp<gene_map_tmp[,"End"]& end_tmp>gene_map_tmp[,"End"])
		 gene_tmp<- gene_map_tmp[region_tmp, ]
			if(length(nrow(gene_tmp))>0)
			{Genes_all[[i]]<- paste(gene_tmp[,"Gene_Symbol"], collapse=", ")}
			else
			{Genes_all[[i]]<- "NA"}	 
		}
	 Genes_all2<- Genes_all[-na_index]
	 predict_SR_cases2<- cbind(predict_SR_cases, Genes_all2)
	 multiple_gene_index<- grep(", ", predict_SR_cases2[,ns_new+1])
	 matrix_multiple_gene<- c()
		for(i in 1: length(multiple_gene_index))
		{
		 index_tmp<- multiple_gene_index[i]
		 gene_tmp<- unlist(strsplit(predict_SR_cases2[index_tmp,ns_new+1], ", "))
		 matrix_tmp<- matrix(NA, nrow=length(gene_tmp), ncol = ns_new+1)
			for(j in 1: length(gene_tmp))
			{
			 matrix_tmp[j, ]<- predict_SR_cases2[multiple_gene_index[i], ]
			 matrix_tmp[j, ns_new+1]<- gene_tmp[j]
			}
		 matrix_multiple_gene<- rbind(matrix_multiple_gene, matrix_tmp)
		}
	 new_summary_matrix_PDX<- rbind(predict_SR_cases2[-c(multiple_gene_index),], matrix_multiple_gene)
	 #new_summary_matrix_PDX<- new_summary_matrix_PDX[-which(new_summary_matrix_PDX[,ns_new+1]=="NA"), ]
	 D<- split(as.data.frame(new_summary_matrix_PDX),as.factor(unlist(new_summary_matrix_PDX[,ns_new+1])))
	 rm_index_tmp<- which(sapply(1: length(D), function(k) nrow(D[[k]]))<=2)
	 D<- D[-c(rm_index_tmp)]
	 MSR_case<- matrix(0, nrow = length(D), ncol = ns_new)
	 rownames(MSR_case)<- names(D)
		for(i in 1: length(D))
		{
		 MSR_case[i, ]<- apply(matrix(as.numeric(as.matrix(D[[i]][,1:ns_new])), ncol = ns_new), 2, median)
		}
	 sMSR_case<- t(t(MSR_case) - apply(MSR_case, 2 , median))
	 ##	Write an output	##
	 summary_table<- matrix(NA, nrow = nrow(sMSR_case), ncol = 1 + ns_new*2)
	 colnames(summary_table)<- c("Gene", paste("sMSR_S", 1:ns_new, sep=""), paste("CNV_S", 1:ns_new, sep=""))
	 summary_table[,"Gene"]<- rownames(sMSR_case)
	 summary_table[,paste("sMSR_S", 1:ns_new, sep="")]<- sMSR_case
		for(i in 1: ns_new)
		{
		 CNV_tmp<- rep("N", nrow(summary_table))
		 CNV_tmp[which(as.numeric(summary_table[,paste("sMSR_S", i, sep="")])< Cut_CNVs[1])]<- "D"
		 CNV_tmp[which(as.numeric(summary_table[,paste("sMSR_S", i, sep="")])> Cut_CNVs[2])]<- "A"
		 summary_table[,paste("CNV_S", i, sep="")]<- CNV_tmp
		}
	 write.table(summary_table, file="CNVs_summary.xls", col.name=T, row.name=F, sep="\t", quote=F)
	}

