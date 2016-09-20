
# Factors can be very annoying so I typically set stringsAsFactors to FALSE
options(stringsAsFactors = FALSE)


ef_import_model_sybil = function(x.version,x.name,x.type = "sybil",x.ext = ".rda",x.folder = "") {
  try({
    return(readRDS(paste0(x.folder,x.version,"_",x.name,"_","sybil.rda")))
  })
  cat("failed to load sybil model: ",paste0(x.folder,x.version,"_",x.name,"_","sybil.rda"))
  return(NULL)
}


ef_import_model_info = function(x.version, x.name, x.folder, x.file.ext = ".tbl.rda") {
  x.types = c("rxn_info","met_info","gene_info")
  x.files = paste0(x.folder, x.version,"_",x.name,"_",x.types,x.file.ext)
  
  if (any(grepl("\\.rda$",x.file.ext))) {
    x.info.list = setNames(x.files, x.types) %>% 
      lapply(function(x) as.tbl(readRDS(x)))
    return(x.info.list)
  } else {
    x.info.list = setNames(x.files, x.types) %>% 
      lapply(function(x) as.tbl(read.table(x, sep = "\t", quote = "", header = T, check.names = F, stringsAsFactors = F, fill = T)))
    return(x.info.list)
  }
}


flux_balance_analysis = function(model,x.digits = 8) {
  x.msg = paste0(model@react_id[model@obj_coef != 0],collapse = " ")
  if (nchar(x.msg) < 50) {
    cat("sybil maximization of",x.msg,": ")
  }
  x.result = optimizeProb(model, algorithm = "fba", retOptSol = T)
  
  if (x.result@lp_stat != 5) {
    cat(".....glpk error code", x.result@lp_stat,"\n")
    return(NA)
  }
  x.value = x.result@lp_obj %>% round(x.digits)
  cat(" ", x.value,'\n')
  return(x.value)
}


# ef_df is a simple function used to clean up character and numeric columns in a data.frame
ef_df = function(x.tbl,x.filter = c(),
                 x.names = names(x.tbl)[sapply(x.tbl,function(x) {
                   !any(class(x) %in% c("numeric","integer","logical"))})]) {
  x.names = x.names[x.names %in% names(x.tbl)]
  if (length(x.names) < 1) return(x.tbl %>% as.tbl)
  x.col = x.names[1]
  x.tbl[[x.col]] = ef_df_cln(x.tbl[[x.col]])
  if (x.col %in% x.filter) {
    x.tbl = x.tbl[nchar(x.tbl[[x.col]])>0,]
  }
  ef_df(x.tbl,x.filter,x.names[-1]) %>% as.tbl
}


# ef_df_clean is called from ef_df to trim excess white space from characters
ef_df_cln = function(x) {
  if (any(class(x) %in% c("numeric","integer","logical"))) return(x)
  x = ifelse(!is.na(x),as.character(x),"")
  if (!any(grepl("\\s",x))) return(x)
  gsub("^\\s+|\\s+$","",x)
}


ef_df_slice = function(x.tbl,x.column = "contrast_id",
                       x.id.list = unique(x.tbl[[x.column]]) %>% 
                         setNames(unique(x.tbl[[x.column]]))) {
  x.tbl = x.tbl %>% mutate_(slice_id = as.name(x.column))
  
  lapply(x.id.list,function(x.id) {
    x.tbl %>% filter(slice_id %in% x.id) %>% select(-slice_id)
  })
}


ef_read_table = function(x.filename, x.folder) {
  read.table(paste0(x.folder, x.filename), header=T, sep = "\t", quote = "", 
             check.names = F,  stringsAsFactors = F, fill = T) %>% as.tbl
}


ef_split_value = function(x.tbl,x.by = ";",x.melt = !all(c("variable","value") %in% names(x.tbl))) {
  x.split.count = sapply(x.tbl,function(x) sum(grepl(x.by,x))) %>% sort
  print(x.split.count)
  if (all(x.split.count == 0)) return(x.tbl)
  x.split.id = names(x.split.count)[x.split.count > 0]
  x.tbl.id = names(x.split.count)[x.split.count == 0]
  if (x.melt) {
    x.tbl.data = x.tbl %>% melt(x.tbl.id) %>% as.tbl
  } else {
    x.tbl.data = x.tbl %>% as.tbl
  }
  x.split = x.tbl.data$value %>% lapply(strsplit, x.by) %>% lapply(unlist)
  x.length = sapply(x.split,length)
  x.tbl.data %>% select(one_of(x.tbl.id %>% c("variable"))) %>%
    lapply(rep.int,x.length) %>% lapply(as.character) %>% 
    c(list(value = x.split %>% unlist)) %>% as_data_frame
}


ef_dataset_file_name = function(x.tbl,x.component = "",x.id = x.tbl$eset_id,
                                x.version = x.tbl$eset_version, x.type = x.tbl$eset_type,
                                x.folder = x.tbl$eset_folder,x.file.ext = x.tbl$eset_file_ext) {
  gsub("\\_\\.",".",paste0(x.folder,x.id,"_",x.version,"_",x.type,"_",x.component,x.file.ext))
}


ef_dataset_file_check = function(x.check) {
  x.check$eset_file_fdata = ef_dataset_file_name(x.tbl = x.check,x.component = "fdata")
  x.check$eset_file_pdata = ef_dataset_file_name(x.tbl = x.check,x.component = "pdata")
  x.check$eset_file_exprs = ef_dataset_file_name(x.tbl = x.check,x.component = "exprs")
  x.check$eset_root = ef_dataset_file_name(x.tbl = x.check, x.component = "", x.file.ext = "", x.folder = "")
  x.check$eset_file_count = file.exists(x.check$eset_file_fdata) + 
    file.exists(x.check$eset_file_pdata) + 
    file.exists(x.check$eset_file_exprs)
  
  x.check$eset_exists = x.check$eset_file_count > 0
  x.check$eset_ok = x.check$eset_file_count == 3
  x.check
}


ef_limma_file_check = function(x.check) {
  x.check$efit_file = ef_limma_file_name(x.tbl = x.check)
  x.check$efit_root = ef_limma_file_name(x.tbl = x.check,x.file.ext = "",x.folder = "")
  x.check$efit_exists = file.exists(x.check$efit_file)
  x.check$eset_file_efits = ef_dataset_file_name(x.tbl = x.check,x.component = "efits")
  x.check %>% group_by(eset_id,eset_type) %>% mutate(efit_ok = all(efit_exists)) %>% ungroup
}

ef_limma_file_name = function(x.tbl,x.lm = x.tbl$efit_lm, x.id = x.tbl$efit_id,
                                   x.version = x.tbl$efit_version, x.type = x.tbl$efit_type,
                                   x.folder = x.tbl$efit_folder, x.file.ext = x.tbl$efit_file_ext) {
  gsub("\\_\\.",".",paste0(x.folder,x.id,"_",x.version,"_",x.type,"_",x.lm,x.file.ext))
}


ef_dataset_preprocess = function(x.dataset,x.pdatabase,x.version = "gene", x.folder) {
  try({
    stopifnot(require(oligo))
    x.pdata = x.pdatabase %>% 
      dplyr::filter(dataset_id %in% x.dataset) %>% dplyr::mutate(pdata_id = sample_id) %>% 
      data.frame(row.names = "pdata_id",stringsAsFactors = F, check.names = F, check.rows = F) %>% 
      AnnotatedDataFrame
    
    x.celfiles = oligo::read.celfiles(filenames = x.pdata$sample_file_full, phenoData = x.pdata, rm.outliers = F)
    x.rma = oligo::rma(x.celfiles, normalize = T, background = T) 
    
    x.core = x.rma %>% ef_dataset_annotation(x.pdata$sample_annotation[1])
    x.eset = x.core %>% genefilter::featureFilter() %>% ef_dataset_annotation()
    x.nset.filtered = x.core %>% genefilter::nsFilter()
    x.nset = x.nset.filtered$eset %>% ef_dataset_annotation()
    fData(x.core)$feature_pass = fData(x.core)$probeset_id %in% fData(x.eset)$probeset_id
    
    x.core %>% ef_dataset_write(x.dataset, x.version = x.version, x.type = "core", x.folder = x.folder)
    x.eset %>% ef_dataset_write(x.dataset, x.version = x.version, x.type = "eset", x.folder = x.folder)
    x.nset %>% ef_dataset_write(x.dataset, x.version = x.version, x.type = "nset", x.folder = x.folder)
    cat(paste0(x.dataset," ESET OK ... "))
    gc()
    cat(" DONE\n")
    return(data_frame(dataset_id = x.dataset,dataset_status = "OK"))
  })
  
  gc()
  cat(paste0("\n\n",x.dataset," FAILED\n\n"))
  data_frame(dataset_id = x.dataset,dataset_status = "FAIL")
}


ef_dataset_annotation = function(x.eset, x.annotation = NULL){
  
  if (is.null(x.annotation)) {
    x.annotation = annotation(x.eset)
  }
  if (x.annotation == "skip") {
    return(x.eset)
  }
  x.feature.id = as.character(featureNames(x.eset))
  x.fdata = data.frame(
    feature_id = as.character(x.feature.id),
    probeset_id = as.character(x.feature.id),
    gene_id = as.character(unlist(annotate::lookUp(x.feature.id, x.annotation, c("ENTREZID"),load=T))),
    gene_symbol = as.character(unlist(annotate::lookUp(x.feature.id, x.annotation, c("SYMBOL"),load=T))),
    gene_name = as.character(unlist(annotate::lookUp(x.feature.id, x.annotation, c("GENENAME"),load=T))),
    check.rows = F,check.names = F,stringsAsFactors = F)
  x.fdata$gene_url = paste0("http://www.genecards.org/cgi-bin/carddisp.pl?gene=",x.fdata$gene_symbol)
  x.fdata$gene_unique = length(unique(x.fdata$feature_id)) == length(unique(x.fdata$gene_id))
  x.fdata$feature_annotation = as.character(x.annotation)
  x.fdata$feature_array = as.character(gsub("\\.db","",x.fdata$feature_annotation))
  x.fdata$feature_id = as.character(ifelse(x.fdata$gene_unique,x.fdata$gene_id,x.fdata$feature_id))
  x.fdata$feature_row_name = as.character(x.fdata$feature_id)
  rownames(x.fdata) = x.fdata$feature_row_name
  featureNames(x.eset) = as.character(x.fdata$feature_id)
  fData(x.eset) = x.fdata
  annotation(x.eset) = x.annotation
  x.eset
}


ef_dataset_write = function(x.eset,x.dataset, x.version, x.type, x.folder) {
  
  x.exprs.file = gzfile(paste0(x.folder, x.dataset,"_",x.version,"_",x.type,"_exprs.txt.gz"))
  x.pdata.file = gzfile(paste0(x.folder, x.dataset,"_",x.version,"_",x.type,"_pdata.txt.gz"))
  x.fdata.file = gzfile(paste0(x.folder, x.dataset,"_",x.version,"_",x.type,"_fdata.txt.gz"))
  
  x.exprs = exprs(x.eset) %>% 
    data.frame(check.names = F, check.rows = F, stringsAsFactors = F)
  x.pdata = pData(x.eset) %>% dplyr::mutate(pdata_id = rownames(pData(x.eset))) %>% 
    data.frame(check.names = F, check.rows = F, stringsAsFactors = F)
  x.fdata = fData(x.eset) %>% dplyr::mutate(fdata_id = rownames(fData(x.eset))) %>% 
    data.frame(check.names = F, check.rows = F, stringsAsFactors = F)
  write.table(x.exprs,file = x.exprs.file,sep = "\t",quote = F, row.names = T)
  write.table(x.pdata,file = x.pdata.file,sep = "\t",quote = F, row.names = F)
  write.table(x.fdata,file = x.fdata.file,sep = "\t",quote = F, row.names = F)
  data_frame(dataset_id = x.dataset, dataset_version = x.version, 
             dataset_type = x.type, dataset_folder = x.folder)
}


ef_dataset_load = function(x.tbl,x.eset.only = F) {
  
  x.dataset = x.tbl$eset_id[1]
  x.version = x.tbl$eset_version[1]
  x.type = x.tbl$eset_type[1]
  x.folder = x.tbl$eset_folder[1]
  x.file.ext = x.tbl$eset_file_ext[1]
  
  cat(paste0("Loading eset at: ",x.dataset,"_",x.version,"_",x.type,"..."))
  x.exprs.file = ef_dataset_file_name(x.tbl %>% dplyr::slice(1),"exprs")
  x.pdata.file = ef_dataset_file_name(x.tbl %>% dplyr::slice(1),"pdata")
  x.fdata.file = ef_dataset_file_name(x.tbl %>% dplyr::slice(1),"fdata")
  try({
    stopifnot(file.exists(x.pdata.file))
    stopifnot(file.exists(x.fdata.file))
    stopifnot(file.exists(x.exprs.file))
    
    x.pdata.load = read.table(file = x.pdata.file, sep = "\t", quote = "", header = T, check.names = F, stringsAsFactors = F)
    x.fdata.load = read.table(file = x.fdata.file, sep = "\t", quote = "", header = T, check.names = F, stringsAsFactors = F)
    x.exprs.load = read.table(file = x.exprs.file, sep = "\t", quote = "", header = T, check.names = F, stringsAsFactors = F)
    
    x.exprs = x.exprs.load %>% data.frame(check.names = F, check.rows = F, stringsAsFactors = F)
    rownames(x.exprs) = as.character(rownames(x.exprs))
    x.pdata = x.pdata.load %>% data.frame(row.names = "pdata_id",check.names = F, check.rows = F, stringsAsFactors = F)
    x.fdata = x.fdata.load %>% data.frame(row.names = "fdata_id",check.names = F, check.rows = F, stringsAsFactors = F)
    
    x.pdata$eset_id = x.dataset
    x.pdata$eset_type = x.type
    x.pdata$eset_version = x.version
    x.pdata$eset_folder = x.folder
    x.pdata$eset_file_ext = x.file.ext
    x.pdata$eset_root = paste0(x.dataset,"_",x.version,"_",x.type)
    
    stopifnot(identical(rownames(x.exprs),rownames(x.fdata)))
    stopifnot(identical(colnames(x.exprs),rownames(x.pdata)))
    
    x.eset = Biobase::ExpressionSet(assayData = as.matrix(x.exprs), 
                                    phenoData = AnnotatedDataFrame(x.pdata), 
                                    featureData = AnnotatedDataFrame(x.fdata))
    cat(" DONE\n")
    if (x.eset.only) {
      return(x.eset)
    }
    
    return(list(eset = x.eset, 
                dataset = x.dataset, version = x.version, type = x.type, folder = x.folder,
                info = x.tbl, status = "ok", 
                exprs = x.exprs %>% as.tbl %>% mutate(fdata_id = as.character(rownames(x.exprs))), 
                pdata = x.pdata %>% as.tbl %>% mutate(pdata_id = as.character(rownames(x.pdata))), 
                fdata = x.fdata %>% as.tbl %>% mutate(fdata_id = as.character(rownames(x.fdata)))))
  })
  if (x.eset.only) {
    return(NULL)
  }
  return(list(dataset = x.dataset, version = x.version, type = x.type, folder = x.folder,
              info = x.tbl, status = "fail"))
}


ef_limma_load = function(x.load) {
  x.load %>% ef_limma_file_check %>% filter(efit_exists) %>% 
    with(setNames(efit_file,efit_root)) %>% lapply(function(x.file) {
      Sys.sleep(0.1)
      cat(paste0("loading expression changes: ", x.file, '\n'))
      try({
        x.tbl = read.table(file = x.file, sep = "\t", quote = "", 
                           header = T, comment.char = "",
                           stringsAsFactors = F, check.names = F) %>% 
          as.tbl %>% mutate(efit_file = x.file) %>% ef_df_check_id("_id$") %>%
          left_join(x.load %>% 
                      select(starts_with("eset_"), starts_with("efit_")) %>% 
                      select(-one_of("efit_root")) %>% distinct)#, by = c("efit_id", "efit_file")
        return(x.tbl)
      })
      return(NULL)
    }) %>% bind_rows(.id = "efit_root")
}


ef_df_check_id = function(x.check,x.id.pattern = "_id$") {
  x.id.match = names(x.check)[grepl(x.id.pattern, names(x.check))]
  if (length(x.id.match) < 1) return(x.check)
  lapply(x.id.match,function(x.id) {
    x.check[[x.id]] <<- as.character(x.check[[x.id]])
  })
  return(x.check)
}


ef_limma_setup = function(x.contrasts, x.design, x.pdata, x.eset.folder, x.efit.folder,
                          x.type = "eset", x.version = "gene", x.ext = ".txt.gz",x.lm = "robust") {
  x.contrasts %>% lapply(as.list) %>% lapply(melt) %>% 
    lapply(as.tbl) %>% bind_rows(.id = "list_id") %>%
    transmute(dataset_id = list_id, contrast_abbrev = gsub("\\.","_",L1), contrast_string = value) %>% 
    mutate(contrast_id = paste0(dataset_id, "_", contrast_abbrev)) %>% 
    inner_join(x.design %>% lapply(as.character) %>% lapply("[",2) %>% melt %>% as.tbl %>% 
                 rename(dataset_id = L1, dataset_formula = value), by = "dataset_id") %>%
    inner_join(x.pdata, by = "dataset_id") %>% 
    mutate(organism_id = "hsa", eset_organism = organism_id,
           eset_id = dataset_id, eset_type = x.type, eset_version = x.version,
           eset_folder = x.eset.folder, eset_file_ext = x.ext) %>%
    ef_dataset_file_check %>%
    mutate(efit_id = contrast_id, efit_abbrev = contrast_abbrev, efit_dataset = eset_id,
           efit_contrast = contrast_string, efit_formula = dataset_formula,
           efit_type = "efit", efit_version = eset_version,
           efit_folder = x.efit.folder, efit_lm = x.lm, efit_file_ext = eset_file_ext) %>%
    ef_limma_file_check %>% 
    group_by(eset_id) %>% 
    mutate(eset_sample_count = length(unique(sample_id)),
           eset_contrast_count = length(unique(efit_id)),
           eset_contrast_list = paste0(unique(sort(efit_id)), collapse = ";")) %>% 
    ungroup %>% arrange(eset_sample_count, eset_contrast_count, eset_id, efit_id)
}


ef_limma_preprocess = function(x.setup) {
  require(limma)
  try({
    x.ebayes = ef_limma_efit(x.setup = x.setup)
    x.status = x.setup %>% select(eset_id, eset_version, eset_type, eset_folder, eset_file_ext,
                                  efit_id, efit_version, efit_type, efit_folder, efit_file_ext, efit_lm) %>%
      ef_limma_file_check() %>% ef_df_slice("efit_root") %>% 
      lapply(function(efit.setup,efit.data,efit.coef = efit.setup$efit_id,efit.digits = 4) {
        stopifnot(nrow(efit.setup) == 1 && length(unique(efit.coef)) == 1)
        efit.top.tbl = topTable(efit.data, coef = efit.coef[1], number = Inf, sort.by = "none") %>% as.tbl %>% 
          mutate(efit_id = efit.coef[1],
                 ftest_pval = signif(ftest_pval, efit.digits),
                 ftest_fdr = signif(ftest_fdr, efit.digits),
                 ftest_ave = signif(ftest_ave, efit.digits),
                 logFC = signif(logFC, efit.digits),
                 AveExpr = signif(AveExpr, efit.digits),
                 P.Value = signif(P.Value, efit.digits),
                 adj.P.Val = signif(adj.P.Val, efit.digits)) %>% 
          select(efit_id, feature_id, 
                 probeset_id, gene_id, gene_symbol, gene_name, 
                 logfc = logFC, ave = AveExpr, pval = P.Value, fdr = adj.P.Val,
                 ftest_pval = ftest_pval, ftest_fdr = ftest_fdr,
                 ftest_ave = ftest_ave)
        cat(paste0("# genes: ", nrow(efit.top.tbl),"\n"))
        try({
          Sys.sleep(.5)
          write.table(efit.top.tbl,file = gzfile(efit.setup$efit_file[1]), sep="\t", quote = F, row.names = F, append = FALSE)
          Sys.sleep(.5)
          return(efit.setup %>% mutate(efit_status = "OK"))
        })
        return(efit.setup %>% mutate(efit_status = "fail"))
      },efit.data = x.ebayes) %>% bind_rows
    return(x.status %>% mutate(eset_status = ifelse(all(efit_status == "OK"), "OK", "fail")))
  })
  return(x.setup %>% mutate(eset_status = "fail", efit_status = "fail"))

}


ef_limma_efit = function(x.setup) {
  
  stopifnot(length(unique(x.setup$eset_id)) == 1)
  stopifnot(length(unique(x.setup$eset_version)) == 1)
  stopifnot(length(unique(x.setup$eset_type)) == 1)
  stopifnot(length(unique(x.setup$eset_folder)) == 1)
  stopifnot(length(unique(x.setup$efit_formula)) == 1)
  stopifnot(length(unique(x.setup$efit_lm)) == 1)
  
  x.contrasts = x.setup$efit_contrast %>% setNames(x.setup$efit_id)
  x.formula = as.formula(paste0("~",x.setup$efit_formula[1]))
  x.eset = ef_dataset_load(x.setup %>% dplyr::slice(1),T)
  x.pdata = pData(x.eset)
  rownames(x.pdata) = x.pdata$sample_id
  x.design = model.matrix(x.formula, data = x.pdata)
  # print(x.design)
  x.attr = names(attr(x.design,"contrasts")) %>% paste0(collapse = "|")
  colnames(x.design) = gsub("\\(|\\)","",gsub("\\.+","_",gsub(x.attr,"",colnames(x.design))))
  rownames(x.design) = x.pdata$sample_id
  
  x.make.contrasts = data.frame(lapply(x.contrasts,function(x) limma::makeContrasts(contrasts = x,levels=x.design)))
  names(x.make.contrasts) = names(x.contrasts)
  
  x.robust = all(x.setup$efit_lm[1] == "robust")
  x.ebayes = limma::eBayes(contrasts.fit(lmFit(
    x.eset, x.design, method = x.setup$efit_lm[1], 
    maxit=1000), x.make.contrasts), robust = x.robust)
  x.ebayes$genes = fData(x.eset)
  x.ebayes$genes$feature_id = as.character(x.ebayes$genes$feature_id)
  x.ftest = limma::topTableF(x.ebayes,number = Inf,sort.by = "none")
  x.ebayes$genes$ftest_pval = x.ftest$P.Value
  x.ebayes$genes$ftest_fdr = x.ftest$adj.P.Val
  x.ebayes$genes$ftest_ave = x.ftest$AveExpr
  x.ebayes
}


ef_timbr_weights = function(gene.expression, gpr.setup, gene.ignore = c(), gene.keep.all) {
  
  stopifnot(length(unique(gene.expression$limma_id)) == 1)
  gpr.filter = gpr.setup %>% ef_df("value") %>% 
    filter(variable %in% gene.expression[["organism_id"]]) 
  gpr.isozyme.value = ef_timbr_summarize_isozyme(gpr.filter, gene.expression, gene.ignore, gene.keep.all, 3, -3)
  gpr.complex.value = ef_timbr_summarize_complex(gpr.filter, gene.expression, gene.ignore, gene.keep.all, 3, -3)
  gpr.subunit.value = ef_timbr_summarize_subunit(gpr.filter, gene.expression, gene.ignore, gene.keep.all, 3, -3)
  gpr.setup %>% select(rxn_id, organism_id = variable,timbr_weight_default) %>% distinct %>%
    inner_join(gene.expression %>% select(limma_id, organism_id, drug_id, time_id, dose_id) %>% distinct) %>% 
    left_join(bind_rows(list(gpr.isozyme.value,gpr.subunit.value,gpr.complex.value)) %>%
                select(organism_id, rxn_id, gpr_significance, gpr_logfc)) %>% #organism_id, drug_id, time_id, dose_id, 
    mutate(gpr_significance = ifelse(!is.na(gpr_significance),gpr_significance,0),
           gpr_logfc = ifelse(!is.na(gpr_logfc),gpr_logfc,0),
           relative_weight_ctl = (2 ^ gpr_logfc),
           relative_weight_trt = (2 ^ -gpr_logfc),
           timbr_weight_ctl = signif(timbr_weight_default * relative_weight_ctl,3),
           timbr_weight_def = timbr_weight_default,
           timbr_weight_trt = signif(timbr_weight_default * relative_weight_trt,3))
}


ef_timbr_gpr_join = function(x.gpr, x.gene,x.keep = F) {
  if (x.keep) {
    x.tbl = x.gpr %>% mutate(gene_id = as.character(gene_id)) %>% 
      left_join(x.gene %>% select(gene_id, organism_id, logfc, fdr, pval, ave) %>%
                  mutate(gene_id = as.character(gene_id)), by = c("organism_id","gene_id")) %>%
      mutate(fdr = ifelse(!is.na(fdr),fdr,1),
             logfc = ifelse(!is.na(logfc),logfc,0))
    x.rule.count = length(unique(x.gpr$gene_id))
    x.data.count = length(unique(x.gene$gene_id))
    x.join.count = length(unique(x.tbl$gene_id))
    cat(paste0(x.join.count, " / ", x.rule.count , " genes (ALL) used to calculate weights (", x.data.count, " with expression)\n"))
  } else {
    x.tbl = x.gpr %>% mutate(gene_id = as.character(gene_id)) %>% 
      inner_join(x.gene %>% select(gene_id, organism_id, logfc, fdr, pval, ave) %>%
                   mutate(gene_id = as.character(gene_id)), by = c("organism_id","gene_id")) %>%
      mutate(fdr = ifelse(!is.na(fdr),fdr,1),
             logfc = ifelse(!is.na(logfc),logfc,0))
    x.rule.count = length(unique(x.gpr$gene_id))
    x.data.count = length(unique(x.gene$gene_id))
    x.join.count = length(unique(x.tbl$gene_id))
    cat(paste0(x.join.count, " / ", x.rule.count , " genes used to calculate weights (", x.data.count, " with expression)\n"))
  }
  if (!is.null(x.tbl$logfc) && !is.null(x.tbl$fdr)) {
    return(x.tbl %>% mutate(mag = sign(logfc) * -log10(fdr)))
  }
  x.tbl
}


ef_timbr_summarize_subunit = function(gpr.info, gene.value, gene.ignore, gene.keep.all, max.value, min.value) {
  gpr.subunit.split = gpr.info %>% filter(grepl("\\(",value)) %>%
    mutate(gpr_rule = gsub("\\;","|",gsub("\\:","&",value))) %>%
    ef_split_value("\\:") %>% ef_df
  gpr.subunit.unique = gpr.subunit.split %>% select(variable,value) %>% distinct %>% 
    mutate(subunit_id = paste0("subunit",1:n()))
    #       subunit_rule = gsub("[ \\(\\)]","",gsub("\\;"," | ",value)))
  gpr.subunit.breakdown = gpr.subunit.split %>% 
    left_join(gpr.subunit.unique, by = c("variable","value")) %>% 
    mutate(value = gsub("\\)$","",gsub("^\\(","",value))) %>% 
    ef_split_value("\\;") %>% ef_df("value") %>% filter(!grepl("^0$",value)) %>% 
    rename(organism_id = variable, gene_id = value) %>%
    mutate(gene_ignore = gene_id %in% gene.ignore) #%>%
    #group_by(rxn_id, organism_id) %>% 
    #mutate(n_gene = length(unique(gene_id)),
    #       n_ignore = length(unique(gene_id[gene_ignore])),
    #       n_subunit = length(unique(subunit_rule))) %>% ungroup
  gpr.subunit.breakdown %>% ef_timbr_gpr_join(gene.value, gene.keep.all) %>% 
    mutate(gene_significance = ifelse(!gene_ignore, ifelse(!is.na(mag),mag, 0), NA)) %>%
    mutate(gene_logfc = ifelse(!gene_ignore, ifelse(!is.na(logfc),logfc, 0), NA)) %>%
    group_by(organism_id, rxn_id, subunit_id) %>% #organism_id, drug_id, time_id, dose_id, n_gene, n_ignore, n_subunit, 
    summarize(gpr_significance = mean(gene_significance, na.rm = T),
              gpr_logfc = mean(gene_logfc, na.rm = T)) %>% ungroup %>%
    arrange(desc(abs(gpr_significance)),desc(abs(gpr_logfc))) %>% 
    group_by(organism_id, rxn_id) %>% #organism_id, drug_id, time_id, dose_id, n_gene, n_ignore, n_subunit, 
    slice(1) %>% ungroup %>%
    mutate(gpr_logfc = pmax(pmin(ifelse(!is.na(gpr_logfc),gpr_logfc, 0), max.value), min.value)) %>%
    mutate(gpr_significance = pmax(pmin(ifelse(!is.na(gpr_significance),gpr_significance, 0), max.value), min.value))
}


ef_timbr_summarize_complex = function(gpr.info, gene.value, gene.ignore, gene.keep.all, max.value, min.value) {
  
  gpr.complex.split = gpr.info %>% 
    filter(grepl("\\:",value), !grepl("\\(",value)) %>% 
    mutate(gpr_rule = gsub("\\;","|",gsub("\\:","&",value))) %>%
    ef_split_value(";") %>% ef_df
  gpr.complex.unique = gpr.complex.split %>% select(variable,value) %>% distinct %>% 
    mutate(complex_id = paste0("complex",1:n()), complex_rule = gsub("\\:"," & ",value))
  gpr.complex.breakdown = gpr.complex.split %>% 
    left_join(gpr.complex.unique, by = c("variable","value")) %>% 
    ef_split_value("\\:") %>% ef_df("value") %>% filter(!grepl("^0$",value)) %>% 
    rename(organism_id = variable, gene_id = value) %>% 
    mutate(gene_ignore = gene_id %in% gene.ignore) %>%
    group_by(rxn_id, organism_id) %>% 
    mutate(n_gene = length(unique(gene_id)),
           n_ignore = length(unique(gene_id[gene_ignore])),
           n_complex = length(unique(complex_rule))) %>% ungroup %>%
    mutate(n_subunit = n_gene)
  gpr.complex.breakdown %>% ef_timbr_gpr_join(gene.value, gene.keep.all) %>% 
    mutate(gene_significance = ifelse(!gene_ignore, ifelse(!is.na(mag),mag, 0), NA)) %>%
    mutate(gene_logfc = ifelse(!gene_ignore, ifelse(!is.na(logfc),logfc, 0), NA)) %>%
    arrange(desc(abs(gene_significance)), desc(abs(gene_logfc))) %>% 
    group_by(organism_id, rxn_id, complex_id) %>%
    slice(1) %>% ungroup %>% 
    group_by(organism_id, rxn_id) %>%
    summarize(gpr_significance = mean(gene_significance, na.rm = T),
              gpr_logfc = mean(gene_logfc, na.rm = T)) %>% ungroup %>%
    mutate(gpr_logfc = pmax(pmin(ifelse(!is.na(gpr_logfc), gpr_logfc, 0), max.value), min.value)) %>%
    mutate(gpr_significance = pmax(pmin(ifelse(!is.na(gpr_significance), gpr_significance, 0), max.value), min.value))
}


ef_timbr_summarize_isozyme = function(gpr.info, gene.value, gene.ignore, gene.keep.all, max.value, min.value) {
  gpr.isozyme.split = gpr.info %>% 
    filter(!grepl("\\:",value), !grepl("\\(",value)) %>% 
    mutate(gpr_rule = gsub("\\;","|",gsub("\\:","&",value))) %>%
    ef_split_value(";") %>% ef_df
  gpr.isozyme.breakdown = gpr.isozyme.split %>% ef_df("value") %>% filter(!grepl("^0$",value)) %>% 
    rename(organism_id = variable, gene_id = value) %>%
    mutate(gene_ignore = gene_id %in% gene.ignore) %>%
    group_by(rxn_id,organism_id) %>% 
    mutate(n_gene = length(unique(gene_id)),
           n_ignore = length(unique(gene_id[gene_ignore])),
           n_subunit = 1) %>% ungroup
  gpr.isozyme.breakdown %>% ef_timbr_gpr_join(gene.value, x.keep = gene.keep.all) %>%
    mutate(gene_significance = ifelse(!gene_ignore, ifelse(!is.na(mag),mag, 0), NA)) %>%
    mutate(gene_logfc = ifelse(!gene_ignore, ifelse(!is.na(logfc),logfc, 0), NA)) %>%
    group_by(organism_id, rxn_id) %>%
    summarize(gpr_significance = mean(gene_significance, na.rm = T),
              gpr_logfc = mean(gene_logfc, na.rm = T)) %>% ungroup %>%
    mutate(gpr_logfc = pmax(pmin(ifelse(!is.na(gpr_logfc), gpr_logfc, 0), max.value), min.value)) %>%
    mutate(gpr_significance = pmax(pmin(ifelse(!is.na(gpr_significance), gpr_significance, 0), max.value), min.value))
}


ef_gpr_classification = function(gpr.rules) {
  gpr.rules %>% 
    mutate(rule_empty = is.na(value) | grepl("^0$", value) | grepl("^$", value)) %>%
    mutate(rule_notempty = !rule_empty & grepl("[0-9]+", value)) %>%
    mutate(rule_parentheses = grepl("\\(|\\)", value)) %>%
    mutate(rule_multiprotein = grepl("\\:", value)) %>% 
    mutate(rule_redundant = grepl("\\;", value)) %>%
    mutate(rule_type = ifelse(rule_notempty, ifelse(
      rule_parentheses, 
      ifelse(rule_multiprotein & rule_redundant, "subunit", "invalid1"),
      ifelse(rule_multiprotein, ifelse(rule_redundant, "invalid2", "complex"),
             ifelse(rule_redundant, "isozyme", "unizyme"))), "empty"))
}


ef_efmin_gpr_join = function(x.gpr, x.gene,x.keep, x.default = 0) {
  if (x.keep[1]) {
    x.tbl = x.gpr %>% mutate(gene_id = as.character(gene_id)) %>% 
      left_join(x.gene %>% select(organism_id, gene_id, logfc, fdr, pval, ave, ctl, trt) %>% 
                  mutate(gene_id = as.character(gene_id)), by = c("organism_id","gene_id")) %>%
      mutate(fdr = ifelse(!is.na(fdr),fdr,1),
             logfc = ifelse(!is.na(logfc), logfc, 0),
             ave = ifelse(!is.na(ave), ave, x.default),
             ctl = ifelse(!is.na(ctl), ctl, x.default),
             trt = ifelse(!is.na(trt), trt, x.default))
    x.rule.count = length(unique(x.gpr$gene_id))
    x.data.count = length(unique(x.gene$gene_id))
    x.join.count = length(unique(x.tbl$gene_id))
    cat(paste0(x.join.count, " / ", x.rule.count , " genes (ALL) used to calculate weights (", x.data.count, " with expression)\n"))
  } else {
    
    x.tbl = x.gpr %>% mutate(gene_id = as.character(gene_id)) %>% 
      inner_join(x.gene %>% select(organism_id, gene_id, logfc, fdr, pval, ave, ctl, trt) %>% 
                   mutate(gene_id = as.character(gene_id)), by = c("organism_id","gene_id"))
    x.rule.count = length(unique(x.gpr$gene_id))
    x.data.count = length(unique(x.gene$gene_id))
    x.join.count = length(unique(x.tbl$gene_id))
    cat(paste0(x.join.count, " / ", x.rule.count , " genes used to calculate weights (", x.data.count, " with expression)\n"))
  }
  x.tbl
}


ef_efmin_weights = function(gene.expression, gpr.setup, gene.ignore = c(), gene.keep.all = F) {
  stopifnot(length(unique(gene.expression$limma_id)) == 1)
  gpr.type = gpr.setup %>% ef_gpr_classification %>%
    filter(rule_notempty) %>%
    filter(variable %in% gene.expression[["organism_id"]]) 
  print(gpr.type %>% count(rule_type))
  gpr.isozyme.value = gpr.type %>% filter(rule_type %in% c("isozyme", "unizyme")) %>% ef_efmin_summarize_isozyme(
    gene.expression, gene.ignore = gene.ignore, gene.keep.all = gene.keep.all, max.value = 10, min.value = 0)
  gpr.complex.value = gpr.type %>% filter(rule_type %in% c("complex")) %>% ef_efmin_summarize_complex(
    gene.expression, gene.ignore = gene.ignore, gene.keep.all = gene.keep.all, max.value = 10, min.value = 0)
  gpr.subunit.value = gpr.type %>% filter(rule_type %in% c("subunit")) %>% ef_efmin_summarize_subunit(
    gene.expression, gene.ignore = gene.ignore, gene.keep.all = gene.keep.all, max.value = 10, min.value = 0)
  # gpr.filter = gpr.setup %>% ef_df("value") %>% 
  #  filter(variable %in% gene.expression[["organism_id"]]) 
  # gpr.isozyme.value = ef_efmin_summarize_isozyme(gpr.filter, gene.expression,gene.ignore = gene.ignore, gene.keep.all = gene.keep.all,3)
  # gpr.complex.value = ef_efmin_summarize_complex(gpr.filter, gene.expression,gene.ignore = gene.ignore, gene.keep.all = gene.keep.all,3)
  # gpr.subunit.value = ef_efmin_summarize_subunit(gpr.filter, gene.expression,gene.ignore = gene.ignore, gene.keep.all = gene.keep.all,3)
  gpr.setup %>% 
    ef_gpr_classification %>%
    select(rxn_id, rule_type, organism_id = variable, efmin_weight_default = timbr_weight_default) %>% distinct %>%
    inner_join(gene.expression %>% select(limma_id, organism_id, drug_id, time_id, dose_id) %>% distinct) %>% 
    left_join(bind_rows(list(
      gpr.isozyme.value %>% mutate(gpr_type = "isozyme"),
      gpr.subunit.value %>% mutate(gpr_type = "subunit"),
      gpr.complex.value %>% mutate(gpr_type = "complex"))) %>%
        select(organism_id, rxn_id, gpr_type, gpr_ave, gpr_ctl, gpr_trt)) %>%
    mutate(gpr_ave = ifelse(!is.na(gpr_ave), gpr_ave, 0),
           gpr_ctl = ifelse(!is.na(gpr_ctl), gpr_ctl, 0),
           gpr_trt = ifelse(!is.na(gpr_trt), gpr_trt, 0),
           efmin_weight_ctl = signif(2 ^ -gpr_ctl, 3),
           efmin_weight_ave = signif(2 ^ -gpr_ave, 3),
           efmin_weight_trt = signif(2 ^ -gpr_trt, 3))
}


ef_efmin_summarize_isozyme = function(gpr.info, gene.value, gene.ignore, gene.keep.all, max.value = 10, min.value = 0) {
  gpr.isozyme.split = gpr.info %>% 
    filter(!grepl("\\:",value), !grepl("\\(",value)) %>% 
    #mutate(gpr_rule = gsub("\\;","|",gsub("\\:","&",value))) %>%
    select(rxn_id, variable, value) %>% 
    ef_split_value(";") %>% ef_df
  gpr.isozyme.breakdown = gpr.isozyme.split %>% ef_df("value") %>% filter(!(value %in% c("0"))) %>%
    rename(organism_id = variable, gene_id = value) %>%
    mutate(gene_ignore = gene_id %in% gene.ignore)
  gpr.isozyme.breakdown %>% ef_efmin_gpr_join(x.gene = gene.value, x.keep = gene.keep.all, x.default = min.value) %>% 
    mutate(gene_ave = ifelse(!gene_ignore, ifelse(!is.na(ave), ave, min.value), NA)) %>% 
    mutate(gene_ctl = ifelse(!gene_ignore, ifelse(!is.na(ctl), ctl, min.value), NA)) %>% 
    mutate(gene_trt = ifelse(!gene_ignore, ifelse(!is.na(trt), trt, min.value), NA)) %>% 
    group_by(organism_id, rxn_id) %>% 
    summarize(gpr_significance = 1,
              gpr_ave = max(gene_ave, na.rm = T),
              gpr_ctl = max(gene_ctl, na.rm = T),
              gpr_trt = max(gene_trt, na.rm = T)) %>% ungroup %>%
    mutate(gpr_ave = pmin(pmax(ifelse(!is.na(gpr_ave), gpr_ave, min.value), min.value), max.value)) %>% 
    mutate(gpr_ctl = pmin(pmax(ifelse(!is.na(gpr_ctl), gpr_ctl, min.value), min.value), max.value)) %>% 
    mutate(gpr_trt = pmin(pmax(ifelse(!is.na(gpr_trt), gpr_trt, min.value), min.value), max.value))
}


ef_efmin_summarize_subunit = function(gpr.info, gene.value, gene.ignore, gene.keep.all, max.value = 10, min.value = 0) {
  gpr.subunit.split = gpr.info %>% filter(grepl("\\(",value)) %>%
    # mutate(gpr_rule = gsub("\\;","|",gsub("\\:","&",value))) %>%
    select(rxn_id, variable, value) %>% 
    ef_split_value("\\:") %>% ef_df
  gpr.subunit.unique = gpr.subunit.split %>% select(variable,value) %>% distinct %>% 
    mutate(subunit_id = paste0("subunit",1:n()))
  gpr.subunit.breakdown = gpr.subunit.split %>% 
    left_join(gpr.subunit.unique, by = c("variable","value")) %>% 
    mutate(value = gsub("\\)$","",gsub("^\\(","",value))) %>% 
    ef_split_value("\\;") %>% ef_df("value") %>% filter(!(value %in% c("0"))) %>%
    rename(organism_id = variable, gene_id = value) %>%
    mutate(gene_ignore = gene_id %in% gene.ignore)
  gpr.subunit.breakdown %>% ef_efmin_gpr_join(x.gene = gene.value, x.keep = gene.keep.all, x.default = min.value) %>% 
    mutate(gene_ave = ifelse(!gene_ignore, ifelse(!is.na(ave), ave, min.value), NA)) %>% 
    mutate(gene_ctl = ifelse(!gene_ignore, ifelse(!is.na(ctl), ctl, min.value), NA)) %>% 
    mutate(gene_trt = ifelse(!gene_ignore, ifelse(!is.na(trt), trt, min.value), NA)) %>% 
    group_by(organism_id, rxn_id, subunit_id) %>% 
    summarize(subunit_significance = 1,
              subunit_ave = max(gene_ave, na.rm = T),
              subunit_ctl = max(gene_ctl, na.rm = T),
              subunit_trt = max(gene_trt, na.rm = T)) %>% ungroup %>%
    mutate(subunit_ave = pmin(pmax(ifelse(!is.na(subunit_ave), subunit_ave, min.value), min.value), max.value)) %>% 
    mutate(subunit_ctl = pmin(pmax(ifelse(!is.na(subunit_ctl), subunit_ctl, min.value), min.value), max.value)) %>% 
    mutate(subunit_trt = pmin(pmax(ifelse(!is.na(subunit_trt), subunit_trt, min.value), min.value), max.value)) %>% 
    group_by(organism_id, rxn_id) %>% 
    summarize(gpr_significance = 1,
              gpr_ave = min(subunit_ave, na.rm = T),
              gpr_ctl = min(subunit_ctl, na.rm = T),
              gpr_trt = min(subunit_trt, na.rm = T)) %>% ungroup %>%
    mutate(gpr_ave = pmin(pmax(ifelse(!is.na(gpr_ave), gpr_ave, min.value), min.value), max.value)) %>% 
    mutate(gpr_ctl = pmin(pmax(ifelse(!is.na(gpr_ctl), gpr_ctl, min.value), min.value), max.value)) %>% 
    mutate(gpr_trt = pmin(pmax(ifelse(!is.na(gpr_trt), gpr_trt, min.value), min.value), max.value))
}


ef_efmin_summarize_complex = function(gpr.info, gene.value, gene.ignore, gene.keep.all, max.value = 10, min.value = 0) {
  
  gpr.complex.split = gpr.info %>% 
    filter(grepl("\\:",value), !grepl("\\(",value)) %>% 
    #mutate(gpr_rule = gsub("\\;","|",gsub("\\:","&",value))) %>%
    ef_split_value(";") %>% ef_df
  gpr.complex.unique = gpr.complex.split %>% select(variable,value) %>% distinct %>% 
    mutate(complex_id = paste0("complex",1:n()))#, complex_rule = gsub("\\:"," & ",value))
  gpr.complex.breakdown = gpr.complex.split %>% 
    left_join(gpr.complex.unique, by = c("variable","value")) %>% 
    ef_split_value("\\:") %>% ef_df("value") %>% filter(!(value %in% c("0"))) %>%
    rename(organism_id = variable, gene_id = value) %>% 
    mutate(gene_ignore = gene_id %in% gene.ignore)
  gpr.complex.breakdown %>% ef_efmin_gpr_join(x.gene = gene.value, x.keep = gene.keep.all, x.default = min.value) %>% 
    mutate(gene_ave = ifelse(!gene_ignore, ifelse(!is.na(ave), ave, min.value), NA)) %>% 
    mutate(gene_ctl = ifelse(!gene_ignore, ifelse(!is.na(ctl), ctl, min.value), NA)) %>% 
    mutate(gene_trt = ifelse(!gene_ignore, ifelse(!is.na(trt), trt, min.value), NA)) %>% 
    group_by(organism_id, rxn_id, complex_id) %>% 
    summarize(complex_significance = 1,
              complex_ave = min(gene_ave, na.rm = T),
              complex_ctl = min(gene_ctl, na.rm = T),
              complex_trt = min(gene_trt, na.rm = T)) %>% ungroup %>%
    mutate(complex_ave = pmin(pmax(ifelse(!is.na(complex_ave), complex_ave, min.value), min.value), max.value)) %>% 
    mutate(complex_ctl = pmin(pmax(ifelse(!is.na(complex_ctl), complex_ctl, min.value), min.value), max.value)) %>% 
    mutate(complex_trt = pmin(pmax(ifelse(!is.na(complex_trt), complex_trt, min.value), min.value), max.value)) %>% 
    group_by(organism_id, rxn_id) %>% 
    summarize(gpr_significance = 1,
              gpr_ave = max(complex_ave, na.rm = T),
              gpr_ctl = max(complex_ctl, na.rm = T),
              gpr_trt = max(complex_trt, na.rm = T)) %>% ungroup %>%
    mutate(gpr_ave = pmin(pmax(ifelse(!is.na(gpr_ave), gpr_ave, min.value), min.value), max.value)) %>% 
    mutate(gpr_ctl = pmin(pmax(ifelse(!is.na(gpr_ctl), gpr_ctl, min.value), min.value), max.value)) %>% 
    mutate(gpr_trt = pmin(pmax(ifelse(!is.na(gpr_trt), gpr_trt, min.value), min.value), max.value))
}


ef_reproduce_figure2 = function(gpr.info,gpr.limit = 9,gpr.legend = F) {
  gpr.color = c(`non-enzymatic` = "#7F7F7F", shared = "#674EA7", 
                `rat-specific` = "#C0504D", `human-specific` = "#4F81BD")
  gpr.size = gpr.info %>% 
    mutate(rno = pmin(n_gene_rno,gpr.limit),hsa = pmin(n_gene_hsa,gpr.limit)) %>%
    mutate(organism = ifelse(rno > 0 & hsa > 0,"shared",ifelse(
      rno > 0,"rat-specific",ifelse(hsa > 0,"human-specific","non-enzymatic"))))
  gpr.size.count = gpr.size %>% count(rno,hsa) %>% ungroup %>% mutate(n_gpr = n) %>%
    mutate(organism = ifelse(rno > 0 & hsa > 0,"shared",ifelse(
      rno > 0,"rat-specific",ifelse(hsa > 0,"human-specific","non-enzymatic")))) %>% 
    mutate(color = gpr.color[organism])
  gpr.size.count %>% 
    ggplot(aes(x = factor(hsa), y = factor(rno)))+
    geom_tile(aes(fill = organism, alpha = log10(n_gpr)), size = 0.9) + #, height = 0.9
    scale_fill_manual(values = gpr.color) + scale_color_manual(values = gpr.color) +
    scale_alpha_continuous(range = c(.2,.9)) + 
    scale_x_discrete(breaks = 0:gpr.limit, labels = gsub(as.character(gpr.limit), paste0(gpr.limit,"+"), 0:gpr.limit)) + 
    scale_y_discrete(breaks = 0:gpr.limit, labels = gsub(as.character(gpr.limit), paste0(gpr.limit,"+"), 0:gpr.limit)) + 
    xlab("Human GPR size") + ylab("Rat GPR size") + 
    theme_bw(base_size = 12) +
    theme(title = element_text(size = 12), legend.position = ifelse(gpr.legend, "top", "none"),
          axis.text.y = element_text(hjust = .7),
          panel.grid.major = element_line(size = NA),
          panel.grid.minor = element_line(size = NA))
}


ef_reproduce_figure3 = function(gpr.info,gpr.limit = 9) {
  gpr.color = c(`non-enzymatic` = "#7F7F7F", shared = "#674EA7", 
                `rat-specific` = "#C0504D", `human-specific` = "#4F81BD")
  gpr.background = matrix(1,nrow = gpr.limit+1,ncol = gpr.limit+1,dimnames = list(
    0:gpr.limit,0:gpr.limit)) %>% melt %>% as.tbl %>% select(hsa = Var1,rno = Var2) %>% 
    mutate(organism = ifelse(rno > 0 & hsa > 0,"shared",ifelse(
      rno > 0,"rat-specific",ifelse(hsa > 0,"human-specific","non-enzymatic"))))
  gpr.size = gpr.info %>% mutate(rno = pmin(n_rno,gpr.limit),hsa = pmin(n_hsa,gpr.limit)) %>%
    mutate(organism = ifelse(rno > 0 & hsa > 0,"shared",ifelse(
      rno > 0,"rat-specific",ifelse(hsa > 0,"human-specific","non-enzymatic"))))
  gpr.size %>% ggplot(aes(x = factor(hsa), y = factor(rno)))+
    geom_point(aes(color = organism),alpha = 0.3, size = 1,
               position = position_jitter(width = 0.38, height = 0.38)) + #
    geom_tile(data = gpr.background,aes(fill = organism),
              color = NA,alpha = 0.1,width = .95,height = .95) +
    scale_fill_manual(values = gpr.color)+ scale_color_manual(values = gpr.color) +
    scale_x_discrete(breaks = 0:gpr.limit,labels = gsub(as.character(gpr.limit),paste0(gpr.limit,"+"),0:gpr.limit)) + 
    scale_y_discrete(breaks = 0:gpr.limit,labels = gsub(as.character(gpr.limit),paste0(gpr.limit,"+"),0:gpr.limit)) + 
    xlab("Human GPR size") + ylab("Rat GPR size") + 
    theme_bw(base_size = 12) +
    theme(legend.position="none",
          axis.text.y = element_text(hjust = .7),
          panel.grid.major = element_line(size = NA),
          panel.grid.minor = element_line(size = NA))
}


ef_biomarker_gene_deletion = function(x.tbl,x.sybil) {
  x.genes = x.tbl %>% select(gene_id) %>% 
    mutate(gene_id = as.character(gene_id)) %>% distinct %>% 
    filter(gene_id %in% x.sybil@allGenes) %>% arrange(gene_id) %>% with(gene_id)
  print(x.genes)
  print(intersect(x.genes,x.sybil@allGenes))
  
  x.rxns = c()
  try({
    x.rxns = x.sybil@react_id[x.sybil %>% geneDel(genes = x.genes,checkId = T)]
  })
  
  x.tbl %>% mutate(rxn_direct = rxn_id %in% x.rxns) %>% 
    select(starts_with("iem"), starts_with("rxn"), starts_with("biomarker")) %>% 
    distinct %>% ungroup %>% arrange(iem_id, rxn_id, biomarker_original) %>% 
    mutate(iem_gene_count = length(unique(x.genes)),
           iem_gene_deleted = paste0(unique(x.genes),collapse = ";"),
           iem_rxn_biomarker = paste0(unique(biomarker_id),collapse = ";"),
           iem_rxn_all = paste0(unique(rxn_id), collapse = ";"),
           iem_rxn_direct = paste0(unique(rxn_id[rxn_direct]), collapse = ";"))
}


ef_factor_sort = function(x.break, x.label) {
  x.tbl = data_frame(id = x.break, name = x.label) %>% 
    distinct %>% arrange(id,name) %>%
    group_by(name) %>% slice(1) %>% ungroup %>% arrange(id, name)
  factor(x.label, levels = x.tbl$name, ordered = T)
}


ggcolor = c("non-enzymatic" = "#7F7F7F", "shared" = "#674EA7", 
            "none" = "#7F7F7F", "Unchanged" = "#7F7F7F", #"none" = "#A6A6A6"
            "both" = "#674EA7",
            "rat-specific" = "#C0504D", "rat" = "#C0504D", "rno" = "#C0504D", 
            "human-specific" = "#4F81BD", "human" = "#4F81BD", "hsa" = "#4F81BD", 
            "mmu" = "#B18E5F", 
            "similar" = "#674EA7", "up" = "#674EA7", "Elevated" = "#674EA7", "elevated" = "#674EA7", 
            "dn" = "#F79646", "opposite" = "#F79646", "Reduced" = "#F79646", "reduced" = "#F79646", 
            "ctl" = "pink", "trt" = "green")

ef_timbr_heatmap = function(x.drugs, x.mets, x.compare = timbr.compare, 
                            x.drug.summary = timbr.drug.summary, x.rxn.summary = timbr.rxn.summary,
                            base.size = 24) {
  x.triangle.size = 0.45
  x.heatmap.data = x.compare %>% 
    inner_join(x.drug.summary %>% filter(drug_id %in% x.drugs) %>% 
                 select(drug_id, drug_cor, drug_fdr, drug_pval, drug_index) %>% 
                 arrange(desc(drug_cor), drug_fdr, drug_pval) %>% 
                 mutate(drug_hindex = 1:n())) %>% 
    inner_join(x.rxn.summary %>% filter(met %in% x.mets) %>% 
                 mutate(biomarker = gsub("L-lactate", "lactate", met, fixed = T)) %>% 
                 select(biomarker, rxn_id, rxn_cor, rxn_fdr, rxn_pval, rxn_ave_rno, rxn_ave_hsa, rxn_index) %>% 
                 arrange(desc(rxn_cor), rxn_fdr, rxn_pval) %>%
                 arrange(rxn_ave_rno + rxn_ave_hsa, rxn_fdr, rxn_pval) %>%
                 # arrange(desc(biomarker))%>% 
                 mutate(rxn_hindex = 1:n()))  %>%
    as.tbl %>% ef_df %>% 
    mutate(compound = drug_name) %>% 
    mutate(biomarker = gsub("L-lactate", "lactate", met, fixed = T)) %>% 
    mutate(biomarker = gsub("(R)-3-hydroxybutanoate", "BHB", biomarker, fixed = T)) %>% 
    mutate(biomarker = gsub("prostaglandin E2", "PGE2", biomarker, fixed = T)) %>%
    mutate(biomarker = gsub("chenodeoxycholic acid", "CDCA", biomarker, fixed = T)) %>%
    mutate(biomarker = gsub("ic acid", "ate", biomarker, fixed = T)) %>% 
    mutate(drug_color = ifelse(drug_fdr < 0.1, ifelse(drug_cor > 0, ggcolor["similar"], ggcolor["opposite"]), "#BBBBBB")) %>% 
    mutate(rxn_color = ifelse(rxn_fdr < 0.1, ifelse(rxn_cor > 0, ggcolor["similar"], ggcolor["opposite"]), "#BBBBBB"))
  
  x.heatmap.triangles = x.heatmap.data %>% 
    mutate(poly_hsa = "hsa", poly_rno = "rno") %>%
    arrange(drug_hindex, rxn_hindex) %>%
    left_join(data_frame(x_hsa = c(-1, 1, 1), y_hsa = -c( 1, 1,-1)) %>%
                mutate(poly_hsa = "hsa")) %>%
    left_join(data_frame(x_rno = c(-1,-1, 1), y_rno = -c(-1, 1,-1)) %>%
                mutate(poly_rno = "rno")) %>%
    mutate(x1 = as.numeric(factor(drug_hindex)) + x_hsa * x.triangle.size, 
           y1 = as.numeric(factor(rxn_hindex)) + y_hsa * x.triangle.size,
           x2 = as.numeric(factor(drug_hindex)) + x_rno * x.triangle.size, 
           y2 = as.numeric(factor(rxn_hindex)) + y_rno * x.triangle.size)
  
  x.heatmap.drugs = x.heatmap.data %>% 
    select(compound, drug_id, drug_hindex, drug_color) %>% distinct %>% 
    arrange(drug_hindex) %>% mutate(xx = ef_factor_sort(drug_hindex, compound))
  x.heatmap.rxns = x.heatmap.data %>% 
    select(rxn_id, biomarker, rxn_hindex, rxn_color) %>% distinct %>% 
    arrange(rxn_hindex) %>% mutate(yy = ef_factor_sort(rxn_hindex, biomarker))
  
  heat.max = 1.5
  x.heatmap.plot = x.heatmap.triangles %>% 
    mutate(hsa_fill = pmin(pmax(hsa, -heat.max), heat.max)) %>% 
    mutate(rno_fill = pmin(pmax(rno, -heat.max), heat.max)) %>%
    mutate(xx = ef_factor_sort(drug_hindex, compound)) %>%
    mutate(yy = ef_factor_sort(rxn_hindex, biomarker)) %>%
    ggplot(aes(x = xx, y = yy)) +
    geom_point(alpha = 0) + 
    scale_x_discrete(labels = x.heatmap.drugs$xx, breaks = x.heatmap.drugs$xx) + 
    scale_y_discrete(labels = x.heatmap.rxns$yy, breaks = x.heatmap.rxns$yy) + 
    geom_polygon(aes(x = x1, y = y1, group = interaction(drug_id, rxn_id), fill = hsa_fill), alpha = 1) + 
    geom_polygon(aes(x = x2, y = y2, group = interaction(drug_id, rxn_id), fill = rno_fill), alpha = 1) + 
    scale_fill_gradient2(low = ggcolor["opposite"], mid = "#FFFFFF", 
                         high = ggcolor["similar"], limits = c(-heat.max, heat.max)) + 
    theme_minimal(base_size = base.size) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),#color = x.heatmap.drugs$drug_color, 
          axis.text.y = element_text(angle = 0, hjust = 1),#axis.text.y = element_text(color = x.heatmap.rxns$rxn_color),
          legend.position = "none")  + xlab("") + ylab("")
  x.heatmap.plot
}
