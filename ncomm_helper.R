
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
  
  
  
  #print(x.files)
  #list(name = x.name,version = x.version,organism = x.name) %>% 
  #  c(x.info.list)
  
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


ef_timbr_weights = function(gene.expression, gpr.setup,gene.ignore = c(), gene.keep.all = F) {
  gpr.filter = gpr.setup %>% ef_df("value") %>% 
    filter(variable %in% gene.expression[["organism_id"]]) 
  gpr.isozyme.value = ef_timbr_summarize_isozyme(gpr.filter, gene.expression,gene.ignore = gene.ignore, gene.keep.all = gene.keep.all,3)
  gpr.complex.value = ef_timbr_summarize_complex(gpr.filter, gene.expression,gene.ignore = gene.ignore, gene.keep.all = gene.keep.all,3)
  gpr.subunit.value = ef_timbr_summarize_subunit(gpr.filter, gene.expression,gene.ignore = gene.ignore, gene.keep.all = gene.keep.all,3)
  gpr.setup %>% select(rxn_id, organism_id = variable,timbr_weight_default) %>% distinct %>%
    inner_join(gene.expression %>% select(limma_id, organism_id, drug_id, time_id, dose_id) %>% distinct) %>% 
    left_join(bind_rows(list(gpr.isozyme.value,gpr.subunit.value,gpr.complex.value)) %>%
                select(limma_id, organism_id, drug_id, time_id, dose_id, 
                       rxn_id, gpr_significance, gpr_logfc)) %>%
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
    x.tbl = x.gpr  %>% 
      left_join(x.gene %>% mutate(gene_id = as.character(gene_id)), 
                by = c("organism_id","gene_id")) %>%
      mutate(fdr = ifelse(!is.na(fdr),fdr,1),
             logfc = ifelse(!is.na(logfc),logfc,0))
  } else {
    x.tbl = x.gpr %>% 
      inner_join(x.gene %>% mutate(gene_id = as.character(gene_id)), 
                 by = c("organism_id","gene_id"))
  }
  if (!is.null(x.tbl$logfc)) {
    x.tbl %>% mutate(mag = sign(logfc) * -log10(fdr))
  }
  x.tbl
}


ef_timbr_summarize_subunit = function(gpr.info, gene.value, gene.ignore, gene.keep.all = F, max.value = 3) {
  gpr.subunit.split = gpr.info %>% filter(grepl("\\(",value)) %>%
    mutate(gpr_rule = gsub("\\;","|",gsub("\\:","&",value))) %>%
    ef_split_value("\\:") %>% ef_df
  gpr.subunit.unique = gpr.subunit.split %>% select(variable,value) %>% distinct %>% 
    mutate(subunit_id = paste0("subunit",1:n()), 
           subunit_rule = gsub("[ \\(\\)]","",gsub("\\;"," | ",value)))
  gpr.subunit.breakdown = gpr.subunit.split %>% 
    left_join(gpr.subunit.unique, by = c("variable","value")) %>% 
    mutate(value = gsub("\\)$","",gsub("^\\(","",value))) %>% 
    ef_split_value("\\;") %>% ef_df("value") %>% filter(value != "0") %>% 
    rename(organism_id = variable, gene_id = value) %>%
    mutate(gene_ignore = gene_id %in% gene.ignore) %>%
    group_by(rxn_id, organism_id) %>% 
    mutate(n_gene = length(unique(gene_id)),
           n_ignore = length(unique(gene_id[gene_ignore])),
           n_subunit = length(unique(subunit_rule))) %>% ungroup
  gpr.subunit.breakdown %>% ef_timbr_gpr_join(gene.value, gene.keep.all) %>% 
    mutate(mag = -log10(fdr) * sign(logfc)) %>%
    arrange(fdr, desc(abs(logfc))) %>% 
    group_by(limma_id, organism_id, drug_id, time_id, dose_id, 
             rxn_id, n_gene, n_ignore, n_subunit, subunit_id, subunit_rule) %>% 
    summarize(gpr_significance = mean(mag[!gene_ignore], na.rm = T),
              gpr_logfc = mean(logfc[!gene_ignore], na.rm = T)) %>% ungroup %>%
    arrange(desc(abs(gpr_significance)),desc(abs(gpr_logfc))) %>% 
    group_by(limma_id, organism_id, drug_id, time_id, dose_id, 
             rxn_id, n_gene, n_ignore, n_subunit) %>% 
    slice(1) %>% ungroup %>%
    mutate(gpr_logfc = pmax(pmin(ifelse(!is.na(gpr_logfc),gpr_logfc, 0), 
                                 max.value), -max.value)) %>%
    mutate(gpr_significance = pmax(pmin(ifelse(!is.na(gpr_significance),gpr_significance, 0), 
                                        max.value), -max.value))
}

ef_timbr_summarize_complex = function(gpr.info, gene.value, gene.ignore,
                                      gene.keep.all = F, max.value = 3) {
  
  gpr.complex.split = gpr.info %>% 
    filter(grepl("\\:",value), !grepl("\\(",value)) %>% 
    mutate(gpr_rule = gsub("\\;","|",gsub("\\:","&",value))) %>%
    ef_split_value(";") %>% ef_df
  gpr.complex.unique = gpr.complex.split %>% select(variable,value) %>% distinct %>% 
    mutate(complex_id = paste0("complex",1:n()), complex_rule = gsub("\\:"," & ",value))
  gpr.complex.breakdown = gpr.complex.split %>% 
    left_join(gpr.complex.unique, by = c("variable","value")) %>% 
    ef_split_value("\\:") %>% ef_df("value") %>% filter(value != "0") %>%
    rename(organism_id = variable, gene_id = value) %>% 
    mutate(gene_ignore = gene_id %in% gene.ignore) %>%
    group_by(rxn_id, organism_id) %>% 
    mutate(n_gene = length(unique(gene_id)),
           n_ignore = length(unique(gene_id[gene_ignore])),
           n_complex = length(unique(complex_rule))) %>% ungroup %>%
    mutate(n_subunit = n_gene)
  gpr.complex.breakdown %>% ef_timbr_gpr_join(gene.value, gene.keep.all) %>% 
    mutate(mag = -log10(fdr) * sign(logfc)) %>%
    arrange(fdr, desc(abs(logfc))) %>% 
    group_by(limma_id, organism_id, drug_id, time_id, dose_id, 
             rxn_id, n_gene, n_ignore, n_subunit, n_complex, complex_id, complex_rule) %>% 
    slice(1) %>% ungroup %>% 
    group_by(limma_id, organism_id, drug_id, time_id, dose_id, 
             rxn_id, n_gene, n_ignore, n_subunit, n_complex) %>% 
    summarize(gpr_significance = mean(mag[!gene_ignore], na.rm = T),
              gpr_logfc = mean(logfc[!gene_ignore], na.rm = T)) %>% ungroup %>%
    mutate(gpr_logfc = pmax(pmin(ifelse(!is.na(gpr_logfc), gpr_logfc, 0), 
                                 max.value), -max.value)) %>%
    mutate(gpr_significance = pmax(pmin(ifelse(!is.na(gpr_significance),gpr_significance, 0), 
                                        max.value), -max.value))
}

ef_timbr_summarize_isozyme = function(gpr.info, gene.value, gene.ignore, gene.keep.all = F, max.value = 3) {
  gpr.isozyme.split = gpr.info %>% 
    filter(!grepl("\\:",value), !grepl("\\(",value)) %>% 
    mutate(gpr_rule = gsub("\\;","|",gsub("\\:","&",value))) %>%
    ef_split_value(";") %>% ef_df
  gpr.isozyme.breakdown = gpr.isozyme.split %>% ef_df("value") %>% filter(value != "0")  %>% 
    rename(organism_id = variable, gene_id = value) %>%
    mutate(gene_ignore = gene_id %in% gene.ignore) %>%
    group_by(rxn_id,organism_id) %>% 
    mutate(n_gene = length(unique(gene_id)),
           n_ignore = length(unique(gene_id[gene_ignore])),
           n_subunit = 1) %>% ungroup
  gpr.isozyme.breakdown %>% ef_timbr_gpr_join(gene.value, x.keep = gene.keep.all) %>% 
    mutate(mag = -log10(fdr) * sign(logfc)) %>%
    group_by(limma_id, organism_id, drug_id, time_id, dose_id, 
             rxn_id, n_gene, n_ignore, n_subunit) %>% 
    summarize(gpr_significance = mean(mag[!gene_ignore], na.rm = T),
              gpr_logfc = mean(logfc[!gene_ignore], na.rm = T)) %>% ungroup %>%
    mutate(gpr_logfc = pmax(pmin(ifelse(!is.na(gpr_logfc), gpr_logfc, 0), max.value), -max.value)) %>%
    mutate(gpr_significance = pmax(pmin(ifelse(!is.na(gpr_significance), gpr_significance, 0), 
                                        max.value), -max.value))
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



