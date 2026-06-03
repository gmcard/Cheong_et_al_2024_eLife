neuprint_get_adjacency_matrix_memo <- memoise::memoise(neuprintr::neuprint_get_adjacency_matrix, ~memoise::timeout(168*60*60))

combine_matrix_rows_cols_by_group_memo <- memoise::memoise(function(adj, groups=NA, row_groups=NA, col_groups=NA, sparse = FALSE) {
	if (all(is.na(groups)) & (all(is.na(row_groups))|all(is.na(col_groups)))) stop("Either groups or both row_groups and col_groups need to be specified")
	if (all(is.na(row_groups))|all(is.na(col_groups))) i = j = as.character(groups) 
	else {
		i = as.character(row_groups)
		j = as.character(col_groups)
	}
	i[is.na(i)] = 'NA'
	j[is.na(j)] = 'NA'
	i = Matrix::sparse.model.matrix(~i+0,transpose=T)
	rownames(i) = substr(rownames(i),2,1000)
	j = Matrix::sparse.model.matrix(~j+0)
	colnames(j) = substr(colnames(j),2,1000)
	adj = i %*% adj %*% j
	if (sparse) adj else as.matrix(adj)
}, ~memoise::timeout(168*60*60))


#get indirect upstream or downstream partners of specified path length (or below) away from given bodyids or groups by matrix multiplication of input fractions (as in e.g. Li et al. 2020) 
#use syn_frac_thres to subset top partners; as values drop exponentially with increasing path length, threshold is likewise reduced by multiplication by syn_hop_mult^(path_length-1)
#use by_group = TRUE to check by group-to-group connectivity
#there are options to remove dn or sn upstream connections (remove_dn_us, remove_sn_us) as input fractions are inaccurate for neurons which have input outside volume
#return_type can be "simple" or "detailed"--simple returns groups or bodyids named by path length, detailed returns a list of data frames with values of indirect_conn_str--the most upstream neurons are placed at the lowest list index
#if both inputids and outputids are specified, find top 'in between' neurons that contribute to the indirect connectivity strength values--for a given path length, each neuron (or group) is considered for their contribution to indirect connectivity strength at every possible position in path 
#(e.g. if overall path length is 3, then each neuron is considered at syn hops 1 and 2)
get_partners_by_multistep_conn <- function(inputids = NULL, outputids = NULL, path_length, syn_frac_thres = 0.01, syn_hop_mult = 0.3, by_group = FALSE, remove_dn_us = TRUE, remove_sn_us = TRUE, return_type = "simple", separate_sides = FALSE) {
	if (is.null(inputids) & is.null(outputids)) stop("need to supply inputids or outputids (or both)")
	if (!is.numeric(syn_frac_thres) | syn_frac_thres > 1 | syn_frac_thres < 0 | length(syn_frac_thres) != 1) stop("syn_frac_thres needs to be a single number between 0 and 1")
	if (!is.numeric(syn_hop_mult) | syn_hop_mult > 1 | syn_hop_mult < 0 | length(syn_hop_mult) != 1) stop("syn_hop_mult needs to be a single number between 0 and 1")
	if (is.na(path_length) | path_length%%1 != 0 | path_length < 1) stop("path_length needs to be integer greater than zero")
	if (!(return_type %in% c("simple","detailed"))) stop("return_type must be 'simple' or 'detailed'")
	if (!is.logical(separate_sides)) stop("separate_sides must be TRUE/FALSE")
	
	#only consider "real" neurons: us or ds connectivity needs to be more than 25, needs to have a valid status that isn't orphan, need to have a class that isn't unknown or glia
	#changed criteria as of 1-19-23, now only considers bodies annotated with a neuron class
	#clio_manc_data = neuprintr::neuprint_list2df(neuprintr::neuprint_fetch_custom(cypher="MATCH (a:Neuron) WHERE (a.upstream > 25 OR a.downstream > 25) AND a.status IN ['Traced','Sensory Anchor','Soma Anchor','Anchor','Primary Anchor','Cervical Anchor'] RETURN a.bodyId AS bodyid, a.class AS class, a.group AS group, a.subclass AS subclass, a.rootSide AS root_side, a.somaSide AS soma_side, a.exitNerve AS exit_nerve, a.entry_nerve AS entry_nerve", timeout=2000))
	clio_manc_data = neuprintr::neuprint_list2df(neuprintr::neuprint_fetch_custom(cypher=paste0("MATCH (a:Neuron) WHERE a.class IS NOT NULL RETURN a.bodyId AS bodyid, a.class AS class, a.group AS group, a.subclass AS subclass, a.type AS type, a.synonyms AS synonyms, a.rootSide AS root_side, a.somaSide AS soma_side, a.exitNerve AS exit_nerve, a.entry_nerve AS entry_nerve"), timeout=2000))
	clio_manc_data$class[is.na(clio_manc_data$class)] = "ND"
	clio_manc_data = clio_manc_data[!(clio_manc_data$class %in% c("Unkown","Unknown", "TBD", "Glia", "glia", "ND")),]
	#clio_manc_data$abbrv = class2abbrv(clio_manc_data$class)
	#if is SN, use type as group, as group not always defined
	clio_manc_data$group = ifelse(clio_manc_data$class %in% c("sensory neuron", "sensory ascending", "sensory descending"), clio_manc_data$type, clio_manc_data$group)
	#if no group, make bodyid group
	clio_manc_data$group = ifelse(is.na(clio_manc_data$group), clio_manc_data$bodyid, clio_manc_data$group)
	
	#sanity check that bodyids/groups given in function arguments are actually in neuprint
	if (!is.null(inputids)) {
		if(by_group) {
			if (!any(inputids %in% clio_manc_data$group)) stop("group(s) ", inputids[!(inputids %in% clio_manc_data$group)], " in inputids not found in neuprint groups")
			if (!all(inputids %in% clio_manc_data$group)) {
				warning("some inputids (by_group = TRUE) not found in neuprint groups. Removing.\n")
				inputids = inputids[inputids %in% clio_manc_data$group]
			}
		}
		if(!by_group & !any(inputids %in% clio_manc_data$bodyid)) stop("inputid(s) ", inputids[!(inputids %in% clio_manc_data$bodyid)], " not found in neuprint groups")
	}
	if (!is.null(outputids)) {
		if(by_group) {
			if (!any(outputids %in% clio_manc_data$group)) stop("group(s) ", outputids[!(outputids %in% clio_manc_data$group)], " in outputids not found in neuprint groups")
			if (!all(outputids %in% clio_manc_data$group)) {
				warning("some outputids (by_group = TRUE) not found in neuprint groups. Removing.\n")
				outputids = outputids[outputids %in% clio_manc_data$group]
			}
		}
		if(!by_group & !any(outputids %in% clio_manc_data$bodyid)) stop("outputid(s) ", outputids[!(outputids %in% clio_manc_data$bodyid)], " not found in neuprint groups")
	}
	
	
	#if separate_sides is TRUE, append root/entry/exit/soma side to group to treat them as separate groups
	#separate_sides is only meaningful when by_group is TRUE
	get_nerve_side <- function(nerve) ifelse(grepl("_L$", nerve), "LHS", ifelse(grepl("_R$", nerve), "RHS", NA))
	if (separate_sides & by_group) {
		clio_manc_data$exit_nerve = get_nerve_side(clio_manc_data$exit_nerve)
		clio_manc_data$entry_nerve = get_nerve_side(clio_manc_data$entry_nerve)
		clio_manc_data$side = ifelse(!is.na(clio_manc_data$root_side), clio_manc_data$root_side,
			ifelse(!is.na(clio_manc_data$exit_nerve), clio_manc_data$exit_nerve,
				ifelse(!is.na(clio_manc_data$entry_nerve), clio_manc_data$entry_nerve,
					ifelse(!is.na(clio_manc_data$soma_side), clio_manc_data$soma_side,
					"ND")
				)
			)
		)
		clio_manc_data$side = sapply(clio_manc_data$side, FUN = function(x) switch(x, LHS = "L", RHS = "R", Midline = "M", ND = "ND", "ND"))
		clio_manc_data$original.group = clio_manc_data$group
		clio_manc_data$group = paste(clio_manc_data$group, clio_manc_data$side, sep = "_")

		if (!is.null(inputids)) inputids = unique(clio_manc_data$group[clio_manc_data$original.group %in% inputids])
		if (!is.null(outputids)) outputids = unique(clio_manc_data$group[clio_manc_data$original.group %in% outputids])
	}


	adj_all_vnc = neuprint_get_adjacency_matrix_memo(bodyids = clio_manc_data$bodyid, chunksize = 1000L, timeout=2000, sparse = TRUE)
	adj_all_vnc = adj_all_vnc[Matrix::colSums(adj_all_vnc)>0, Matrix::colSums(adj_all_vnc)>0]
	#combine by group
	if (by_group) {
		adj_all_groups = clio_manc_data$group[match(colnames(adj_all_vnc), clio_manc_data$bodyid)]
		adj_all_vnc = combine_matrix_rows_cols_by_group_memo(adj_all_vnc, adj_all_groups, sparse = TRUE)
		adj_all_groups = rownames(adj_all_vnc)
		adj_all_abbrv = clio_manc_data$class[match(adj_all_groups, clio_manc_data$group)]
	} else {
		adj_all_abbrv = clio_manc_data$class[match(rownames(adj_all_vnc), clio_manc_data$bodyid)]
	}

	#convert to input fractions
	adj_all_vnc = adj_all_vnc %*% Matrix::Diagonal(x = 1/Matrix::colSums(adj_all_vnc))
	#why do colnames go missing here? reassign
	colnames(adj_all_vnc) = rownames(adj_all_vnc)

	if (remove_dn_us) {
		#input fractions are inaccurate for DNs and SNs, since most input comes from outside volume
		#set input of these neurons to zero for now--if we want to find inputs to these neurons, we need to consider them in a separate analysis
		adj_all_vnc[, adj_all_abbrv %in% "descending neuron"] = 0
	}
	if (remove_sn_us) adj_all_vnc[, adj_all_abbrv %in% c("sensory neuron", "sensory ascending", "sensory descending")] = 0
	#set diagonal (within group synapses) to zero
	diag(adj_all_vnc) = 0


	layer_list_no_thres = vector(mode = "list", length = path_length)
	layer_list_no_thres[[1]] = adj_all_vnc
	#get matrix for path length 2 onwards--for each syn hop, multiply input fraction matrix with matrix from current path length - 1
	if (path_length > 1) for (i in 2:path_length) {
		layer_list_no_thres[[i]] = layer_list_no_thres[[i-1]]%*% adj_all_vnc
	}
	
	#now get partners above threshold
	if (!is.null(inputids) & !is.null(outputids)) {
		inputids = as.character(unique(inputids))
		outputids = as.character(unique(outputids))
		# path_length has to be more than 1
		if (path_length == 1) stop("if inputids and outputids are supplied, path length has to be >1")
		layer_list = vector(mode = "list", length = path_length)
		for (curr_path_length in 2:path_length) {
			for (i in 1:(curr_path_length - 1)) {
				#loop through each candidate neuron
				#this would be more efficient as an apply
				temp_plot_list = c()
				for (j in inputids) {
					for (k in outputids) {
						temp_plot_list = (layer_list_no_thres[[i]][, k] * t(layer_list_no_thres[[curr_path_length - i]][j, ]) > syn_frac_thres*syn_hop_mult^(curr_path_length-1))
						temp_plot_list = colnames(layer_list_no_thres[[i]])[as.vector(temp_plot_list)]
						layer_list[[i]] = c(temp_plot_list, layer_list[[i]])
					}
				}
				layer_list[[i]] = unique(layer_list[[i]])
			}
		}
		layer_list[[path_length]] = outputids
		
		bodyids_above_thres = inputids
		names(bodyids_above_thres) = rep(0, length(bodyids_above_thres))
		for (i in 1:path_length) {
			layer_list[[i]] = layer_list[[i]][!(layer_list[[i]] %in% bodyids_above_thres)]
			if (length(layer_list[[i]]) == 0) warning("partners not found at path length ", i, ". Try lowering syn_frac_thres or syn_hop_mult")
			names(layer_list[[i]]) = rep(i, length(layer_list[[i]]))
			bodyids_above_thres = c(bodyids_above_thres, layer_list[[i]])
		}
		
		
		
	} else if (!is.null(inputids)) {
		inputids = as.character(unique(inputids))
		layer_list = vector(mode = "list", length = path_length)
		#get direct upstream or ds--simple frac threshold
		for (i in 1:path_length) {
			neurons_above_thres = apply(layer_list_no_thres[[i]][inputids, , drop = FALSE] >= syn_frac_thres*syn_hop_mult^(i-1), MARGIN = 2, FUN = any)
			if (return_type == 'simple') layer_list[[i]] = colnames(layer_list_no_thres[[i]])[neurons_above_thres] else
				layer_list[[i]] = layer_list_no_thres[[i]][inputids, neurons_above_thres, drop = FALSE]
		}
		
		#if simple output, place bodyids above thres in the lowest path length that they appear in (closest to inputids)
		if (return_type == 'simple') {
			bodyids_above_thres = inputids
			names(bodyids_above_thres) = rep(0, length(bodyids_above_thres))
			for (i in 1:path_length) {
				layer_list[[i]] = layer_list[[i]][!(layer_list[[i]] %in% bodyids_above_thres)]
				if (length(layer_list[[i]]) == 0) warning("partners not found at path length ", i, ". Try lowering syn_frac_thres or syn_hop_mult")
				names(layer_list[[i]]) = rep(i, length(layer_list[[i]]))
				bodyids_above_thres = c(bodyids_above_thres, layer_list[[i]])
			}
		}
		
		
	} else if (!is.null(outputids)) {
		outputids = as.character(unique(outputids))
		layer_list = vector(mode = "list", length = path_length + 1)
		#get direct upstream or ds--simple frac threshold
		for (i in 1:path_length) {
			neurons_above_thres = apply(layer_list_no_thres[[i]][, outputids, drop = FALSE] >= syn_frac_thres*syn_hop_mult^(i-1), MARGIN = 1, FUN = any)
			if (return_type == 'simple') layer_list[[path_length - i + 1]] = colnames(layer_list_no_thres[[i]])[neurons_above_thres] else
				layer_list[[path_length - i + 1]] = layer_list_no_thres[[i]][neurons_above_thres, outputids, drop = FALSE]
		}
		
		#if simple output, place bodyids above thres in the highest path length that they appear in (closest to outputids)
		if (return_type == 'simple') {
			bodyids_above_thres = outputids
			names(bodyids_above_thres) = rep(path_length, length(bodyids_above_thres))
			for (i in path_length:1) {
				layer_list[[i]] = layer_list[[i]][!(layer_list[[i]] %in% bodyids_above_thres)]
				if (length(layer_list[[i]]) == 0) warning("partners not found at path length ", i, ". Try lowering syn_frac_thres or syn_hop_mult")
				names(layer_list[[i]]) = rep(i - 1, length(layer_list[[i]]))
				bodyids_above_thres = c(layer_list[[i]], bodyids_above_thres)
			}
		}
	}
	
	
	#prepare to return data
	if (by_group & return_type == 'simple') {
		#convert group back to bodyids
		temp_bodyids_above_thres = clio_manc_data$bodyid[clio_manc_data$group %in% bodyids_above_thres]
		temp_bodyids_above_thres_groups = clio_manc_data$group[match(temp_bodyids_above_thres, clio_manc_data$bodyid)]
		names(temp_bodyids_above_thres) = names(bodyids_above_thres)[match(temp_bodyids_above_thres_groups, bodyids_above_thres)]
		bodyids_above_thres = temp_bodyids_above_thres
	}
	
	if (return_type == "simple") return(bodyids_above_thres) else {
		if (!is.null(inputids) & !is.null(outputids)) {
			warning("code for return type 'detailed' when inputids and outputids are defined is not yet written.\nReturning simple output\n")
			return(bodyids_above_thres)
		} else {
			return(layer_list)
		}
	}
}




#roi is a vector of rois to search for. If multiple, then sum of rois must be equal or larger than syn_frac_thres
#class is a vector of neuron types (as defined in neuprint class) to search for.
get_neurons_in_roi <- function(roi, syn_frac_thres, class = NA, prepost, by_group = FALSE, ...) {
	if (!is.numeric(syn_frac_thres)) stop("syn_frac_thres needs to be numeric")
	if (syn_frac_thres <= 0 | syn_frac_thres > 1 ) stop("syn_frac_thres needs to be between 0 and 1")
	if (!is.logical(by_group)) stop("by_group needs to be TRUE/FALSE")
	if (!(prepost %in% c("PRE", "POST"))) stop("prepost should be either 'PRE' or 'POST'")
	
	prepost = ifelse(prepost == "PRE", "pre", "post")
	
	roi_enum = letters[1:length(roi)]
	temp_roi_string = paste0("apoc.convert.fromJsonMap(a.roiInfo)['", roi, "'].", prepost," AS ", roi_enum, collapse = ",")
	if (all(is.na(class))) {
		roi_info = neuprintr::neuprint_list2df(neuprintr::neuprint_fetch_custom(cypher=paste0("MATCH (a:Neuron) WHERE (a.upstream > 25 OR a.downstream > 25) AND a.class IS NOT NULL RETURN a.bodyId as bodyid, a.group as group, a.", prepost, " AS prepost,", temp_roi_string), timeout=2000, ...))
		roi_info = roi_info[!(roi_info$class %in% c("Unkown","Unknown", "TBD", "Glia", "glia")),]
	} else
		roi_info = neuprintr::neuprint_list2df(neuprintr::neuprint_fetch_custom(cypher=paste0("MATCH (a:Neuron) WHERE (a.upstream > 25 OR a.downstream > 25) AND a.class IN ['", paste0(class,collapse="','"), "'] RETURN a.bodyId as bodyid, a.group as group, a.", prepost, " AS prepost,", temp_roi_string), timeout=2000, ...))
	roi_info[roi_enum][is.na(roi_info[roi_enum])] = 0
	roi_info$roi_over_sum = rowSums(roi_info[roi_enum])/roi_info$prepost
	roi_info$roi_over_sum[is.na(roi_info$roi_over_sum)] = 0
	
	if (!by_group) return(roi_info[roi_info$roi_over_sum >= syn_frac_thres, c("bodyid", "group", "roi_over_sum")])
	else {
		group_avgs = stats::aggregate(roi_over_sum ~group, data = roi_info, FUN = mean)
		return(roi_info[roi_info$group %in% group_avgs$group[group_avgs$roi_over_sum >= syn_frac_thres],  c("bodyid", "group", "roi_over_sum")])
	}
}