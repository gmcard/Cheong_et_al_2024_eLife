library(neuprintr)
library(malevnc)
#library(gplots)
#library(viridis)
library(ggplot2)
library(dendsort)
library(nat)

#change path to this file as necessary
source(file.path(getwd(), "get_partners_by_multistep_conn_fun.R"))


simple_row_hierarchical_clustering <- function(mat, Rowv = NULL) {
	d=dist(mat,method="euclidean")
	fit = hclust(d, method="ward.D2")
	if (is.null(Rowv)) {
		fit = dendsort(fit)
		return(mat[fit$order,])
	} else {
		fit = reorder(as.dendrogram(fit),Rowv)
		return(mat[labels(fit),])
	}
}

manc_roi_groups = list(NTct = c("NTct(UTct-T1)(L)","NTct(UTct-T1)(R)"), WTct = c("WTct(UTct-T2)(L)","WTct(UTct-T2)(R)"), HTct=c("HTct(UTct-T3)(L)","HTct(UTct-T3)(R)"), IntTct="IntTct", LTct="LTct", 
	LegNp.T1=c("LegNp(T1)(L)", "LegNp(T1)(R)"), LegNp.T2=c("LegNp(T2)(L)", "LegNp(T2)(R)"), LegNp.T3=c("LegNp(T3)(L)", "LegNp(T3)(R)"),
	Ov=c("Ov(L)", "Ov(R)"), ANm="ANm", mVAC=c("mVAC(T1)(L)", "mVAC(T1)(R)", "mVAC(T2)(L)", "mVAC(T2)(R)", "mVAC(T3)(L)", "mVAC(T3)(R)"))

#DN groups that are DNa04, a05 or have similar morphology, to highlight for their role in steering in plots below
all_a04_05 = list(a04=11123, a05=12275, a04_a05_like = c(prior_best_cand_a04=10633, 12155, 12683, 10621, a09=12259))

#subset DNs which have top % input to WTct and HTct
pickRoi=c("WTct", "HTct")
dn_top_percent_wtct = get_neurons_in_roi(unlist(manc_roi_groups[pickRoi]), syn_frac_thres = 0.1, class = "descending neuron", prepost="PRE", by_group = TRUE)
#also needs 100 synapses in target ROI per group
dn_top_percent_wtct_roi_info = neuprint_get_roiInfo(bodyids = dn_top_percent_wtct$bodyid)
dn_top_percent_wtct_roi_info$group = dn_top_percent_wtct$group[match(dn_top_percent_wtct_roi_info$bodyid, dn_top_percent_wtct$bodyid)]
if (length(unlist(manc_roi_groups[pickRoi])>1)) dn_top_percent_wtct_roi_info[[paste0(paste0(pickRoi,collapse='_'),".pre")]] = rowSums(dn_top_percent_wtct_roi_info[paste0(unlist(manc_roi_groups[pickRoi]),".pre")], na.rm = TRUE)
dn_top_percent_wtct_wtct_grp_mean = aggregate(as.formula(paste0(paste0(pickRoi,collapse='_'),".pre ~ group")), data = dn_top_percent_wtct_roi_info, FUN = sum)
dn_top_percent_wtct_roi_info = dn_top_percent_wtct_roi_info[dn_top_percent_wtct_roi_info$group %in% dn_top_percent_wtct_wtct_grp_mean$group[dn_top_percent_wtct_wtct_grp_mean[paste0(paste0(pickRoi,collapse='_'),".pre")] >= 100], ]
dn_top_percent_wtct = dn_top_percent_wtct[dn_top_percent_wtct$group %in% dn_top_percent_wtct_wtct_grp_mean$group[dn_top_percent_wtct_wtct_grp_mean[paste0(paste0(pickRoi,collapse='_'),".pre")] >= 100], ]

#get direct connectivity and indirect conn str of DNs to w-chins of interest and steering MNs
mns_of_interest = c("b1 MN", "b2 MN", "b3 MN", "i1 MN", "i2 MN", "iii1 MN", "iii3 MN", "MNwm35", "hg1 MN", "hg2 MN", "hg3 MN", "hg4 MN", "MNhm42", "MNhm43", "MNhm03")
wchin=c(10073, 10170, 10510, 10667, 10147) #w-chin groups of interest
clio_manc_data = neuprint_list2df(neuprint_fetch_custom(cypher=paste0("MATCH (a:Neuron) WHERE a.group IN [", paste0(wchin, collapse=","), "] OR a.type IN ['", paste0(mns_of_interest, collapse="','"), "'] RETURN a.bodyId AS bodyid, a.class AS class, a.group AS group, a.subclass AS subclass, a.type AS type, a.synonyms AS synonyms, a.rootSide AS root_side, a.somaSide AS soma_side, a.exitNerve AS exit_nerve, a.entry_nerve AS entry_nerve"), timeout=2000))
path_length = 3
clio_manc_data$temp_label = ifelse(clio_manc_data$group %in% wchin,
	paste("w-cHIN", clio_manc_data$group),
	clio_manc_data$type)

#get indirect connectivity upstream of mns/w-chins, then pick out DNs of interest later
mn_us = get_partners_by_multistep_conn(outputids = unique(clio_manc_data$group), path_length = path_length, syn_frac_thres = 0.01, by_group = TRUE, return_type = "detailed")
mn_us_lr = get_partners_by_multistep_conn(outputids = unique(clio_manc_data$group), path_length = path_length, syn_frac_thres = 0, by_group = TRUE, separate_sides = TRUE, return_type = "detailed")

wing_mn_order = c("b1 MN", "b2 MN", "b3 MN", "i1 MN", "i2 MN", "iii1 MN", "iii3 MN", "MNwm35", "hg1 MN", "hg2 MN", "hg3 MN", "hg4 MN", "tp1 MN", "tp2 MN", "tpn MN", "ps1 MN", "ps2 MN", "DLMn a, b", "DLMn c-f", "DVMn 1a-c", "DVMn 2a, b", "DVMn 3a, b", "10178", "46457", "12009", "11372", "11745", "13024", "MNhm42", "MNhm43", "MNhm03")
wing_mn_and_wchin_order = clio_manc_data$group[match(wing_mn_order, clio_manc_data$type)]
wing_mn_and_wchin_order = c(wing_mn_and_wchin_order, wchin)
	
#get meta info for all neurons
manc_all_neurons = neuprint_list2df(neuprint_fetch_custom(cypher=paste0("MATCH (a:Neuron) WHERE a.class IS NOT NULL RETURN a.bodyId AS bodyid, a.class AS class, a.group AS group, a.subclass AS subclass, a.type AS type, a.synonyms AS synonyms, a.rootSide AS root_side, a.somaSide AS soma_side, a.exitNerve AS exit_nerve, a.entry_nerve AS entry_nerve, a.origin AS origin"), timeout=2000))
manc_all_neurons = manc_all_neurons[!(manc_all_neurons$class %in% c("Unkown","Unknown", "TBD", "Glia", "glia")),]
#if is SN, use type as group, as group not always defined--this is also done in the get_partners_by_multistep_conn function
manc_all_neurons$group = ifelse(manc_all_neurons$class %in% c("sensory neuron", "sensory ascending", "sensory descending"), manc_all_neurons$type, manc_all_neurons$group)
#if there are multiple groups per type (i.e. likely subtypes), append group to distinguish them
#this is required for the get_matrix_neuron_names function to work properly!
manc_all_neurons$groups_per_type = sapply(manc_all_neurons$type, FUN=function(x) length(unique(manc_all_neurons$group[manc_all_neurons$type %in% x])))
manc_all_neurons$type = ifelse(manc_all_neurons$groups_per_type>1, 
                               paste(manc_all_neurons$type, manc_all_neurons$group),
                               manc_all_neurons$type)
#manc_all_neurons$type[manc_all_neurons$group %in% wchin] = paste(manc_all_neurons$type[manc_all_neurons$group %in% wchin], manc_all_neurons$group[manc_all_neurons$group %in% wchin])

#function for properly naming neurons in adj matrix for heat area map plotting
#depends on manc_all_neurons being properly retrieved in lines above
get_matrix_neuron_names = function(mat, lr = FALSE) { #lr--is matrix side-separated?
	if (lr) {
		temp_row_group = gsub("_([LRM]|ND)$", "", rownames(mat))
		temp_row_side = ifelse(grepl("_L$", rownames(mat)), "L", ifelse(grepl("_R$", rownames(mat)), "R", "ND"))
		temp_col_group = gsub("_([LRM]|ND)$", "", colnames(mat))
		temp_col_side = ifelse(grepl("_L$", colnames(mat)), "L", ifelse(grepl("_R$", colnames(mat)), "R", "ND"))
	} else {
		temp_row_group = rownames(mat)
		temp_col_group = colnames(mat)
	}
	
	row_names = manc_all_neurons$type[match(temp_row_group, manc_all_neurons$group)]
	col_names = manc_all_neurons$type[match(temp_col_group, manc_all_neurons$group)]
	if(lr) {
		row_names = paste(row_names, temp_row_side, sep = "_")
		col_names = paste(col_names, temp_col_side, sep = "_")
	}
	rownames(mat) = row_names
	colnames(mat) = col_names
	return(mat)
}

heat_area_plot = function(size_mat, color_mat, xlab, ylab, max_point_size = 10) {
	row_order = rownames(size_mat)
	col_order = colnames(size_mat)
	size_mat = as.data.frame(as.matrix(size_mat))
	size_mat$input = rownames(as.matrix(size_mat))
	size_mat = tidyr::pivot_longer(data = as.data.frame(size_mat), cols = -"input", names_to ="output", values_to = "conn_str")
	size_mat$ipsi_contra = apply(size_mat[,c('input','output')], MARGIN = 1, FUN = function(x) {
		if(!(x['input'] %in% rownames(color_mat)) | !(x['output'] %in% colnames(color_mat))) return(NA) else return(color_mat[x['input'],x['output']])
		})
	size_mat$input = factor(size_mat$input, levels = unique(row_order))
	size_mat$output = factor(size_mat$output, levels = unique(col_order))
	
	ggplot(size_mat, aes(x = output, y = input, size = conn_str, color = ipsi_contra)) +
		geom_point(shape = 15, stroke = 0) +
		scale_size_area(max_size = max_point_size) +
		scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0, na.value = 'black') +
		scale_y_discrete(position = "right") +
		scale_x_discrete(guide = guide_axis(angle = 90)) +
		theme_minimal()
}

mn_us_dn = vector(mode = "list", length = path_length)
mn_us_dn_ic = vector(mode = "list", length = path_length)
#make plots
#NOTE THAT HIGHEST i IS ACTUALLY LOWEST PATH LENGTH
for (i in path_length:1) {
	mn_us_dn[[i]] = mn_us[[i]][manc_all_neurons$class[match(rownames(mn_us[[i]]), manc_all_neurons$group)] %in% "descending neuron", ]
	#rearrange by hierarchical clustering
	mn_us_dn[[i]] = simple_row_hierarchical_clustering(mn_us_dn[[i]])
	
	#make ipsi-contra index, for each side, take ipsilateral connectivity minus contralateral over sum
	temp_group = gsub("_([LRM]|ND)$", "", rownames(mn_us_lr[[i]]))
	#check if dn groups have both L and R sides, if not exclude from ipsi-contra analysis
	temp_row_groups_of_interest = sapply(rownames(mn_us_dn[[i]]), FUN = function(x) all(paste0(x, c("_L","_R")) %in% rownames(mn_us_lr[[i]])))

	#reverse side for chINs, as we care about the output side which is contralateral to soma
	temp_chin_groups = gsub("_([LRM]|ND)$", "", colnames(mn_us_lr[[i]]))
	temp_chin_groups_side = ifelse(grepl("_L$", colnames(mn_us_lr[[i]])), "L", ifelse(grepl("_R$", colnames(mn_us_lr[[i]])), "R", NA))
	temp_chin_groups_side[temp_chin_groups %in% wchin] = 
		ifelse(temp_chin_groups_side[temp_chin_groups %in% wchin] == "L", "R", "L")
	temp_mn_us_lr = mn_us_lr[[i]]
	colnames(temp_mn_us_lr) = paste(temp_chin_groups, temp_chin_groups_side, sep = "_")

	#calculate ipsilateral vs contralateral connectivity bias
	#fails if any one group does not have both L and R sides
	temp_ipsi_contra_mat = sapply(names(temp_row_groups_of_interest)[temp_row_groups_of_interest], FUN = function(x) sapply(colnames(mn_us_dn[[i]]), function(y)
			(temp_mn_us_lr[paste0(x,"_L"), paste0(y,"_L")] - temp_mn_us_lr[paste0(x,"_L"), paste0(y,"_R")]
			+ temp_mn_us_lr[paste0(x,"_R"), paste0(y,"_R")] - temp_mn_us_lr[paste0(x,"_R"), paste0(y,"_L")])/
			(temp_mn_us_lr[paste0(x,"_L"), paste0(y,"_L")] + temp_mn_us_lr[paste0(x,"_L"), paste0(y,"_R")]
			+ temp_mn_us_lr[paste0(x,"_R"), paste0(y,"_R")] + temp_mn_us_lr[paste0(x,"_R"), paste0(y,"_L")])
		))
	mn_us_dn_ic[[i]] = t(temp_ipsi_contra_mat)
	
	#reorder rows to put a04/05 and similar at top, reorder columns to specified order
	mn_us_dn[[i]] = mn_us_dn[[i]][order(rownames(mn_us_dn[[i]]) %in% unlist(all_a04_05)),
		order(match(colnames(mn_us_dn[[i]]), wing_mn_and_wchin_order))]
	#color dns that aren't within WTct DN set red
	temp_row_color = ifelse(rownames(mn_us_dn[[i]]) %in% dn_top_percent_wtct$group, "black", "red")
	mn_us_dn[[i]] = get_matrix_neuron_names(mn_us_dn[[i]])
	mn_us_dn_ic[[i]] = get_matrix_neuron_names(mn_us_dn_ic[[i]])
	#png(file=paste0("dn_mn_adj_ipsi_contra_path_len_",i-path_length+1,"_", Sys.Date(),".png"), height=800, width=600)
	pdf(file=paste0("dn_mn_adj_ipsi_contra_path_len_",path_length-i+1,"_", Sys.Date(),".pdf"), height = 800/72, width = 500/72)
	print(heat_area_plot(size_mat = mn_us_dn[[i]], color_mat = mn_us_dn_ic[[i]], xlab = "output", ylab = "input", max_point_size = 8))
	dev.off()
}

