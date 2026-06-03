library(neuprintr)
library(malevnc)
library(igraph)
library(nat)
library(RCy3)

#find neurons that make up an indirect pathway between two neurons of interest, 
#by their contribution to indirect connectivity strength
#uses the get_partners_by_multistep_conn function with inputids and outputids specified
#for a given path length, each neuron (or group) is considered for their contribution to indirect connectivity strength at every possible position in path 
#(e.g. if overall path length is 3, then each neuron is considered at syn hops 1 and 2)
#note that neurons are not guaranteed to form a connected pathway (e.g. a neuron can be slightly above threshold but upstream or downstream neurons are below threshold)
#recommended to be generous (lower) with threshold then manually curate, syn_frac_thres = 0.001 is a decently stringent start

#export graph to cytoscape at end for manual curation and node arrangement, requires cytoscape open

#change path to this file as necessary
source(file.path(getwd(), "get_partners_by_multistep_conn_fun.R"))

#dn groups of interest
all_a04_05 = list(a04=11123, a05=12275, a04_a05_like = c(prior_best_cand_a04=10633, 12155, 12683, 10621, a09=12259))
edge_combine_method = mean #name of a function, used to combine neuron synapse counts (either mean or sum)
#ds neuron groups of interest (MNs and wchins), obtained after inspecting plots from get_partners_by_multistep_conn_fun.R
neurons_of_interest = c(wchin=10073, wchin=10170, wchin=10510, wchin=10667, wchin=10147,
                        b1=10013, b2=10064, b3=10055, iii3=10287, #mns with electrical contact with wchins, even though chemical synapse count is low
                        hm03=10825, hm42=10571, hm43=17216)

#colors
delacroix = c("#C70E7B", "#FC6882", "#007BC3", "#54BCD1",
                       "#009F3F", "#8FDA04", "#AF6125", "#B25D91",
                       "#EF7C12", "#F4B95A", "#C23A4B", "#FBBB48", "#EFEF46", "#31D64D",
                       "#132157","#EE4244", "#D72000", "#1BB6AF", "#8B008B", "#551A8B")
                       
names(delacroix) = c("magenta", "pink", "blue", "cyan", "darkgreen", "green",
                              "brown",
                              "mauve",
                              "darkorange", "orange",
                              "darkred", "darkyellow",
                              "yellow", "palegreen",
                              "navy","cerise",
                              "red", "marine",
                              "purple","darkpurple")

#Neurotransmitter predictions
nt_colors = c("acetylcholine" = delacroix[["darkorange"]],
              "gaba" = delacroix[["blue"]],
              "glutamate" = delacroix[["marine"]],
              "unknown" = delacroix[["darkpurple"]])

#superclass
superclass_colors = c("sensory neuron" = delacroix[["green"]],
                        "sensory ascending"  = delacroix[["darkgreen"]],
                        "descending neuron" = delacroix[["pink"]],
                        "sensory descending" = delacroix[["pink"]],
                        "intrinsic neuron" = delacroix[["cyan"]],
                        "ascending neuron" = delacroix[["blue"]],
                        "motor neuron" = delacroix[["purple"]],
                        "efferent ascending" = delacroix[["brown"]],
                        "efferent neuron" = delacroix[["orange"]])

a04_05_ds_to_mn = get_partners_by_multistep_conn(inputids = unlist(all_a04_05), outputids = neurons_of_interest, path_length=3, syn_frac_thres = 0.001, syn_hop_mult = 0.3, by_group = TRUE, return_type="simple")

#now make the graph
neurons_of_interest = neuprint_list2df(neuprint_fetch_custom(cypher=paste0("MATCH (a:Neuron) WHERE a.group IN [", paste0(a04_05_ds_to_mn, collapse=","), "] RETURN a.bodyId AS bodyid, a.class AS class, a.group AS group, a.subclass AS subclass, a.type AS type, a.synonyms AS synonyms, a.rootSide AS root_side, a.somaSide AS soma_side, a.exitNerve AS exit_nerve, a.entry_nerve AS entry_nerve, a.predictedNt AS pred_nt, a.predictedNtProb AS pred_nt_prob"), timeout=2000))
#append group to w-chins to separate them
neurons_of_interest$type[neurons_of_interest$type %in% "w-cHIN"] = paste(neurons_of_interest$type[neurons_of_interest$type %in% "w-cHIN"], neurons_of_interest$group[neurons_of_interest$type %in% "w-cHIN"])

#consensus nt
neurons_of_interest$pred_nt = c("acetylcholine"="ACh", "gaba"="GABA", "glutamate"="GLU")[neurons_of_interest$pred_nt] #recode
#if multiple nts predicted per group, consider not determined
consensus_nt = aggregate(pred_nt~type, data = neurons_of_interest, FUN = function(x) length(unique(x)))
consensus_nt$consensus_nt = ifelse(consensus_nt$pred_nt > 1, "ND", neurons_of_interest$pred_nt[match(consensus_nt$type, neurons_of_interest$type)])
#if min pred nt prob is below 0.7, consider not determined
consensus_nt$min_prob = sapply(consensus_nt$type, FUN = function(x) min(neurons_of_interest$pred_nt_prob[neurons_of_interest$type==x]))
consensus_nt$consensus_nt = ifelse(consensus_nt$min_prob < 0.7, "ND", consensus_nt$consensus_nt)

neuron_matrix = neuprint_get_adjacency_matrix(neurons_of_interest$bodyid)
temp_types = neurons_of_interest$type[match(colnames(neuron_matrix), neurons_of_interest$bodyid)]

#apply edge combining function of choice
neuron_matrix = t(apply(t(neuron_matrix), 2, function(x) tapply(x, temp_types, edge_combine_method, na.rm = TRUE)))
neuron_matrix = apply(neuron_matrix, 2, function(x) tapply(x, temp_types, edge_combine_method, na.rm = TRUE))

neuron_graph = graph_from_adjacency_matrix(adjmatrix = neuron_matrix, mode = "directed", weighted = TRUE, diag = FALSE)
#set edge symbols and node colors
V(neuron_graph)$nt = consensus_nt$consensus_nt[match(V(neuron_graph)$name, consensus_nt$type)]
V(neuron_graph)$cell_count = sapply(V(neuron_graph)$name, FUN = function(x) sum(neurons_of_interest$type==x))
V(neuron_graph)$display_name = paste0(V(neuron_graph)$name, "(", V(neuron_graph)$cell_count, ")")
E(neuron_graph)$source_nt = V(neuron_graph)$nt[ match(tail_of(neuron_graph, E(neuron_graph))$name, V(neuron_graph)$name) ]
nt_arrowhead_shape = c("ND" = "Circle", ACh = "Delta", GABA = "T", GLU = "Square")
E(neuron_graph)$source_nt_arrow_shape = nt_arrowhead_shape[E(neuron_graph)$source_nt]
V(neuron_graph)$class = neurons_of_interest$class[match(V(neuron_graph)$name, neurons_of_interest$type)]
V(neuron_graph)$class_color = superclass_colors[V(neuron_graph)$class]
E(neuron_graph)$source_class = V(neuron_graph)$class[ match(tail_of(neuron_graph, E(neuron_graph))$name, V(neuron_graph)$name) ]
E(neuron_graph)$source_class_color = superclass_colors[E(neuron_graph)$source_class]
#igraph::V(neuron_graph)$node_shape = node_shape
V(neuron_graph)$node_width = nchar(V(neuron_graph)$display_name)*7

#needs cytoscape open
createNetworkFromIgraph(neuron_graph, title = "a04_05_like_ds")
