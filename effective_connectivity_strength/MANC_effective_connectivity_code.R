library(neuprintr)
library(malevnc)
library(bit64)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Matrix)
library(cowplot)

choose_malevnc_dataset("MANC")
#get the whole maleVNC matrix, these first steps will take quite some time
mdf = manc_dvid_annotations()
query.ids = as.integer64(mdf$bodyid[!(mdf$status %in% c("Orphan","PRT Orphan", "RT Orphan", "Unimportant",NA, ""))])
query_info = manc_neuprint_meta(query.ids)
query_info = query_info[!(query_info$class %in% c("TBD", "glia",NA)),]
query.ids2 = na.omit(as.integer64(query_info$bodyid))
np_adj = neuprint_get_adjacency_matrix(query.ids2, sparse=T)

colScale <- function(A, na.rm=TRUE) {
  scalefac = 1 / Matrix::colSums(A)
  if(na.rm) scalefac[!is.finite(scalefac)]=0
  B = A %*% Matrix::Diagonal(x = scalefac)
  B
}

#now calculate the input percent using this funtion
np_adj_per_in <- colScale(np_adj)
colnames(np_adj_per_in) <- colnames(np_adj) 
#set the threshold and number of layers
desired_layers = 10

#get your input IDs and info from clio
all_info = manc_body_annotations()
#for example some query or a list of IDs (example here descending neurons of the front legs)
# input_ids <- c(13491, 14742, 14844, 15695) #using specific IDs
input_ids <- na.omit(all_info$bodyid[all_info$subclass=="fl"&all_info$class=="descending neuron"])
input_info = all_info[all_info$bodyid %in% input_ids,]
input_info_sel <- input_info[,c("bodyid","class","group","serial","type","systematic_type","soma_side","root_side","origin","target","predicted_nt","predicted_nt_prob")]
input_info_sel$soma_or_root_side = ifelse(!is.na(input_info_sel$soma_side)&input_info_sel$soma_side!="",input_info_sel$soma_side,input_info_sel$root_side)
input_info_sel %>% mutate(soma_or_root_side=recode(soma_or_root_side,"LHS" = "L","RHS" = "R","Midline" = "M")) ->input_info_sel

#get your output IDs and info from clio
#for example some query or a list of IDs (example here motor neurons of the T1 leg)
# output_ids <- c(10090, 10104, 10360, 10379, 10442, 10460, 10543, 10672, 10760, 10837, 10884, 10924, 10971, 10975, 11002, 11066, 11233, 11258, 11314, 11392, 11563, 11610, 11625, 11739) #using specific IDs
output_ids <- na.omit(all_info$bodyid[all_info$subclass=="fl"&all_info$class=="motor neuron"])
output_info = all_info[all_info$bodyid %in% output_ids,]
output_info_sel <- output_info[,c("bodyid","class","group","serial","type","systematic_type","soma_side","root_side","origin","target","predicted_nt","predicted_nt_prob")]
output_info_sel$soma_or_root_side = ifelse(!is.na(output_info_sel$soma_side)&output_info_sel$soma_side!="",output_info_sel$soma_side,output_info_sel$root_side)
output_info_sel %>% mutate(soma_or_root_side=recode(soma_or_root_side,"LHS" = "L","RHS" = "R","Midline" = "M")) ->output_info_sel

#precalculate all layers for all input neurons, this step will take a while depending on how many input neurons you are using
#set up a list of vectors for each iteration through the network and make the first vector considering only your starting neurons
layer_vec = paste0("v",seq(from=0,to=desired_layers))
layer_list <- vector("list", length(layer_vec))
names(layer_list) <- layer_vec
layer_list[["v1"]] = np_adj_per_in[rownames(np_adj_per_in) %in% input_info_sel$bodyid,]
for (layer in 2:desired_layers) {
  layer_list[[layer+1]] <-  layer_list[[layer]] %*% np_adj_per_in
}

#now pull out data from multiplied matrices for input neurons, this can be based on a lot of things, choosing type as example
input_types <- unique(input_info$type)
for (s in 1:length(input_types)) {
  assign(input_types[s], na.omit(unique(input_info$bodyid[grepl(input_types[s],input_info$type)])))
}
#for output neurons, this can be based on a lot of things, choosing type as example
output_order <- sort(unique(output_info$type))


in_out_scores <- setNames(data.frame(matrix(ncol = 18, nrow = 0)),c("bodyid","conn_score","ds_layer","class","group","serial","type","systematic_type",
                                                                    "soma_side","root_side","origin","target","predicted_nt","predicted_nt_prob",
                                                                    "soma_or_root_side","input","in_side","ipsi_contra"))
for (n in 1:length(input_types)) {
  startn_r = as.character(na.omit(input_info_sel$bodyid[input_info_sel$type == input_types[n] & input_info_sel$soma_or_root_side == "R"]))
  startn_l = as.character(na.omit(input_info_sel$bodyid[input_info_sel$type == input_types[n] & input_info_sel$soma_or_root_side == "L"]))
  if (length(startn_r)>0 & length(startn_l)>0) {
    input = input_types[n]
    
    #set up a list of vectors for each iteration through the network
    norm_list_l <- vector("list", length(layer_vec))
    names(norm_list_l) <- layer_vec
    norm_list_r <- vector("list", length(layer_vec))
    names(norm_list_r) <- layer_vec
    #get data from precalculated matrices and normalize
    for (layer in 1:desired_layers) {
      if (length(startn_l)>1) {
        temp_layer_conn_str = Matrix::colSums(layer_list[[layer+1]][na.omit(match(startn_l, rownames(layer_list[[layer+1]]))),])
      } else {
        temp_layer_conn_str = layer_list[[layer+1]][na.omit(match(startn_l, rownames(layer_list[[layer+1]]))),]  
      }
      norm_list_l[[layer+1]] <- setNames(as.vector(temp_layer_conn_str/mean(temp_layer_conn_str[temp_layer_conn_str>0])),
                                         colnames(np_adj_per_in))
      if (length(startn_r)>1) {
        temp_layer_conn_str = Matrix::colSums(layer_list[[layer+1]][na.omit(match(startn_r, rownames(layer_list[[layer+1]]))),])
      } else {
        temp_layer_conn_str = layer_list[[layer+1]][na.omit(match(startn_r, rownames(layer_list[[layer+1]]))),]  
      }
      norm_list_r[[layer+1]] <- setNames(as.vector(temp_layer_conn_str/mean(temp_layer_conn_str[temp_layer_conn_str>0])),
                                         colnames(np_adj_per_in))
    }
    
    #now get the highest connectivity score for the input to each output type
    out_r_scores_df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),c("bodyid","conn_score","ds_layer"))
    for (layer in 2:desired_layers) {
      out_r_scores <- norm_list_r[[layer]][names(norm_list_r[[layer]])%in% output_ids]
      df = data.frame(bodyid=as.double(names(out_r_scores)), conn_score=out_r_scores, ds_layer = layer-1, row.names=NULL)
      out_r_scores_df = rbind(out_r_scores_df, df)
    }
    out_r_scores_df %>%
      group_by(bodyid) %>%
      dplyr::slice(which.max(conn_score)) -> out_r_scores_df_max
    
    out_r_scores_df_max_info <- left_join(out_r_scores_df_max, output_info_sel, by = "bodyid")
    out_r_scores_df_max_info$input <- input
    out_r_scores_df_max_info$in_side <- "R"
    
    out_l_scores_df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),c("bodyid","conn_score","ds_layer"))
    for (layer in 2:desired_layers) {
      out_l_scores <- norm_list_l[[layer]][names(norm_list_l[[layer]])%in% output_ids]
      df = data.frame(bodyid=as.double(names(out_l_scores)), conn_score=out_l_scores, ds_layer = layer-1, row.names=NULL)
      out_l_scores_df = rbind(out_l_scores_df, df)
    }
    out_l_scores_df %>%
      group_by(bodyid) %>%
      dplyr::slice(which.max(conn_score)) -> out_l_scores_df_max
    
    out_l_scores_df_max_info <- left_join(out_l_scores_df_max, output_info_sel, by = "bodyid")
    out_l_scores_df_max_info$input <- input
    out_l_scores_df_max_info$in_side <- "L"
    
    out_scores <- rbind(out_r_scores_df_max_info,out_l_scores_df_max_info)
    out_scores$ipsi_contra <- NA
    out_scores$ipsi_contra[out_scores$soma_or_root_side == out_scores$in_side] <- "ipsi"
    out_scores$ipsi_contra[out_scores$soma_or_root_side != out_scores$in_side] <- "contra"
    
    in_out_scores <- rbind(in_out_scores,out_scores)
  }
}

in_out_scores %>% group_by(input,type,ipsi_contra) %>%
  summarise_at(vars("conn_score", "ds_layer"), mean) -> in_out_scores_mean
#cluster result on ipsi connectivity
in_out_scores_mean %>% filter(ipsi_contra=="ipsi") -> in_out_scores_mean_ipsi
in_out_scores_mean_ipsi_m <- tidyr::pivot_wider(in_out_scores_mean_ipsi[,c("input","conn_score","type")], names_from = "input", values_from = "conn_score")
in_out_scores_mean_ipsi_m <- as.matrix(in_out_scores_mean_ipsi_m[, -1]) # -1 to omit categories from matrix

clust <- hclust(dist(t(in_out_scores_mean_ipsi_m)))
my_cols_fun <- colorRampPalette(c("#3E1F4B","#8B008B","#F2F2F2"))

#plot as heatmap for all input to output neurons
figure_heatmap <- ggplot(as_tibble(in_out_scores_mean), aes(x = factor(type, level = output_order), y = input)) +
  geom_point(aes(col = ds_layer, size = conn_score), shape = 15) +
  theme_minimal() +
  theme(
    legend.position = 'right',
    text = element_text(color = 'grey40')
  ) +
  scale_size_area(max_size = 5) +
  scale_colour_gradientn(colours = my_cols_fun(desired_layers), limits=c(1, desired_layers)) +
  scale_y_discrete(limits = colnames(in_out_scores_mean_ipsi_m)[clust$order]) +
  guides(colour = guide_legend(override.aes = list(size=10)))  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size = 14, family = "Arial"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  facet_grid(. ~ factor(ipsi_contra, levels=c("ipsi","contra")), scales = "fixed")
ggsave("output/effective_connectivity_heatmap.pdf", figure_heatmap, dev=cairo_pdf, width=5000, height=length(input_types)*70+1000, units= "px")

#boxplots for the same data above
for (nth in 1:length(input_types)) {
  in_out_scores_type <- in_out_scores[in_out_scores$input == input_types[nth] & in_out_scores$type %in% output_order,]
  
  type_score <- ggplot(in_out_scores_type, aes(x=factor(type, level = output_order), y=conn_score)) +
    geom_boxplot(fill="#8B008B") + ylab("Score") + xlab("output neurons") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid = element_line(colour = "grey", linewidth = 0.2)) + ylim(0, ceiling(max(in_out_scores_type$conn_score))) +
    facet_grid(. ~ factor(ipsi_contra, levels=c("ipsi","contra")), scales = "fixed")
  
  type_layer <- ggplot(in_out_scores_type, aes(x=factor(type, level = output_order), y=ds_layer)) +
    geom_boxplot(fill="#8B008B") + ylab("Layer") + xlab("output neurons") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid = element_line(colour = "grey", linewidth = 0.2)) + ylim(1, desired_layers) +
    facet_grid(. ~ factor(ipsi_contra, levels=c("ipsi","contra")), scales = "fixed")
  
  legend <- get_legend(
    type_score +
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  title <- ggdraw()+
    draw_label(paste0(input_types[nth]," connectivity to leg MNs"), fontface='bold')
  
  figure_nwh <- plot_grid(type_score + theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position="none",plot.margin = margin(0.1,0.1,0,0.3, "cm")),
                          type_layer + theme(legend.position="none",strip.background = element_blank(),strip.text.x = element_blank(),plot.margin = margin(0.1,0.1,0,0.3, "cm")),
                          ncol = 1, common.legend = FALSE,
                          rel_heights = c(1.3, 0.8), align = "v", vjust=0, axis = "rlbt")
  figure2_nwh <- plot_grid(title, legend, nrow=1, rel_widths =c(1,1))
  figure3_nwh <- plot_grid(figure2_nwh, figure_nwh, ncol=1, rel_heights =c(0.1,1))
  ggsave(paste0("output/",input_types[nth], "_combined_plot.pdf"), figure3_nwh, dev=cairo_pdf, width=4715, height=3295, units= "px")
}


