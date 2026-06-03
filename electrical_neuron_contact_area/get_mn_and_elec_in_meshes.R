library(malevnc)
#library(googlesheets4)
library(neuprintr)
library(Rvcg)


wing_haltere_mns = neuprint_list2df(neuprint_fetch_custom(cypher=paste0("MATCH (a:Neuron) WHERE a.class = 'motor neuron' AND a.subclass IN ['wm','hm'] RETURN a.bodyId AS bodyid, a.class AS class, a.group AS group, a.subclass AS subclass, a.type AS type, a.synonyms AS synonyms, a.rootSide AS root_side, a.somaSide AS soma_side, a.exitNerve AS exit_nerve, a.entry_nerve AS entry_nerve"), timeout=2000))

#read pre-curated spreadsheet of potential electrical interneurons that were ordered by low pre site count divided by volume
#change path to wherever its saved
putative_elec_in_ss = read.csv(file.path(getwd(), "sup_file5_putative_elec_in_by_pre_over_vol.csv"))

#create folders to save neuron mesh
manc_meshes_path = file.path(getwd(), "manc_meshes")
if (!dir.exists(manc_meshes_path)) dir.create(manc_meshes_path)
manc_meshes_subpath = file.path(getwd(), "manc_meshes", "mn")
if (!dir.exists(manc_meshes_subpath)) dir.create(manc_meshes_subpath)

for (i in wing_haltere_mns$bodyid) {
	obj = read_manc_meshes(i)
	#scale to microns
	obj[[1]]$vb[1:3,] = obj[[1]]$vb[1:3,]/1000
	vcgObjWrite(obj[[1]], filename = paste0(manc_meshes_subpath, '/', i))
}


manc_meshes_subpath = file.path(getwd(), "manc_meshes", "elec_in")
if (!dir.exists(manc_meshes_subpath)) dir.create(manc_meshes_subpath)

#get only top putative electrical IN in upper tectulum by pre site count over volume, that have no DCVs and no tracing issues
putative_elec_in_ss$'manually_inspected'[is.na(putative_elec_in_ss$'manually_inspected')] = FALSE
putative_elec_in_ss$dont_consider_electrical[is.na(putative_elec_in_ss$dont_consider_electrical)] = FALSE
putative_elec_in_ss$in_utct[is.na(putative_elec_in_ss$in_utct)]=FALSE
putative_elec_in_ss = putative_elec_in_ss[putative_elec_in_ss$manually_inspected, ]
putative_elec_in_ss = putative_elec_in_ss[!(putative_elec_in_ss$has_dcvs | putative_elec_in_ss$dont_consider_electrical) & putative_elec_in_ss$in_utct,]
#putative_elec_in_ss = putative_elec_in_ss[1:30,]
temp_elec_in_bodyids = neuprint_list2df(neuprint_fetch_custom(cypher=paste0("MATCH (a:Neuron) WHERE a.group IN [", paste0(putative_elec_in_ss$group, collapse = ","), "] RETURN a.bodyId AS bodyid"), timeout=2000))


for (i in temp_elec_in_bodyids$bodyid) {
	obj = read_manc_meshes(i)
	#scale to microns
	obj[[1]]$vb[1:3,] = obj[[1]]$vb[1:3,]/1000
	vcgObjWrite(obj[[1]], filename = paste0(manc_meshes_subpath, '/', i))
}

manc_meshes_subpath = file.path(getwd(), "manc_meshes", "mn_in_intersection")
if (!dir.exists(manc_meshes_subpath)) dir.create(manc_meshes_subpath)