#
# # code to simulate data?
#
#
# memory.limit(size = 250000)
# a <- data.frame(l = letters[1:26])
# b <- bind_rows(data.frame(l = letters[2:26]), data.frame(l = letters[1]))
# c <- bind_rows(data.frame(l = letters[3:26]), data.frame(l = letters[1:2]))
# d <- bind_rows(data.frame(l = letters[4:26]), data.frame(l = letters[1:3]))
# A1 <- data.frame(Let = LETTERS[1:26])
# A2 <- A1
# A3 <- A1
# A4 <- A1
# LA <- dplyr::bind_rows(A1, A2, A3, A4)
# l_abcd <- dplyr::bind_rows(a, b, c, d)
#
# Ll <- dplyr::bind_cols(LA, l_abcd) %>% dplyr::transmute(Ll = paste0(Let, l))
#
# # Using Gaussian Data
# #
# # Load predefined data
# n <- nrow(PC_data_SC)
# V <- colnames(data.frame(PC_data_SC)) # labels aka node names
# # estimate CPDAG
# # PC_skel <-  skeleton(suffStat = list(C = cor(PC_data), n = n),
# #                     indepTest = gaussCItest,labels = V, method ="stable", # indep.test: partial correlations
# #                     alpha=0.05, fixedGaps = NULL, fixedEdges = NULL,verbose = FALSE)
# # PC_DAG <-  pc(suffStat = list(C = cor(PC_data_SC), n = n),
# #              indepTest = gaussCItest,labels = V, skel.method ="stable", # indep.test: partial correlations
# #              alpha=0.05, fixedGaps = NULL, fixedEdges = NULL,verbose = FALSE,maj.rule = F, solve.confl = F)
#
# PC_ren <- c(Ll$Ll[1:length(colnames(PC_data_SC))])
# RENAME_PC <- dplyr::bind_cols(data.frame(colnames(PC_data_SC)), data.frame(PC_ren))
# colnames(RENAME_PC) <- c("Metabolite", "Exposure")
# colnames(PC_data_SC) <- PC_ren
# colnames(PC_data) <- PC_ren
#
# LABEL <- data.frame(Metabolite = RENAME_PC$Metabolite, Met_label = sapply(strsplit(as.character(RENAME_PC$Metabolite), split = "_", fixed = TRUE), function(x) (paste0(x[3], sep = "/", x[4]))))
# LABEL$Metabolite <- as.character(LABEL$Metabolite)
#
# #
# # Using Gaussian Data
# #
# # Load predefined data
# n <- nrow(PC_data_SC)
# V <- colnames(data.frame(PC_data_SC)) # labels aka node names
# PC_skel <- skeleton(
#   suffStat = list(C = cor(PC_data_SC), n = n),
#   indepTest = gaussCItest, labels = V, method = "stable", # indep.test: partial correlations
#   alpha = 0.05, fixedGaps = NULL, fixedEdges = NULL, verbose = FALSE
# )
#
# # Generate networks: DAG; Skeleton; Adjacency matrix
# n <- nrow(PC_data_SC)
# V <- colnames(data.frame(PC_data_SC)) # labels aka node names
#
# NW_PC_EPIC <- initEdgeAttribute(NW_PC_EPIC,
#   attribute.name = "pCor",
#   attribute.type = "numeric",
#   default.value = 1
# )
# NW_PC_EPIC <- initEdgeAttribute(NW_PC_EPIC,
#   attribute.name = "Est_range",
#   attribute.type = "char",
#   default.value = 1
# )
#
# for (i in 1:length(NW_PC_EPIC@edgeData@data))
# {
#   Y <- unlist(strsplit(attributes(NW_PC_EPIC@edgeData@data)[[1]][[i]], "[|]"))
#   print(Y)
#   edgeData(NW_PC_EPIC, from = Y[[1]], to = Y[[2]], "pCor") <- pCor_PC$estimate[Y[[1]], Y[[2]]]
# }
# for (i in 1:length(NW_PC_EPIC@edgeData@data))
# {
#   Y <- unlist(strsplit(attributes(NW_PC_EPIC@edgeData@data)[[1]][[i]], "[|]"))
#   print(Y)
#   edgeData(NW_PC_EPIC, from = Y[[1]], to = Y[[2]], "Est_range") <- as.character(format(round(pCor_PC$estimate[Y[[1]], Y[[2]]], 2), nsmall = 2))
# }
# Mets <- NW_PC_EPIC@nodes
