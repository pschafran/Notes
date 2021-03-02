library(UpSetR)
upsetMatrix <- read.csv(file.choose())
View(upsetMatrix)
upsetMatrix2 <- upsetMatrix[ -c(5, 11, 22, 23, 26)]
View(upsetMatrix2)
upset(upsetMatrix, nsets = 100, order.by = "degree")
upset(upsetMatrix, nsets = 100, nintersects = NA, order.by = "freq", 
  queries = list(
  list(query = intersects, params = list("vali", "engl.N"), color = "red", active = T), 
  list(query = intersects, params = list("vali", "engl.S"), color = "orange", active = T), 
  list(query = intersects, params = list("echi", "engl.N"), color = "green", active = T), 
  list(query = intersects, params = list("Uwhr", "mlpd.3", "silv"), color = "purple", active = T),
  list(query = intersects, params = list("Uwhr", "viri"), color = "blue", active = T), 
  list(query = intersects, params = list("mlpd.3", "silv"), color = "light blue", active = T), 
  list(query = intersects, params = list("mlpd.1", "mlpd.2", "snow"), color = "brown", active = T),
  list(query = intersects, params = list("engl.S", "matt", "vali"), color = "magenta", active = T),
  list(query = intersects, params = list("bola", "engl.S", "matt"), color = "gold", active = T),
  list(query = intersects, params = list("bola", "Uwhr", "silv", "viri"), color = "light green", active = T),
  list(query = intersects, params = list("bola", "silv"), color = "pink", active = T)
   ),
 # sets = c("bola","chap","echi","engl.N","engl.S","flac","Lear","matt","miss","mlpd.1","mlpd.2","mlpd.3","mlpd.4","mlsp","prot","silv","snow","Uwhr","vali","viri"),
  sets = c("viri","vali","Uwhr","snow","silv","prot","mlsp","mlpd.4","mlpd.3","mlpd.2","mlpd.1","miss","matt","Lear","flac","engl.S","engl.N","echi","chap","bola"),
  keep.order = TRUE,
  sets.x.label = "Num. of Diploid OTUs",
  mainbar.y.label = "Num. of Samples",
  text.scale = c(2, 1.7, 2, 1.7, 2, 2),
  point.size = 3.5,
  line.size = 1
  )


