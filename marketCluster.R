
##########################################IMPORTAZIONE LIBRERIE
library(FactoMineR) ##grafici clustering
library(factoextra) ##grafici clustering
library(cluster)    ##silhouette, pam/medoids
library(fpc)        ##Calinski-Harabasz, dbscan
library(clusterSim) ##Davis-Bouldin
library(DescTools)  ##Gini
library(dendextend) ##get_branches_heights
library(qgraph)     ##qgraph




############################################DEFINIZIONE FUNZIONI


####decimaliCSV(): dato un vettore, sostituisce le virgole con i punti e lo trasforma in una variabile continua

decimaliCSV = function(vettore){
for (i in 1:length(vettore)){
vettore[i] = sub(',','.',vettore[i])
}
vettore=as.double(vettore)
vettore
}



####filtra(): un dataset in base a una condizione specificata:
filtra = function(dataset, opzione='=', nomeVariabile, termine){
dati = dataset

if (opzione == '='){
id = which(dati[[nomeVariabile]]==termine)
dati = dati[id,]
}

if (opzione == '>'){
id = which(dati[[nomeVariabile]]>termine)
dati = dati[id,]
}

if (opzione == '<'){
id = which(dati[[nomeVariabile]]<termine)
dati = dati[id,]
}

dati
}



###memoria(): restituisce il peso di calcolo di un oggetto
memoria= function(dati){
nomi = c('Bytes','kiloBytes','megaBytes','gigaBytes')
bit = object.size(dati)
pesi = c(bit, bit/1000, bit/1000000, bit/1000000000)
mat = cbind(pesi, nomi)
mat
}

####pecoreCapri(): dato un dataset, restituisce la matrice delle variabili continue e quella delle variabili categoriche
###per raccogliere separatamente i due output 
###bisogna fare VarCat = pecoreCapri(mat)[[1]] e VarNum = pecoreCapri(mat)[[2]]

pecoreCapri = function(mat){

VarCat = rep(0, nrow(mat))
VarNum = rep(0, nrow(mat))
nomiCat = c()
nomiNum = c()

for (i in 1:ncol(mat)){
if(typeof(mat[,i])=='double'){
VarNum = cbind(VarNum, mat[[i]])
nomiNum = c(nomiNum, colnames(mat)[i])
}else{
VarCat = cbind(VarCat, mat[,i])
nomiCat = c(nomiCat, colnames(mat)[i])
}}
VarCat = VarCat[,-1]
VarNum = VarNum[,-1]

colnames(VarCat) = nomiCat
colnames(VarNum) = nomiNum

list(VarCat, VarNum)
}



#####idOutZ(): trova gli indici degli outliers in un vettore (secondo il metodo Zscore

idOutZ = function(vettore){
  vet = vettore 
  vetNorm = as.numeric(scale(vet))     #normalizza  il vettore
  idMaleInf = which(vetNorm<(-3))       #trova gli indici degli outliers, nel vettore normalizzato
  idMaleSup = which(vetNorm>3) 
  idMale = c(idMaleSup, idMaleInf)
  idMale
}



#####outCluser(): rimuove gli elementi appartenenti al cluster più piccolo da un dataset
#####restituisce il dataset filtrato

outCluster = function(dati, vettoreCluster){
conto = table(vettoreCluster)
categorie = data.frame(conto)[,1]
outClust = categorie[which(conto==min(conto))]
outClust = which(vettoreCluster %in% outClust)
dati = dati[-outClust,]
dati
}

#####getmode(): resituisce il nome del cluster più grande, dato il vettore degli indici
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}



#####getstrange(): resituisce il nome del cluster più piccolo, dato il vettore degli indici
getstrange <- function(v) {
   uniqv <- unique(v)
   uniqv[which.min(tabulate(match(v, uniqv)))]
}



####getCluster(): restituisce gli elementi appartenenti al cluster specificato
getCluster = function(dati, vettoreCluster, clu=1, nomi=T){
id = which(vettoreCluster==clu)
if (nomi){
cluster = rownames(dati)[id]
}else{
cluster = dati[id,]
}
cluster
}


#####elbowPCA(): crea il grafico elbow e trova un numero di componenti principali ottimo

elbowPCA = function(pca, outK=F){
sdev = pca$sdev
Koptions = 1:length(sdev)
tab = cbind(Koptions,sdev)

distanze = c()
a = tab[1,]
b = tab[nrow(tab),]
m = -(abs(b[2]-a[2]))/abs((b[1]-a[1])) 
q = a[2] - m * a[1]
mp = -1/m

for (r in 1:nrow(tab)){
c = tab[r,]
d = c(0,0)
qp = c[2] - mp * c[1]
d[1] = (qp - q) / (m - mp)
d[2] = m * d[1] + q
distanze = c(distanze, dist(rbind(c,d)))
}
maxc = tab[which(distanze==max(distanze)),]
qp = maxc[2] - mp * maxc[1]
maxd = c(0,0)
maxd[1] = (qp - q) / (m - mp)
maxd[2] = m * maxd[1] + q

bestK = maxc[1]
if (outK){
return(bestK)
}else{
propvar = sum(sdev[1:bestK]^2)/sum(sdev^2)

plot(Koptions, sdev, main='Elbow method to find optimal k (PCA)',
 sub='where k is the number of components',
 type='o', col='blue')
lines(rbind(a,b), lty=2, col='red')
lines(rbind(maxc,maxd), lty=2, col='red')
lines(maxc[1],maxc[2], type='p', pch=8, col='red')
text(quantile(Koptions, 0.75),max(sdev), 
labels=paste('good k =',as.character(bestK),'\n', 'cum var =', 
as.character(round(propvar*100, digits=2)), '%'))
}}




#####elbowK(): crea il grafico elbow per il kmeans clustering, se maxk non è un numero è il massimo possibile
#####calcola anche un valore ottimo di k, corrispondente all'angolo del gomito

elbowK = function(dati, maxk=15, outK=F){
maxk = ifelse(is.double(maxk), maxk, nrow(dati)-1)
Koptions=1:maxk
wss = c()
for (k in Koptions){
model = kmeans(dati, k)
wss= c(wss, model$tot.withinss)
}

tab = cbind(Koptions,wss)
distanze = c()
a = tab[1,]
b = tab[nrow(tab),]
m = -(abs(b[2]-a[2]))/abs((b[1]-a[1])) 
q = a[2] - m * a[1]
mp = -1/m

for (r in 1:nrow(tab)){
c = tab[r,]
d = c(0,0)
qp = c[2] - mp * c[1]
d[1] = (qp - q) / (m - mp)
d[2] = m * d[1] + q
distanze = c(distanze, dist(rbind(c,d)))
}
maxc = tab[which(distanze==max(distanze)),]
qp = maxc[2] - mp * maxc[1]
maxd = c(0,0)
maxd[1] = (qp - q) / (m - mp)
maxd[2] = m * maxd[1] + q

bestK = maxc[1]

if (outK){
return(bestK)
}else{
plot(Koptions, wss, main='Elbow method to find optimal k (kmeans)', type='o', col='blue')
lines(rbind(a,b), lty=2, col='red')
lines(rbind(maxc,maxd), lty=2, col='red')
lines(maxc[1],maxc[2], type='p', pch=8, col='red')
text(quantile(Koptions, 0.75),max(wss), labels=paste('good k =',as.character(bestK)))
}}



#####elbowH(): crea il grafico elbow per il clustering gerarchico, se maxk non è un numero è il massimo possibile
#####calcola anche un valore ottimo di k, corrispondente all'angolo del gomito
#####dati è un dataset di variabili continue, 
#####outK è un valore bool che dice se fare il grafico o restituire k
#####metrica e metodo sono i criteri di costruzione del modello gerarchico

elbowH = function(dati, maxk=nrow(dati), outK=F, metrica='euclidean', metodo='complete'){
#maxk = nrow(dati)
Koptions=2:maxk

dist_matrix <- dist(dati, method=metrica)
hc <- hclust(dist_matrix, method=metodo)
heights = rev(get_branches_heights(as.dendrogram(hc)))

tab = cbind(Koptions,heights)
distanze = c()
a = tab[1,]
b = tab[nrow(tab),]
m = -(abs(b[2]-a[2]))/abs((b[1]-a[1])) 
q = a[2] - m * a[1]
mp = -1/m

for (r in 1:nrow(tab)){
c = tab[r,]
d = c(0,0)
qp = c[2] - mp * c[1]
d[1] = (qp - q) / (m - mp)
d[2] = m * d[1] + q
distanze = c(distanze, dist(rbind(c,d)))
}
maxc = tab[which(distanze==max(distanze)),]
qp = maxc[2] - mp * maxc[1]
maxd = c(0,0)
maxd[1] = (qp - q) / (m - mp)
maxd[2] = m * maxd[1] + q

bestK = maxc[1]
bestH = maxc[2]

if (outK){
return(bestK)
}else{
par(mfrow=c(1,2))
plot(Koptions, heights, main='Elbow method to find optimal k (hclust)', type='o', col='blue')
lines(rbind(a,b), lty=2, col='red')
lines(rbind(maxc,maxd), lty=2, col='red')
lines(maxc[1],maxc[2], type='p', pch=8, col='red')
text(quantile(Koptions, 0.75),max(heights), labels=paste('good k =',as.character(bestK)))
plot(hc)
lines(c(0,10000),c(bestH, bestH), lty=2, col='red')
}}






########bestCluster: produce una tabella comparativa dei metodi di clustering principali

bestCluster = function(dati, maxk=nrow(dati)){
hmetriche = c('euclidean', 'manhattan', 'maximum')
hmetodi = c('complete', 'single', 'average')

tab = rep(0,5)
modelli=c('kmeans')

maxg = maxk   #ifelse((nrow(dati)<100),T,100)

k = elbowK(dati, maxk=maxg, outK=T)
Kmean = kmeans(dati, k)
sil = silhouette(Kmean$cluster, dist(dati, method='euclidean'))
msil = mean(sil[,3])
sizeDev = Gini(table(Kmean$cluster))
davbou = index.DB(dati,Kmean$cluster)$DB
calhar = cluster.stats(dati,Kmean$cluster)$ch

tab = rbind(tab, c(k, msil,sizeDev, calhar, davbou))

for (d in hmetriche){
for (m in hmetodi){
nome = paste('hierarchic-',d,'-',m, sep='')
modelli = c(modelli, nome)
distmat = dist(dati, method=d)
model = hclust(distmat, method=m)
k = elbowH(dati, outK=T, metrica=d, metodo=m)
model = cutree(model, k=k)
sil = silhouette(model, distmat)
msil = mean(sil[,3])
sizeDev = Gini(table(model))
davbou = index.DB(dati,model)$DB
calhar = cluster.stats(dati, model)$ch

tab = rbind(tab, c(k, msil, sizeDev, calhar, davbou))
}}

tab = data.frame(tab[-1,])
rownames(tab) = modelli

for (c in 1:ncol(tab)){
for (r in 1:nrow(tab)){
tab[r,c] = as.numeric(tab[r,c])
}}

colnames(tab) = c('k', 'mean silhouette','Gini-ClusterSize', 'Calinski-Harabasz', 'Davis-Bouldin')
tab
}


###heatmapDM(): costruisce la heatmap della distance matrix
heatmapDM = function(dati, titolo = '', ordina=T , metrica='euclidean', metodo='complete',
colori = colorRampPalette(c('white','yellow','orange','red','darkred'))(100)){

distmat = dist(scale(dati),method=metrica)

if (ordina){
hc <- hclust(distmat, method=metodo)
distmat = as.matrix(distmat)
matrice_distanze_ordinate <- (distmat)[hc$order, hc$order]
nomi_righe <- rownames(matrice_distanze_ordinate)
dim <- nrow(distmat)
image(1:dim, 1:dim, matrice_distanze_ordinate, 
axes = F, xlab="", ylab="", main=titolo, col=colori)
axis(1, 1:dim, nomi_righe[hc$order], cex.axis = 0.6, las=3)
axis(2, 1:dim, nomi_righe[hc$order], cex.axis = 0.6, las=1)

}else{
matrice_distanze_ordinate = as.matrix(distmat)
nomi_righe <- rownames(matrice_distanze_ordinate)
dim <- nrow(distmat)
image(1:dim, 1:dim, matrice_distanze_ordinate, 
axes = F, xlab="", ylab="", main=titolo, col=colori)
axis(1, 1:dim, nomi_righe, cex.axis = 0.4, las=3)
axis(2, 1:dim, nomi_righe, cex.axis = 0.4, las=1)
}}


###grafoDM(): costruisce il grafo della distance matrix

grafoDM = function(dati, titolo='', colore=c('green', 'darkgreen')
,size=5, metrica='manhattan', cex=3, str=12){

nomi = substr(rownames(dati),1,str)
distmat = dist(scale(dati), method=metrica)
dist_mi <- 1/as.matrix(distmat)
qgraph(dist_mi, layout='spring', vsize=size, posCol=colore, 
borders=T, legend=F, title=titolo, label.cex=cex, labels = nomi)
}



####freeClu(): consente di produrre rapidamente i grafici di fviz_cluster, per ogni metodo possibile

freeClu = function(dati ,dist='euclidean', mod='kmean', k=5, link='complete', eps=5, MinPts=14,
 geom='point', size = 8, main = 'freeClu', axes=c(1,2)){

dati = scale(dati)

if (mod=='kmean'){
clu = kmeans(dati,k)
cluster = clu$cluster
}

if (mod=='kmedoid'){
distmat = dist(dati, method=dist)
clu = pam(distmat,k)
cluster = clu$clustering
}

if (mod=='dbscan'){
clu = dbscan(dati, eps=eps, MinPts = MinPts, method='hybrid')
cluster = clu$cluster
}

if (mod=='tree'){
distmat = dist(dati, method=dist)
clu = hclust(distmat, method=link)
clu = cutree(clu, k=k)
cluster = clu
}

fviz_cluster(list(data = dati, cluster = cluster),
stand = FALSE, geom = geom, labelsize=size, main = main,
axes = axes)
}

















###############################ANALISI DELLE IMPRESE ITALIANE


##importazione dati
link = 'https://raw.githubusercontent.com/Fede-Rausa/Dati-per-tesi/main/impreseItaliane.csv'
dati = read.csv(link, sep=';')
dim(dati)
colnames(dati)

#########data cleaning
summary(dati)

correggi = c('ROI', 'TurnoverRatio','IndLiquidit','IndiceIndipFinanz')
for (nome in correggi){
dati[,nome] = decimaliCSV(dati[,nome])
}

summary(dati)

################################ANALISI DATI

##########ANALISI FATTORIALE
dati_numerici = dati[,4:ncol(dati)]
dati_scale = scale(dati_numerici)
pca = prcomp(dati_numerici, scale=T)

summary(pca)

#elbowPCA(pca)
#PCA(dati_scale)



#######RIMOZIONE OUTLIER

out1=idOutZ(pca$x[,1])
length(out1)
out2=idOutZ(pca$x[,2])
length(out2)

outliers = c(out1, out2)
length(outliers)
dati_clean = dati[-outliers,]





###########ANALISI DI BASE


##conteggio settori
sort(table(dati_clean[,'ATECO_2007']))  
###il settore più numeroso è 
###412000 COSTRUZIONE DI EDIFICI RESIDENZIALI E NON RESIDENZIALI

##conteggio macro settori
sort(table(substr(dati_clean[,'ATECO_2007'],1,2)))
###il macrosettore più numeroso è
###46 COMMERCIO ALL'INGROSSO (ESCLUSO QUELLO DI AUTOVEICOLI E DI MOTOCICLI)

##c'è correlazione tra settori e provincie?
tab = table(dati_clean[,c('ATECO_2007', 'Provincia')])
chisq.test(tab)





##########RIDUZIONI DATASET
##questo codice richiede molto tempo, se fatto correre su 173000 righe
##In alcuni casi sembra impossibile fare calcoli sull'intero dataset
##come ad esempio per la distance matrix 
##perciò limito l'analisi solo ad alcuni settori ateco
##e ai valori medi di ciascun settore/macro_settore



######costruzione dataset SETTORI (1500 righe)

ateco = unique(dati_clean[,'ATECO_2007'])

continue = colnames(dati)[4:ncol(dati_clean)]

dati_cleanII = dati_clean[,continue]
dati_settori = 1:ncol(dati_cleanII)

for (code in ateco){
dati_code = dati_cleanII[which(dati_clean[,'ATECO_2007'] == code),]
row_code = c()

for (c in 1:ncol(dati_code)){
med = mean(dati_code[,c])
row_code = c(row_code, med)
}
dati_settori = rbind(dati_settori, row_code)
}
dati_settori = dati_settori[-1,]
rownames(dati_settori) = ateco
colnames(dati_settori) = continue



######costruzione dataset MACRO_SETTORI (87 righe)

macroVet = substr(dati_clean[,'ATECO_2007'],1,2)
macro = unique(macroVet)

dati_cleanII = dati_clean[,continue]
dati_macro = 1:ncol(dati_cleanII)

for (code in macro){
dati_code = dati_cleanII[which(macroVet == code),]
row_code = c()

for (c in 1:ncol(dati_code)){
med = mean(dati_code[,c])
row_code = c(row_code, med)
}
dati_macro = rbind(dati_macro, row_code)
}
dati_macro = dati_macro[-1,]
rownames(dati_macro) = macro
colnames(dati_macro) = continue


######analisi univariata


#boxplot(dati_clean[,continue], las=2)
#boxplot(dati_settori, las=2)
boxplot(dati_macro, las=2)




############confronto del costo computazionale
memoria(dati)  
memoria(dati_settori)
memoria(dati_macro)   ##utilizzabile per la distance matrix



##########IMPORTAZIONE ETICHETTE
###importo una tabella contenente i nomi dei macrosettori
link = 'https://raw.githubusercontent.com/Fede-Rausa/Dati-per-tesi/main/macro_settori.csv'
etichette = read.csv(link, sep=';')

idMacro = sort(as.double(sub(',', '',rownames(dati_macro))))
for (i in 1:length(idMacro)){
if (idMacro[i] %in% etichette[,"macroID"]){
id = which(etichette[,"macroID"]==idMacro[i])
idMacro[i] = etichette[id,"macroNOME"]
}}

##si possono utilizzare i nomi dei macro settori 
##al posto delle prime 2 cifre ateco
##al fine dei grafici delle distance matrix però non conviene farlo
##mentre può essere utile per i grafici dei cluster
nomiMacro = idMacro
idMacro = rownames(dati_macro)




##############DISTANCE MATRIX DEI MACROSETTORI

colori = colorRampPalette(c('blue','yellow','orange','red','darkred'))(100)
title='somiglianze macrosettori ateco 2007'
metodi = c('complete', 'average', 'single')
metriche = c('manhattan','euclidean','maximum')

rownames(dati_macro)= idMacro

par(mfrow=c(3,3))
for (d in metriche){
for (m in metodi){
heatmapDM(dati_macro, metodo=m, metrica=d, titolo=paste(d,'-',m), colori= colori)
}}




##metriche adottate
#par(mfrow=c(2,2))
heatmapDM(dati_macro, ordina=F, titolo=title, metrica='manhattan')

heatmapDM(dati_macro, metodo='complete', titolo=title, metrica='manhattan')

grafoDM(dati_macro, titolo=title, colore=c('red', 'darkred'), size=4, cex=2, metrica='manhattan')




##########CLUSTERING DEI MACROSETTORI

bestCluster(scale(dati_macro))
elbowK(scale(dati_macro), maxk=nrow(dati_macro))
elbowH(scale(dati_macro))

##kmeans
set.seed(123)
rownames(dati_macro) = nomiMacro  #imposto etichette
kmean = kmeans(scale(dati_macro), 4)
fviz_cluster(list(data = scale(dati_macro), cluster = kmean$cluster), 
stand = FALSE, geom=c("text"), labelsize=6, main=title)
table(kmean$cluster)


rownames(dati_macro) = idMacro
distmat = dist(scale(dati_macro),method='manhattan')
albero = hclust(distmat, method='average')
plot(albero)
alberoClust = cutree(albero, k=21)
table(alberoClust)
fviz_cluster(list(data = dati_macro, cluster = alberoClust), 
stand = FALSE, geom=c("text"), labelsize=8, main = 'macro_settori_ateco')





##########CLUSTERING DEI SETTORI
title='cluster settori ateco 2007'

elbowK(scale(dati_settori), maxk=500)

kmean = kmeans(scale(dati_settori), 10)
fviz_cluster(list(data = scale(dati_settori), cluster = kmean$cluster),
 stand = FALSE, geom = "point")


distmat = dist(scale(dati_settori))
kmedoid <- pam(distmat, 10)
fviz_cluster(list(data = scale(dati_settori), cluster = kmedoid$clustering),
 stand = FALSE, geom = "point")


dbscan = dbscan(scale(dati_settori), eps=5, MinPts = 14, method='hybrid')
table(dbscan$cluster)
fviz_cluster(list(data = scale(dati_settori), cluster = dbscan$cluster),
stand = FALSE, geom = "point")


elbowH(scale(dati_settori))
distmat = dist(scale(dati_settori),method='manhattan')
albero = hclust(distmat, method='complete')
plot(albero, cex=0.2, labels=F)
alberoClust = cutree(albero, k=88)
table(alberoClust)
fviz_cluster(list(data = scale(dati_settori), cluster = alberoClust),
stand = FALSE, geom=c("point"))






##########CLUSTERING DI UN SETTORE SPECIFICO

title = 'Produzione di aeromobili in Italia'
settore = 303009

dati_aerei = filtra(dati_clean,'=','ATECO_2007',settore)
rownames(dati_aerei) = dati_aerei[,"RagioneSociale"]
provincie = dati_aerei[,"Provincia"]
dati_aerei = dati_aerei[,continue]


heatmapDM(dati_aerei, metodo='complete',metrica='manhattan', titolo=title)
grafoDM(dati_aerei, metrica='manhattan')

########CLUSTERING

###kmeans
elbowK(dati_aerei, maxk=nrow(dati_aerei))

kmean = kmeans(scale(dati_aerei), 10)
fviz_cluster(list(data = scale(dati_aerei), cluster = kmean$cluster),
 stand = FALSE, geom = "text", labelsize=7)

silkmean = silhouette(kmean$cluster,distmat)
plot(silkmean)




###hierarchical
elbowH(dati_aerei, metodo='complete', metrica='maximum')

distmat = dist(scale(dati_aerei),method='maximum')
albero = hclust(distmat, method='single')
plot(albero, cex=0.5, main=title)
alberoClust = cutree(albero, k=10)
fviz_cluster(list(data = scale(dati_aerei), cluster = alberoClust),
stand = FALSE, geom=c("point"))

sil = silhouette(alberoClust,distmat)
plot(sil)


###kmedoids
distmat = dist(scale(dati_aerei))
kmedoid <- pam(distmat, 10)
fviz_cluster(list(data = scale(dati_aerei), cluster = kmedoid$clustering),
 stand = FALSE, geom = "point")

sil = silhouette(kmedoid$clustering,distmat)
plot(silkmean)


###dbscan
dbscan = dbscan(scale(dati_aerei), eps=5, MinPts = 14, method='hybrid')
fviz_cluster(list(data = scale(dati_aerei), cluster = dbscan$cluster),
stand = FALSE, geom = "point", labelsize=8, main ='dbscan')

sil = silhouette(dbscan$cluster,distmat)
plot(sil)




