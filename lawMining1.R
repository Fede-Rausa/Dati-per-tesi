
###################ATTI
##importa librerie:
library(quanteda)   ##text mining
library(FactoMineR) ##grafici clustering
library(factoextra) ##grafici clustering
library(cluster)    ##silhouette, pam/kmedoids
library(fpc)        ##Calinski-Harabasz, dbscan
library(clusterSim) ##Davis-Bouldin
library(DescTools)  ##Gini
library(dendextend) ##get_branches_heights
library(tidyverse)           ##str_length
library(quanteda.textplots)  ##textplot_wordcloud
library(qgraph)     ##qgraph


##DEFINISCI FUNZIONI



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


######DTMatrix(): costruisce la document term matrix di un vettore di testi. 
######Richiede anche il vettore dei relativi titoli. Fa anche la wordcloud se richiesta

DTMatrix = function(testi, nomi, wcloud=F, min_freq=10){
dtm <- testi |>
  tokens(remove_punct = TRUE, remove_numbers = TRUE, remove_symbols = TRUE) |>
  tokens_tolower() |>
  tokens_remove(stopwords('it')) |>
  tokens_wordstem() |>
  dfm()

# Riduci la DTM alle parole con una frequenza minima di 10
dtm <- dfm_trim(dtm, min_termfreq = min_freq)

if (wcloud){
textplot_wordcloud(dtm, max_words = 50,color = c('blue', 'orange', 'red'))
}

dtm = convert(dtm, to = "data.frame")

idBadWords = unique(c(1,which(str_length(colnames(dtm))<2)))   #,which(grepl('.',colnames(dtm)))))
dtm = dtm[,-idBadWords]

rownames(dtm) = nomi
dtm
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

plot(Koptions, sdev, main='Elbow method to find optimal k (PCA)', sub='where k is the numer of components',
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

elbowH = function(dati, maxk=nrow(dati), outK=F, metrica='euclidean', metodo='complete', cex=1){
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
plot(hc, cex = cex)
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

























##############ANALISI TESTI DELLE NORME IN ITALIA


##importa dati:

path = 'https://raw.githubusercontent.com/Fede-Rausa/Dati-per-tesi/main/codici040.csv'
dati = read.csv(path, sep=';')  ##dati per articolo


#######dataset dei codici e dei testi unici
atti = unique(dati[,'nomeAtto'])
lenStr = str_length(atti)
atti = atti[-which(lenStr<5)]

##ricostruzione atti: questo codice richiede qualche minuto di calcolo 
testi_codici = c()
for (atto in atti){
id = which(dati[,'nomeAtto']==atto)
articoli = dati[id,'testoArt']
testo=''
for (t in articoli){
testo=paste(testo,t)
}
testi_codici = c(testi_codici, testo)
}

##creazione DTM: questo codice richiede qualche minuto di calcolo 
dtm_atti = DTMatrix(testi_codici, atti, wcloud=F, min_freq=20)
dim(dtm_atti)

##prima visualizzazione
textplot_wordcloud(as.dfm(dtm_atti), max_words=300)

###conteggio articoli per ciascun codice
table(dati[,'nomeAtto'])

#vedo che la maggior parte dei testi contiene più o meno 200 articoli
#infatti il cluster che conterrà più testi sarà determinato 
#avrà in media quel numero di articoli


#######Principal Component Analysis: estraggo solo le componenti principali nel fare clustering
plot(prcomp(dtm_atti), type='l')

pca = prcomp(dtm_atti, scale=T)
dtm_pca = pca$x[,1:8]


#######INTERPRETAZIONI RISULTATI PCA
PCA(dtm_atti)
summary(pca)
elbowPCA(pca)



####DISTANCE MATRIX DEI CODICI


heatmapDM(dtm_pca, metrica = 'manhattan')
heatmapDM(dtm_pca, metrica = 'euclidean')
heatmapDM(dtm_pca, metrica = 'maximum')

grafoDM(dtm_pca, metrica = 'manhattan')
grafoDM(dtm_pca, metrica = 'euclidean')
grafoDM(dtm_pca, metrica = 'maximum')


######ANALISI DI CLUSTERING DEI CODICI

bestCluster(dtm_pca, maxk=nrow(dtm_pca))


###kmeans
elbowK(dtm_pca, maxk=nrow(dtm_pca))

kmean = kmeans(scale(dtm_pca), 10)
fviz_cluster(list(data = scale(dtm_pca), cluster = kmean$cluster),
 stand = FALSE, geom = "text", labelsize=7)

sil = silhouette(kmean$cluster,distmat)
plot(sil)



###hierarchical
elbowH(dtm_pca, metodo='complete', metrica='maximum')

distmat = dist(scale(dtm_pca),method='maximum')
albero = hclust(distmat, method='single')
plot(albero, cex=0.5)
alberoClust = cutree(albero, k=10)
fviz_cluster(list(data = scale(dtm_pca), cluster = alberoClust),
stand = FALSE, geom=c("point"))

sil = silhouette(alberoClust,distmat)
plot(sil)


###kmedoids
distmat = dist(scale(dtm_pca))
kmedoid <- pam(distmat, 10)
fviz_cluster(list(data = scale(dtm_pca), cluster = kmedoid$clustering),
 stand = FALSE, geom = "point")

sil = silhouette(kmedoid$clustering,distmat)
plot(sil)


###dbscan
distmat = dist(scale(dtm_pca))
dbscan = dbscan(scale(dtm_pca), eps=5, MinPts = 14, method='hybrid')

fviz_cluster(list(data = scale(dtm_pca), cluster = dbscan$cluster),
stand = FALSE, geom = "text", labelsize=8, main ='dbscan',
axes = c(1, 2))

sil = silhouette(dbscan$cluster,distmat)
plot(sil)




