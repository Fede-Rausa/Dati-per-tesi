# Dati-per-tesi

Dati per la tesi in Fondamenti di Data Analytics

Il file marketCluster.R comprende parte del codice inerente la tesi sui settori Ateco in Italia.
Il file lawMining.R comprende il codice, quasi completo, per la tesi sul clustering delle leggi in Italia.

Il dataset impreseItaliane.csv contiene, per ogni impresa per cui sono disponibili questi valori:
RagioneSociale;Provincia;ATECO_2007;Ricavi;SalariStipendi;TotImposte;Costi_ricerca_pubb;Diritti_brevetto;Dipendenti;ROI;EBITDA;CostiMateriePrime;CostiServizi;
CostiGodimentoBeniDiTerzi;CostiImpiantoAmpl;ImmoMateriali;TurnoverRatio;IndLiquidit;IndiceIndipFinanz

Il dataset comprende oltre 173000 imprese (circa il 13% dell'economia in Italia). Alcune delle variabili continue sono valori espressi in migliaia di euro, altre sono tassi. La colonna Dipendenti comprende i conteggi dei dipendenti dell'azienda.
Questo dataset proviene dal database Aida, accessibile dal portale insuBRE per gli studenti.

Il file province-regioni include una semplice tabella con la corrispondente regione per ogni provincia. 
E' utile per filtrare i dati in base alla regione su tabelle conteneti solo le province.
latlon.csv contiene latitudine e longitudine di ogni capoluogo di provincia.
I file codici_ateco e macro_settori contengono i nomi corrispondenti ad ogni attività produttiva in base al codice ateco.

Il dataset codici040 contiene i testi di 40 codici di legge, suddivisi per articoli, con i relativi testi, per un totale di 18000 articoli
