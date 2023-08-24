# using the pipeline

I'm planning how to use which pipeline.

I should start from the data i need to use

# data

La run di Aionostat che ha generato i dati aveva 2 provette con dentro 2 fagi ciascuna (in realta' c'erano 3 provette ma la terza provetta non aveva fagi), Valentin ha sequenziato i seguenti campioni:
1. phage populations: 2 populations, 4 timepoints each. Campioni delle provette durante l'esperimento.
2. phage isolates: lest timepoint, 2 populations, 4 isolates each. Singoli fagi isolati dalle popolazioni dell'ultimo giorno, 4 fagi per popolazione.
Quello che vogliamo fare con i dati delle phage populations e' metterli nella pipeline che abbiamo usato per lo scorso esperimento, ma usando due references (i due fagi contenuti nella provetta) per ogni sample. Valentin vorrebbe usare la tua pipeline e voleva sapere come si possono mettere due references.
Per i phage isolates il piano e' assemblare ogni isolate e poi mapparci sopra le references dei fagi iniziali della provetta, in questo modo se c'e stata della ricombinazione dovremmo vederla facilmente. Per fare cio' posso usare la parte iniziale della mia pipeline che e' gia' predisposta per l'assembly e i mapping con le references.
Valentin ha anche sequenziato vari batteri dell'esperimento:
1. initial bacteria (wbbl+)
2. bacterial culture samples: samples dalle provette che contengono solo i batteri (3 samples), alla fine dell'esperimento. Ha raccolto questi samples perche' i batteri hanno iniziato a replicarsi piu' velocemente e vuole capire perche'.
3. resistant bacteria: 3 clones. Alcuni batteri sono sopravvissuti nella provetta con i fagi, ha isolato 3 colonie e le ha sequenziate.
Il piano che avevamo con questi dati e' assemblare il genoma iniziale dei batteri (1) e poi usarlo come reference per mapparci sopra i dati dei batteri delle culture (2). Per quanto riguarda i dati dei batteri resistenti, andrebbero assemblati e poi confrontati con i batteri iniziali.

In questo caso il dubbio che ho riguarda l'assembly dei genomi batterici, conviene usare trycycler (magari usando una vostra pipeline) o continuo con flye?
Poi pensavo di tenere i dati sui fagi e sui batteri completamente separati e fare due repositories, ma non so se e' una buona idea.

the idea is to keep bacteria and phage data separated and to set up two repositories.



