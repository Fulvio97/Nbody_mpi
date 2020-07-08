# Nbody_mpi
Progetto di PCPC A. A. 19/20
## Obiettivo
Risolvere il problema Nbody attraverso la computazione parallela. Il progetto sfrutta MPI per la comunicazione fra processi.
Sia il numero di iterazioni che i corpi sono scelti dall'utente tramite riga di comando.
## Implementazione
Il file Nbody.c può essere compilato ed eseguito grazie ai comandi **mpicc** e **mpirun**, nell'ultimo caso da riga di comando bisogna aggiungere sia il numero di iterazioni che il numero di corpi da genereare
## Test correttezza
Per verificare la correttezza dell'output si possono sfruttare SequentialTest e ParallelTest, per fare ciò basta lanciare prima SequentialTest per generare tutti i valori e l'output di una esecuzione sequenziale e poi far partire ParallelTest il qule legge i corpi dal file generato da ParallelTest ed esegue la computazione parallela; i file da controllare sono log.txt e TestLog.txt
