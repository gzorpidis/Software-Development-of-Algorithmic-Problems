# Ανάπτυξη Λογισμικού για δυσεπίλυτα Αλγοριθμικά Προβλήματα
## Εργασία 2



---


##### Εκπόνηση από: 

1) Ζορπίδης Γεώργιος(1115202000055)

2) Κωνσταντίνος Φράγκος (1115202000207)

---
##### Λίγα πράγματα για την υλοποίηση

Οι αλγόριθμοι για τους nearest neighbours μέσω γράφων (με GNN, MRNG)  δημιουργήθηκαν για χρήση με εικόνες διάστασης 10, και έχουν το πρότυπο του MNIST dataset, για αυτό και γίνεται εκτενής χρήση του τύπου `uint8_t`.

Επίσης, υπάρχει χρήση του flag `-g <file_name> ` για κατασκευή των γράφων προς χρήση για το GNN και MRNG search, μέσω του διαβάσματος απο το αρχείο αυτό. Αν θέλουμε η κατασκεύη να γίνει απτην αρχή, τότε δεν πρέπει να υπάρξει εισαγωγή της παραμέτρου `-g`. Αν το `-g` flag δεν έχει δωθεί, βγάζουμε το αποτέλεσμα ανάλογα με την μέθοδο και σε αρχείο. Ωστόσο, επειδή το MRNG χρησιμοποιεί threads για αποδοτικότητα, ένα αρχείο εξάγεται για κάθε thread.Στο τέλος τα αρχεία μπορούν να ενωθούν σε 1 εννοιαίο για χρήση μετά, αλλά πρέπει να γίνει από τον χρήστη χειροκίνητα, δηλαδή δεν υλοποιείται στον κώδικά μας αυτή η λειτουργία.




Για την διαφοροποίηση των γράφων από gnn και mrng, έχουμε ορίσει: Ο γράφος που διαβάζεται από τον GNN να περιέχει εγγραφές της μορφής:  ` <node_id> : <neighbor1> , <neigbor2> , ...` (comma seperated οι γείτονες), ενώ ο γράφος για τον MRNG να περιέχει εγγραφές της μορφής: `<node_id> : <neighbor1> <neighbor2>` δηλαδή να είναι space seperated.



#### Graph Based KNN/ANN

Για το similarity search μέσω γράφων, έχουμε υλοποιήσει 2 αλγορίθμους, το GNN και το MRNG

###### GNN

Για τον GNN έχουμε δημιουργία του γράφου όπως φαίνεται και στις διαφάνειες δηλαδή πέρασμα από κάθε σημείο του γράφου και επιλογή των πλησιέστερων (`-k` το πλήθος) γειτόνων του που θα ενωθούν με το εκάστοτε σημείο. Μετά την κατασκευη του γράφου, η αναζήτηση γίνεται όπως ορίζουν πάλι οι διαφάνειες, δηλαδή ξεκινάμε `-R` φορές από ένα τυχαίο σημείο και κινούμαστε προς το query point, κάνοντας expand από τον κόμβο προς τους γείτονές του. Ο κώδικας μας συμπεριλαμβάνει κριτήριο τελειώματος στα T=50 extensions ή όταν δεν υπάρχει γείτονας ο οποίος να μας πηγαίνει πιο κοντά στο σημείο μας (Με ό,τι αυτό μπορεί να επιφέρει στο αποτέλεσμα αναζήτησης, αφού μπορεί να διακόψουμε το search νωρίς και να μην αφήσουμε την δυνατότητα να προχωρήσουμε σε έναν ενδεχόμενα καλύτερο κόμβο πιο μετά). Ο αριθμός των κόμβων που γίνονται extended απτους συνολικά διαθέσιμους είναι `-E`. Στο τέλος επιλέγουμε τους `-N` πιο κοντινούς με χρήση priority queue.


###### MRNG

Για τον MRNG ακολουθήσαμε ό,τι λένε οι διαφάνειες με μικρές αλλαγές που σημειώθηκαν και στο μάθημα για επιτάχυνση του αλγορίθμου. Για τον construction του γράφου στο `mrng_construction` χρησιμοποιήσαμε threads για να επισπεύσουμε την διαδικασία, και LSH για τον nearest neighbor του συνόλου $L_p$. Η αναζήτηση είναι η ίδια ακριβώς με αυτήν που ορίζεται στις διαφάνειες.


#### Μέρος Γ - Αποτελέσματα

Παρακάτω υπάρχουν τα αποτελέσματα από κάποια πειράματα που κάναμε (το AAF είναι το Average Approximation Factor από όλα τα queries). Τα αποτελέσματα είναι αυτά που λάβαμε μετά από την διεξαγωγή 25 queries.

| GNN(T=50,E=30) | R | Av. Time (ms)  | MAF     | AAF     |
| --- | - | -------- | ------- | ------- |
|  | 1 | 2.41667  | 3.82757 | 1.96841 |
| | 10 | 27.0417 | 2.02649 | 1.43181 |
|       | 100 | 243.667 | 1.70311 | 1.2058  |

| GNN(T=50,E=50) |    R   | Av.Time (ms)         | MAF     | AAF      |
| --- | ----- | -------- | ---- | ---- |
|  | 1   | 3.83333  | 3.28107 | 1.6981  |
|       | 10  | 46.25    | 1.74784 | 1.31598 |
|       | 100 | 382.583  | 1.4411  | 1.20215 |

| GNN(T=50,E=1) |    R   | Av.Time (ms)         | MAF     | AAF      |
| --- | ----- | -------- | ---- | ---- |
|  | 1   | 0  | 9.04017 | 3.26938 |
|       | 10  | 0.875    | 5.5801 | 3.26938 |
|       | 100 | 11.7917  |2.33704  | 1.56634 |

| MRNG | l         | Av. Time (ms) | MAF     | AAF     |
| ---- | --------- | ------- | ------- | ------- |
| | 20   | 0.0416667 | 6.79958 | 2.77082
| | 50   | 0.25      | 7.34335 | 2.40441 |         
| | 100  | 1.125     | 3.28764 | 1.54908 |        
| | 200  | 2.54167   | 1.54523 | 1.18795 |       
| | 500 | 9.45833 | 1.04168 | 1.00174
| | 1000 | 50.7917   | 1       | 1       |        

| LSH | w    | k  | l | Av. Time (ms) | MAF     | AAF     |
| --- | ---- | -- | - | ------- | ------- | ------- |
|     | 4    | 4  | 5 | 0       | 2.61661 | 1.91133 |
|     | 4    | 10 | 5 | 0       | 2.34288 | 2.34288 |
|     | 4    | 10 | 8 | 0.125   | 3.26639 | 2.13256 |
|     | 15   | 10 | 8 | 1.16667 | 2.16314 | 1.55292 |
|     | 30   | 10 | 8 | 1.875   | 1.68566 | 1.44499 |
|     | 100  | 10 | 8 | 6.25    | 1.52667 | 1.28124 |
|     | 1000 | 10 | 8 | 63.875  | 1.32238 | 1.11055 |

| HC     | w | k  | probes | M     | Av. Time (ms)  | MAF     | AAF     |
| ------- | ------- | -- | ------ | ----- | ------- | ------- | ------- |
|         | 500     | 4  | 2      | 10000 | 48.1667 | 1.16698 | 1.04131 |
|         | 4       | 4  | 2      | 10000 | 45.4583 | 1.6304  | 1.16797 |
|         | 4       | 10 | 2      | 10000 | 0.55    | 3.06635 | 1.72826 |
|         | 4       | 14 | 2      | 10000 | 1.95    | 6.03217 | 2.65682 |
|         | 4       | 8  | 4      | 10000 | 5.35    | 1.57319 | 1.30445 |


###### GNN

Όπως αναμέναμε, τα αποτελέσματα για τον GNN αλλάζουν ανάλογα με τις τιμές του R και του E. Βλέπουμε ότι για ένα σταθερό T, όταν το R είναι μικρό (δηλαδή 1 ή ακόμα και 10) έχουμε σχετικά κακά approximations σε σχέση με ένα πιο μεγάλο R, το οποίο είναι απόλυτα λογικό αφού δεν δίνουμε την δυνατότητα στην αναζήτηση να ξεκινήσει από αρκετά τυχαία σημεία ώστε να μπορέσουμε να αποφύγουμε την περίπτωση να πέσουμε σε ένα κακό αρχικό κόμβο. Όπως επίσης περιμέναμε ένα μικρό Ε δίνει επίσης κακά αποτελέσματα, αφού επιλέγουμε με αυτόν τον τρόπο πάντα τον 1ο KNN για expansion. Ακόμα και έτσι ωστόσο, για ένα μεγάλο R τα αποτελέσματα δεν είναι τόσο απογοητευτικά αφού πάλι δίνουμε την δυνατότητα στον αλγόριθμο να πέσει κοντά στο κόμβο αναζήτησης, ακόμα και αν κάνει expand λίγους κόμβους κάθε φορά.


###### MRNG

Όπως πάλι αναμέναμε, όσο το L, δηλαδή ο αριθμός των υποψήφιων μεγαλώνει τόσο καλύτερο γίνεται το search μας αφού δίνουμε την δυνατότητα στο μονοπάτι να προχωρήσει αρκετά, εκμεταλλευόμενοι ταυτόχρονα το ότι είναι και μονοτονικό.

###### LSH
Όπως φαίνεται καθαρά, το LSH και ο χρόνος του επηρεάζεται αρκετά από την τιμή του W (window size). Όταν αυτό είναι χαμηλά ο αλγόριθμος βγάζει αρκετά καλά αποτελέσματα, και αρκετά γρήγορα. Όσο μεγαλώνει, τόσο η απόδοση του γίνεται καλύτερη και ο χρόνος αναζήτησης μεγαλώνει. Όπως επίσης φαίνεται, η προσθήκη και άλλων hash functions ή hash tables για το LSH δεν επιφέρει πάντα καλύτερα αποτελέσματα.

###### Hypercube
Ο υπερκύβος αν και επηρεάζεται και αυτός απτό W, επηρεάζεται ίσως πιο σημαντικά από την επιλογή του k και των probes. Το αξιοσημείωτο είναι ότι το AAF για τον υπερκύβο δεν ξεφεύγει πολύ, αν και ο χρόνος μεταβάλλεται αρκετά. Ένα μεγάλο w με μικρά k και probes έχει οριακά ειδανικά αποτελέσματα αν και είναι αρκετά πιο αργό απτους υπόλοιπους συνδυασμούς. Βλέπουμε όμως ότι το AAF παραμένει σε σχετικά χαμηλές τιμές και για άλλους γρήγορους συνδυασμούς.

###### Συμπέρασμα

Ο MRNG φαίνεται να έχει τα καλύτερα αποτελέσματα από όλους τους αλγόριθμους σε συνδυασμό χρόνου και MAF/AFF. Σε μόλις 50ms προβαίνει στο να βρίσκει πάντα τον καλύτερο nearest neighbor. Οι LSH/HYPERCUBE/GNN είναι αρκετά κοντά του σε παρόμοιο time ωστόσο. Αξιοσημείωτο ωστόσο είναι το ότι για μικρούς χρόνους, οι αλγόριθμοι LSH/Hypercube φαίνεται να αποφέρουν καλύτερα αποτελέσματα.


#### Δομή


Στο directory `includes` υπάρχουν όλα τα header files που χρησιμοποιούνται.

 - `lsh_includes` υπάρχουν τα Includes για το LSH και τον Υπερκύβο για χρήση εντός των modules
 - `graph.hpp` περιέχει τον ορισμό και τον κώδικα για τον γράφο μας, ο οποίος χρησιμοποιεί templates για να είναι generic.
 - `graph_search.hpp` περιέχει την υλοποιήση για τους αλγορίθμους GNN και MRNG μέσω του Class Graph_Search_Ann.

 

Στο directory `modules/gnn` υπάρχουν εσωτερικά όλες οι υλοποιήσεις των modules.

Στο directory `obj` τοποθετούνται τα Object files.

Στο directory `outputs` τοποθετούνται by default αρχεία εξόδου.

Στο directory `data` τοποθετούνται τα αρχεία που διαβάζονται απτον χρήστη (αν θέλει) ή τα άλλα δεδομένα.

Στο directory `lib` υπάρχει η compiled σε library version του LSH/HYPERCUBE όπως παραδόθηκε για την Εργασία 1.

---

#### Compilation και εκτέλεση

Για compilation:

1) Για το compilation GNN/MRNG: `make graph_search`

Για quick runs:

1) Για το MRNG: `make run_mrng`

2) Για το GNN: `make run`

