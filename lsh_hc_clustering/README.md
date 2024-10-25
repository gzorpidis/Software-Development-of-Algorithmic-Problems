# Ανάπτυξη Λογισμικού για δυσεπίλυτα Αλγοριθμικά Προβλήματα
## Εργασία 1

---


##### Εκπόνηση από: 

1) Ζορπίδης Γεώργιος(1115202000055)

2) Κωνσταντίνος Φράγκος (1115202000207)

---

##### Αλλαγές απτό πρωτότυπο/επεκτάσεις:

Προσθέσαμε στο LSH,Hypercube το `-latent` flag το οποίο χρησιμεύει για να βγάλουμε αποτελέσματα με βάσει ένα αρχικό dataset αν τρέχουμε τους αλγόριθμους στο latent space. Χρειάζεται να προσδιορίσουμε επίσης το `-id` και `-iq` , τα αρχικά δηλαδή dataset και query sets. Για να δουλέψει ορθά πρέπει να υπάρχει πλήρης αντιστοίχιση των latent και initial datasets.


Προσθέσαμε στο Clustering την μέθοδο `_compute_latent_silhouettes(<initial_dataset>)` για τον υπολογισμό του silhouette με βάσει ένα αρχικό πάλι dataset.

##### Λίγα πράγματα για την υλοποίηση

Οι αλγόριθμοι για τους nearest neighbours (με LSH, Hypercube) και το clustering δημιουργήθηκαν για κύρια χρήση με το MNIST database, δηλαδή για εικόνες 28*28 grayscale, για αυτό και γίνεται εκτενής χρήση του τύπου `uint8_t`.

#### Similarity Search

Για το similarity search, έχουμε υλοποιήσει 2 αλγορίθμους, το Locality Sensitive Hashing (LSH) και το Randomized Projection σε υπερκύβο (Hypercube).

###### LSH

Η λογική του LSH βασίζεται στην διαμέριση του ευκλείδιου χώρου σε πολλά και διαφορετικά τυχαία (υπερ)επίπεδα. Για την υλοποίηση του αλγορίθμου αυτού, αρχικά δημιουργούμε τις βασικές L το πλήθος G hash functions `vector<g_hash> g_hash_functions`. Έπειτα έχουμε την δημιουργία μιας οικογένειας από τυχαίες h hash functions. Κάθε G hash function, επιλέγει k το πλήθος από αυτές τις h hash functions, τις οποίες χρησιμοποιεί για να παράξει ένα τυχαίο νούμερο ως γραμμικό συνδυασμό των h hash functions αυτών (σε συνδυασμό με ένα τυχαίο $r_i$) $g(p) = r_1h_1(p) + r_2h_2(p) + \cdots + r_kh_k(p) $ για κάθε input (εικόνα=vector). Στον υπολογισμό αυτό έχουμε όπως αναφέρεται και στις διαφάνειες, 2 modulo πράξεις (mod M, και ID(p) mod TableSize).  Για το search μέσα στα bucket, επειδή 2 items τα οποία είναι φαινομενικά πολύ μακριά μπορεί να συμπέσουν (εξαιτίας του mod TableSize), κάνουμε έλεγχο για να έχουνε ίδιο ID.

###### Hypercube

Βασική σκέψη για το randomized projection στον υπερκύβο είναι να κρατάμε τις τιμές που παίρνουν οι διάφορες h(p) και να τις κάνουμε map στο ${0,1}$ με σωστό τρόπο. Για να το επιτύχουμε αυτό χρησιμοποιούμε ένα mapping από integer σε ένα δεύτερο mapping (hypercube.hpp: line 58: `map<int, map<int,int>> fs`). Βασική της λειτουργία είναι να κρατάει τις αντιστοιχίες των τιμών ανάμεσα στις f. Για παράδειγμα το `fs[0]` αναφέρεται στο mapping το οποίο υπάρχει για τις τιμές της $h_1$. Όταν μια νέα τιμή της $h(p)$ η οποία δεν υπάρχει στο map εντοπιστεί, τότε κάνουμε flip ένα coin για να δούμε αν πρέπει να αντιστοιχιστεί σε 0 ή 1. Ωστόσο, αν η τιμή αυτή υπάρχει στο mapping, δηλαδή έχει πάρει τιμή, τότε δεν ξαναυπολογίζουμε τίποτα.

Επίσης, ως σχεδιαστική επιλογή πήραμε το να κάνουμε construct από τα 0 και 1 των επιμέρους $f_i$ τον ακέραιο που μας δείχνει την γωνία του υπερκύβου. Για παράδειγμα, αν το $k=3$ αυτό σημαίνει ότι ο υπερκύβος μας έχει $2^k = 8$ πιθανές τιμές για γωνίες. Έπειτα ανάλογα με τα επιμέρους 0 και 1, και τα αντίστοιχα shiftings, δημιουργούμε ένα ακέραιο που αντιστοιχείται σε μια γωνία του υπερκύβου ( Η υλοποίηση της παραπάνω λογικής γίνεται στο hypercube.hpp : line 105: `_evaluate()` σε άρρηκτο συνδυασμό με την `_map_point_to_f_values()` ενώ η ολοκληρωμένη λογική hashing γίνεται στην `_train_point()` )


#### Centroid-based Clustering

Για το clustering έχουμε 3 αλγόριθμους υλοποίησης:

1) Lloyd's Assignment

2) Reverse Assignment μέσω LSH

3) Reverse Assignment μέσω Hypercube

Για την υλοποίησή του, χρησιμοποιούμε 3 βασικές κλάσεις. Η `Cluster` είναι η βασική κλάση, την οποία ο χρήστης πρέπει να χρησιμοποιήσει για να κάνει cluster τα δεδομένα του. Η `Cluster_center` είναι μια κλάση η οποία χρησιμοποιείται για να επιδείξει το κέντρο ενός cluster. Η `cluster_item` είναι μια κλάση η οποία χρησιμοποιείται για το κάθε σημείο το οποίο μπαίνει στο αντίστοιχο `cluster_center`.

Η `cluster_item` κρατάει τα δεδομένα του κάθε σημείου (`_vector`), το id του σημείου(εικόνας) (`_id`) και 2 flags για το αν έχει μπεί σε κάποιο cluster (`_marked`) ή αν είναι να μπεί σε κάποιο cluster (`_queued`) (χρησιμοποιείται για την επίλυση των conflicts στo reverse assignment).

Το `Cluster_center` κρατάει ένα id το οποίο είναι της αρχικής εικόνας που έγινε assigned ως κέντρο, μέσω του kmeans++, τις συντεταγμένες του κεντροειδούς (σε double και uint8_t), πληροφορίες για την Silhouette, και πρωτίστως ένα vector από `cluster_items`, για τα σημεία τα οποία βρίσκονται assigned εκεί.

---

#### Δομή

Στο directory `includes` υπάρχουν όλα τα header files που χρησιμοποιούνται.

 - `Cluster.hpp` περιέχει τις κλάσεις του clustering
 - `common.hpp` περιέχει συχνά χρησιμοποιούμενο κώδικα
 - `get_long_opts_single.hpp` μια βιβλιοθήκη που αναπτύξαμε για να διαβάζει τα options απτην γραμμή εντολών
 - `headers.hpp` Ένα συνολικό header file αρχείο, που περιέχει πολλά includes, για να μην γίνονται συνέχεια
 - `hypercube.hpp` περιέχει τις κλάσεις για τον υπερκύβο
 - `lsh.hpp` περιέχει τις κλάσεις για το LSH
 - `structures.hpp` περιέχει συχνά χρησιμοποιούμενες κλάσεις, κοινές μεταξύ των υλοποιήσεων
 

Στο directory `modules` υπάρχουν εσωτερικά όλες οι υλοποιήσεις των modules και των αλγορίθμων.

Στο directory `objects` τοποθετούνται τα Object files

Στο directory `outputs.txt` τοποθετούνται by default αρχεία εξόδου

Στο directory `data` τοποθετούνται τα αρχεία που διαβάζονται απτον χρήστη (αν θέλει)

---

#### Compilation και εκτέλεση

Για compilation:

1) Για το LSH: ``` make lsh ```

2) Για το Hypercube: ```make cube ```

3) Για το Clustering: ```make cluster```

Για quick runs:

1) Για το LSH: `make run_lsh`

2) Για το Hypercube: `make run_cube`

3) Για το Clustering: `make run_clustering_classic ` ή `make run_clustering_lsh` ή `make run_clustering_hc`

#### Παρατηρήσεις:
Έστω ε το σφάλμα της απόστασης του πραγματικού πλησιέστερου γείτονα, με αυτό που βρίσκουν οι αλγόριθμοί μας.
Όσον αφορά το lsh, ανεβάζοντας την τιμή του window, παρατηρούμε ότι παίρνουμε μικρότερο ε, λογικό καθώς τόσο μεγαλύτερο
είναι το w, τόσο περισσότερο συμπιέζοντε τα σημεία πάνω σε μία ευθεία h με αποτέλεσμα να έρχοντε πιο κοντά και να μπαίνουν
περισσότερα σημεία σε ένα κουβά. Μεγαλώνοντας το window όμως χάνουμε σε χρόνο καθώς πρέπει να ελέγχουμε περισσότερα σημεία.
Το ίδιο ισχύει και για το L όπου για μεγαλύτερο L αυξάνουμε την πιθανότητα να βρούμε τον πλησιέστερο γείτονα και
ελαχιστοποιούμε το ε, όμως χάνουμε σε χρόνο καθώς κάνουμε περισσότερη δουλειά, αφού αυξάνοντε οι g.
Ο χρόνος για την εξαντλητική μέθοδο είναι 287ms, ενώ για την lsh ο χρόνος με χαμηλές τιμές των L και w, μπορεί
να φτάσει και τo 1ms. Βέβαια όσο μεγαλώνουμε αυτές τις παραμέτρους μεγαλώνει και ο χρόνος.

Παρόμοια, στο hypercube με την παράμετρο Μ εξετάζουμε παραπάνω σημεία, δηλαδή μικρότερο ε και μεγαλύτερος χρόνος,
η παράμετρος k όσο μεγαλύτερη είναι, χωρίζει σε περισσότερα τμήματα τα σημεία μας με αποτέλεσμα όσο μεγαλύτερη είναι
τόσο πιο γρήγορη είναι η μέθοδός μας, όμως λιγότερο αποτελεσματική. Η probes προφανώς όσο μεγαλύτερη είναι, τόσο 
πιο πιθανόν να έχουμε μικρότερο ε καθώς εξετάζουμε περισσότερες γωνίες του υπερκύβου, χάνοντας στο χρόνο όμώς.
Ο χρόνος μπορεί να πέσει κάτω από το 1ms σε αυτή την υλοποίηση, όμως για να έχουμε αρκετά αποτελεσματική 
υλοποίηση, τροποποιούμε τις παραμέτρους μας με τέτοιο τρόπο ώστε να να έχουμε χαμηλό χρόνο το πολύ 80ms
πχ( m = 10000, k = 4, w = 500, probes = 2).