"""
Densitometry του Delaunay et al. 2000 (EMBO J), Figure 2C.
Καμπύλη δόσης οξείδωσης Yap1 in vivo, 5 λεπτά έκθεσης, 25-800 uM H2O2.

Το panel C μετράει ΜΕΤΑΤΟΠΙΣΗ ΜΠΑΝΤΑΣ, όχι ένταση: η οξειδωμένη μορφή έχει
ενδομοριακούς δισουλφιδικούς δεσμούς, είναι πιο συμπαγής, και τρέχει χαμηλότερα
στο non-reducing gel. Η συνολική ένταση ανά λωρίδα μετράει μόνο τη φόρτωση
(εδώ κυμαίνεται 2950-6580, δηλαδή >2x), οπότε ΠΡΕΠΕΙ να δουλέψουμε με λόγο
εντός λωρίδας.

Έξοδος: οξειδωμένο κλάσμα ανά δόση, με δύο ανεξάρτητες μεθόδους.
"""
import numpy as np
from PIL import Image
from scipy.optimize import nnls, curve_fit

IMG    = 'fig2c_crop.png'          # crop μόνο του gel, χωρίς τους αριθμούς δόσης
BORDER = 8                          # pixels πλαισίου προς αφαίρεση
DOSES  = np.array([25, 50, 100, 150, 200, 300, 800])
GAPS   = [145, 290, 440, 570, 710, 850]   # κατά προσέγγιση κέντρα των κενών
RED_W  = slice(55, 100)             # παράθυρο ανηγμένης μορφής (πάνω μπάντα)
OX_W   = slice(105, 152)            # παράθυρο οξειδωμένης μορφής (κάτω μπάντα)
INSET  = 6                          # περιθώριο εντός κάθε λωρίδας

img  = np.array(Image.open(IMG).convert('L')).astype(float)
core = (255.0 - img)[BORDER:-BORDER, BORDER:-BORDER]   # αντιστροφή: σκούρο = σήμα

# --- εντοπισμός ορίων λωρίδων από τα ελάχιστα του προφίλ στηλών ---
px = core.sum(axis=0)
sm = np.convolve(px, np.ones(9) / 9, mode='same')
edges = [5] + [int(np.argmin(sm[c-25:c+25]) + c - 25) for c in GAPS] + [len(px) - 5]

# --- προφίλ καθ' ύψος ανά λωρίδα, με τοπικό background ---
prof = {}
for i, (l, r) in enumerate(zip(edges[:-1], edges[1:])):
    lane = core[:, l + INSET : r - INSET]
    bg   = np.median(np.r_[lane[:25, :], lane[-25:, :]])
    prof[DOSES[i]] = (lane - bg).clip(0).mean(axis=1)

# --- Μέθοδος 1: σταθερά παράθυρα γραμμών (ίδια για ΟΛΕΣ τις λωρίδες) ---
m1 = np.array([prof[d][OX_W].sum() / (prof[d][RED_W].sum() + prof[d][OX_W].sum())
               for d in DOSES])

# --- Μέθοδος 2: αποσύνθεση σε δύο πρότυπα (ακραίες λωρίδες ως βάση) ---
T  = np.c_[prof[DOSES[0]], prof[DOSES[-1]]]
m2 = np.array([(lambda w: w[1] / w.sum())(nnls(T, prof[d])[0]) for d in DOSES])

# --- Hill fit ---
hill = lambda S, K, n: (S / K) ** n / (1 + (S / K) ** n)
print(' dose   ox_frac(M1)  ox_frac(M2)')
for d, a, b in zip(DOSES, m1, m2):
    print(f'{d:5d}      {a:.3f}        {b:.3f}')
for name, y in (('Μέθοδος 1', m1), ('Μέθοδος 2', m2)):
    yy = (y - y.min()) / (y.max() - y.min())
    K, n = curve_fit(hill, DOSES, yy, p0=[150, 2], maxfev=20000)[0]
    print(f'{name}: K = {K:.0f} uM, n = {n:.2f}')

# ΠΡΟΣΟΧΗ: αυτά είναι ΕΝΕΡΓΑ μεγέθη του A_ox στα 5 min, ΟΧΙ οι παράμετροι
# K_act_ox / n_ox της συνάρτησης Hill μέσα στο rate rule. Ο βρόχος A_ox/X_ox
# παρεμβάλλεται. Χρησιμοποίησέ τα ως αρχική εικασία και προσάρμοσε το ΜΟΝΤΕΛΟ
# στα επτά σημεία (βλ. fit_Kact_from_fig2c.m).
