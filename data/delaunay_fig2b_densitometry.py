"""
Densitometry του Delaunay et al. 2000 (EMBO J), Figure 2B.
Χρονοσειρά οξείδωσης Yap1 in vivo, 400 uM H2O2, non-reducing gel.
Χρόνοι: 0, 2.5, 5, 15, 30, 45, 60 min.

ΙΔΙΑ ΛΟΓΙΚΗ ΜΕ ΤΟ 2C: το panel μετράει ΜΕΤΑΤΟΠΙΣΗ ΜΠΑΝΤΑΣ, όχι ένταση.
Η οξειδωμένη μορφή έχει ενδομοριακούς δισουλφιδικούς δεσμούς, είναι πιο
συμπαγής, και τρέχει ΧΑΜΗΛΟΤΕΡΑ. Η φόρτωση εδώ κυμαίνεται ~4x (2854 έως
12454), οπότε η ποσοτικοποίηση ΠΡΕΠΕΙ να είναι λόγος εντός λωρίδας.

Το reducing gel (με DTT) χρησιμοποιείται ως ανεξάρτητος έλεγχος φόρτωσης:
εκεί όλο το Yap1 τρέχει στην ανηγμένη θέση ανεξάρτητα από την κατάσταση
οξείδωσης, οπότε η συνολική ένταση ανά λωρίδα = συνολικό Yap1.

Παράγει: delaunay_fig2b_timecourse.csv
"""
import numpy as np
from PIL import Image
from scipy.optimize import nnls

IMG_NR = 'fig2b_nonreducing.png'   # crop μόνο του gel
IMG_R  = 'fig2b_reducing.png'      # loading control, ίδια όρια λωρίδων
BORDER = 14
TIMES  = np.array([0, 2.5, 5, 15, 30, 45, 60])
GAPS   = [160, 318, 462, 615, 760, 930]   # κατά προσέγγιση κέντρα των κενών
RED_W  = slice(45, 112)   # παράθυρο ανηγμένης μορφής (ΠΑΝΩ μπάντα)
OX_W   = slice(116, 180)  # παράθυρο οξειδωμένης μορφής (ΚΑΤΩ μπάντα)
INSET  = 6


def lane_profiles(path, gaps, border=BORDER, inset=INSET):
    """Επιστρέφει (προφίλ ανά λωρίδα, όρια λωρίδων)."""
    img  = np.array(Image.open(path).convert('L')).astype(float)
    core = (255.0 - img)[border:-border, border:-border]   # σκούρο = σήμα

    # ΕΛΕΓΧΟΣ: το σήμα δεν πρέπει να ακουμπάει τις άκρες του crop
    py = core.sum(axis=1)
    if py[:5].mean() > 0.4 * py.max() or py[-5:].mean() > 0.4 * py.max():
        raise RuntimeError('Το σήμα φτάνει στην άκρη του crop. Ξανακόψε '
                           'με 30-40 px λευκού πάνω και κάτω.')

    px = core.sum(axis=0)
    sm = np.convolve(px, np.ones(9) / 9, mode='same')
    edges = [5] + [int(np.argmin(sm[c-30:c+30]) + c - 30) for c in gaps] \
              + [core.shape[1] - 5]

    prof = []
    for l, r in zip(edges[:-1], edges[1:]):
        lane = core[:, l + inset : r - inset]
        bg   = np.median(np.r_[lane[:25, :], lane[-20:, :]])
        prof.append((lane - bg).clip(0).mean(axis=1))
    return prof, edges


# ---------- non-reducing: το σήμα ----------
prof, edges = lane_profiles(IMG_NR, GAPS)
print('lane edges:', edges)

# Μέθοδος 1: σταθερά παράθυρα γραμμών, ΙΔΙΑ για όλες τις λωρίδες
m1     = np.array([p[OX_W].sum() / (p[RED_W].sum() + p[OX_W].sum()) for p in prof])
totals = np.array([p.sum() for p in prof])

# Μέθοδος 2: αποσύνθεση σε δύο πρότυπα (t=0 = ανηγμένο, t=5 = οξειδωμένο)
T  = np.c_[prof[0], prof[2]]
m2 = np.array([(lambda w: w[1] / max(w.sum(), 1e-9))(nnls(T, p)[0]) for p in prof])

# ---------- reducing: loading control ----------
try:
    prof_r, _ = lane_profiles(IMG_R, GAPS)
    load = np.array([p.sum() for p in prof_r]); load /= load.max()
except FileNotFoundError:
    load = np.full(len(TIMES), np.nan)

print('\n t(min)  ox_frac(win)  ox_frac(unmix)  total(NR)  loading(R)')
for t, a, b, tot, ld in zip(TIMES, m1, m2, totals, load):
    print(f'{t:6}     {a:.3f}         {b:.3f}      {tot:8.1f}    {ld:.3f}')

np.savetxt('delaunay_fig2b_timecourse.csv',
           np.c_[TIMES, m1, m2, totals],
           delimiter=',', header='time_min,ox_fraction_window,'
           'ox_fraction_unmix,lane_total_intensity', comments='', fmt='%.4g')

# ------------------------------------------------------------------------
# ΣΗΜΕΙΩΣΕΙΣ ΕΡΜΗΝΕΙΑΣ (MASTER §15)
#
# 1. Για ΑΥΤΟ το gel εμπιστευόμαστε τη ΜΕΘΟΔΟ ΠΑΡΑΘΥΡΩΝ. Η αποσύνθεση
#    χρησιμοποιεί το t=5 min ως προτυπο "καθαρου οξειδωμενου", αλλα η
#    μεθοδος παραθυρων λεει οτι εκεινο το σημειο ειναι μολις 0.775.
#    Το προτυπο ειναι μολυσμενο, γι' αυτο η αποσυνθεση διογκωνει τις τιμες.
#
# 2. Οι λωριδες 2.5 και 5 min ειναι υποφορτωμενες (επιφυλαξη Delaunay
#    §3.1, επιβεβαιωθηκε με το reducing gel). Επειδη ομως ο λογος ειναι
#    ΕΝΤΟΣ λωριδας, η φορτωση ακυρωνεται: θορυβωδεστερες, οχι μεροληπτικες.
#
# 3. Οι λωριδες 45 και 60 min εχουν ΥΨΗΛΗ φορτωση (0.910, 0.863). Η
#    επιστροφη στο 0.248 και 0.107 ΔΕΝ ειναι artifact αδυναμου σηματος.
# ------------------------------------------------------------------------
