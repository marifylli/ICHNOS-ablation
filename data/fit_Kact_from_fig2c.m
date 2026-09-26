% ------------------------------------------------------------------------
% Προσαρμογή των K_act_ox και n_ox στα δεδομένα του Delaunay Fig. 2C.
%
% Στόχος: το A_ox στα t = 5 min (όχι η μόνιμη κατάσταση) - το panel C
% μετρήθηκε μετά από 5 λεπτά έκθεσης.
%
% Τα k_off_ox και d_x_ox κρατιούνται σταθερά. ΠΡΟΣΟΧΗ: οι τιμές 150 / 0.6
% βαθμονομήθηκαν υπό K_act_ox = 271. Μετά από αυτό το fit πρέπει να
% ξανατρέξει το sweep R x d_x.
%
% Το CSV πρέπει να είναι στον τρέχοντα φάκελο.
% ------------------------------------------------------------------------

COLUMN = 'ox_fraction_window';   % τρέξε ΚΑΙ με 'ox_fraction_unmix'
kOff   = 150;                    % από A.9
dX     = 0.6;                    % από A.9
p0     = [300, 2];               % αρχική εικασία από τη densitometry

% ------------------------------------------------------------------------
m  = sbioselect('Type','sbiomodel','Name','OxidativeModule');
if isempty(m)
    error('Το μοντέλο δεν βρέθηκε. Φόρτωσε πρώτα το SimBiology project.');
end

sf = createSimFunction(m, ...
        {'K_act_ox','n_ox','S_ox','k_off_ox','d_x_ox'}, {'A_ox'}, []);

T     = readtable('delaunay_fig2c_oxfraction.csv');
doses = T.dose_uM;
y     = T.(COLUMN);
yn    = (y - min(y)) / (max(y) - min(y));

t5  = 5/60;
obj = @(p) objective(p, sf, doses, kOff, dX, t5, yn);

opt = optimset('Display','iter','TolX',1e-3,'TolFun',1e-6);
pf  = fminsearch(obj, p0, opt);

fprintf('\n--- %s ---\n', COLUMN);
fprintf('K_act_ox = %.1f uM\n', pf(1));
fprintf('n_ox     = %.2f\n',    pf(2));

A5 = predictA5(pf, sf, doses, kOff, dX, t5);
figure; semilogx(doses, yn, 'o', 'MarkerFaceColor','k'); hold on
semilogx(doses, (A5-min(A5))/(max(A5)-min(A5)), 's-');
xlabel('H_2O_2 (\muM)'); ylabel('normalised oxidised fraction');
legend('Delaunay Fig. 2C','model, A\_ox at 5 min','Location','southeast');
grid on; title(sprintf('K_{act} = %.0f \\muM, n = %.2f', pf(1), pf(2)));

fprintf('\ndose  data   model\n');
An = (A5-min(A5))/(max(A5)-min(A5));
for i = 1:numel(doses)
    fprintf('%4d  %.3f  %.3f\n', doses(i), yn(i), An(i));
end

% ------------------------------------------------------------------------
% Τοπικές συναρτήσεις - ΠΡΕΠΕΙ να είναι στο τέλος του αρχείου
% ------------------------------------------------------------------------
function A5 = predictA5(p, sf, doses, kOff, dX, t5)
    A5 = zeros(numel(doses),1);
    for i = 1:numel(doses)
        [~, out] = sf([p(1) p(2) doses(i) kOff dX], [], [], [0; t5]);
        A5(i) = out{1}(end);
    end
end

function e = objective(p, sf, doses, kOff, dX, t5, yn)
    if p(1) <= 0 || p(2) <= 0, e = 1e6; return; end
    A5 = predictA5(p, sf, doses, kOff, dX, t5);
    r  = max(A5) - min(A5);
    if r < 1e-9, e = 1e6; return; end
    e  = norm( (A5 - min(A5))/r - yn );
end
