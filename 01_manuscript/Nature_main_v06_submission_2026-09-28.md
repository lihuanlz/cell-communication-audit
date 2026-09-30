# Scale-free statistics are the currency of cellular communication
Huan Li

Shanghai Liangta Biotech Co., Ltd., Shanghai 201899, P. R. China

Corresponding author: HL@liangtabio.com (H.L)


## Abstract

Cells compute with biochemical circuits whose component concentrations vary from cell to cell and cannot be calibrated against any shared reference. We recently showed that under such multiplicative scale uncertainty, assays can transmit only scale-free statistics (counts, ratios, event times); here we ask whether cells themselves communicate this way. We analyse four signalling systems under one discipline: derive each model's degeneracy group, identify the statistics invariant under it, and locate those statistics in data. A covalent-modification sensor is structurally blind to all absolute molecular scales. The p53 pulse encoder places 94.2% of its parameter information on timing rather than amplitude (a calibration-bounded decomposition whose ordering, not its exact value, is the robust result), and its dose response is a counting law. In bacterial chemotaxis the amplitude description is internally inconsistent while the scale-free Weber relation is exact and prospectively predictive. In GPCR pharmacology, 21.5% of archived potency parameters fail the field's own algebraic identities. Both known routes to calibration-freeness, chemical saturation and circuit topology, collapse the likelihood onto counting statistics in the same way. Across four systems, absolute-scale descriptions are unidentifiable, inconsistent or unused; scale-free statistics are the currency of cellular communication.

---

## Introduction

A cell that reports a number to itself or to its neighbours faces a measurement problem no instrument engineer would accept. The components of its signalling circuits (kinases, receptors, transcription factors) are present in counts that vary by tens of per cent from cell to cell and drift over hours; no two cells share a calibration standard, and no external reference exists against which an intracellular concentration can be read out absolutely. Yet cells routinely transmit dose, identity and urgency with quantitative precision. How?

We recently formalised the abstract version of this problem for bioanalytical assays. In the Target-Capture Scaling (TCS) framework, a measurement whose observable is a multiplicative function of uncalibratable scale factors possesses a degeneracy group: the family of parameter transformations that leave every observable invariant. What survives the group, and is therefore reportable without calibration, are scale-free statistics: counts, ratios, ordinals and event times. Assay platforms that report absolute analogue scale must pay for calibration; platforms that count (digital PCR being the limiting case^1^) are calibration-free by construction. The companion preprint^2^ establishes this taxonomy and its payment theorem on the assay side.

The present paper transfers the question from assays to cells: **if the intracellular environment is uncalibratable in absolute scale, what statistics can cellular communication actually use, and do real signalling systems in fact use them?** We answer with analysis, not argument. The analysis is a fixed five-step discipline, applied identically to every system (Fig. 1a): write the input–output relation in dimensionless form; derive the degeneracy group exactly; stratify parameter directions by Fisher information into identifiable, sloppy and degenerate tiers; replace point estimates with one-sided bounds where likelihoods plateau; and name the orthogonal measurement that would break each degeneracy. The analysis returns, for each system, a stratification of what any experiment on its standard observables can in principle return, and therefore which statistics are candidates for carrying information through an uncalibratable interior.

We apply the analysis to four systems that together span the organisational levels of cellular signalling: a static covalent-modification sensor (the phosphorylation–dephosphorylation cycle, PdPC), a dynamic pulse encoder (the p53 DNA-damage response), a chemosensory receptor layer (bacterial chemotaxis, in the best public single-cell dose–response data), and a pharmacological meta-level (GPCR concentration–response parameters archived across a field database). Two unifying layers complete the programme: a theorem-level analysis showing that the two known routes to calibration-freeness (chemical saturation and circuit topology) collapse the likelihood onto counting statistics in the same way (the digital limit); and a set of seven pre-registered blind adjudications on six independent public single-cell datasets that test the framework's extrapolations and record its falsifications verbatim. All analyses were executed under one bookkeeping discipline (pre-registration, frozen verdict cards, and double-recording of successes and failures) described in Methods.

The result is a single pattern with four independent instances: in every system, the absolute-scale description is either structurally invisible, internally inconsistent, or demonstrably unused, while scale-free statistics (counts, event times, ratios, midpoints) are identifiable, internally consistent, and carry the information. Cellular communication, like assay engineering, runs on what cannot be miscalibrated.

## Results

### The framework: degeneracy groups, surviving statistics, and the price of absolute scale

All four analyses below are instances of one mathematical object, which we state once here (full formalism in SI S1). Model a signalling pathway as a cascade of layers $L_1 \to \dots \to L_n$, where layer $i$ with parameters $\theta_i$ maps its upstream input to an output symbol $Y_i$, with biochemical noise absorbed into the layer's conditional distribution. Each layer carries a **degeneracy group** $G_i$: the set of parameter transformations that leave the output distribution invariant for every input. What any downstream measurement can in principle learn about the layer is exactly the content invariant under $G_i$: the group's orbit invariants; everything else is gauge. A layer-boundary rule makes this operational: parameters that are indistinguishable from the same output belong to the same layer, so information structure, not biochemical topology, defines the layers (feedback included) (Extended Data Fig. 1).

Within a layer, an input statistic $X$ (a fold-change, a sign, a count, an event time, an amplitude) either survives or leaks. Formally, $X$ is in the layer's **survival set** if the distribution of the transmitted alphabet $a_i(Y_i)$ it induces is invariant under the group action up to the noise floor (total-variation distance below a noise-calibrated $\varepsilon$ over the physiological working range). The survival set is computable, and for the generic case of an unidentified affine gain it has a sharp algebraic content: the unconditional survivors are difference ratios (fold-changes with zero baseline), signs, counts and event times; amplitude survives only if the gain is pinned from outside; increments survive only where ratios survive and the gain is pinned to the precision the downstream consumer requires. A statistic crosses an interface if and only if three legs hold simultaneously: algebraic survival ($X$ in the survival set), symbolic reachability ($I(X; Y_i) > 0$), and temporal bandwidth (the statistic's timescale lies inside the interface's measured perturbation–response bandwidth). Cascades compose by an exact information-loss identity: writing $L_i = I(S; Y_{i-1}) - I(S; Y_i)$ for the loss at layer $i$ of a Markov chain $S \to Y_1 \to \dots \to Y_n$, one has $H(S) = I(S; Y_n) + \sum_i L_i$ with monotone mutual information, so the end-to-end transferable content of any statistic is limited by the weakest layer, the layer with the largest loss, which the identity locates without any capacity concept. The framework's status is stated precisely: the information-loss identity and monotonicity are theorems given the model; the survival set and three-leg criterion are operational definitions; no tightness or capacity claim is made.

Payment enters through references. Quantities a pathway can use fall into three classes: class I, calibration-free survivors (ratios, signs, counts, event times; no reference needed); class II, structural constants set by evolution and maintained by physics ($K_d$, $K_m$: dimensional converters that trigger no signature); and class III, paid references (set points, internal standards; actively maintained by the cell at continuing energetic cost, diagnosed by operational markers: negative feedback acting on the reference itself, recovery of the reference after direct perturbation, drift of the reference when energy supply is inhibited). The markers are scored independently of the framework's own predictions, so a class-III assignment cannot be invoked after the fact to rescue a failed prediction. The **payment theorem** inherited from the assay side then reads: wherever a function references absolute scale, a class-III reference must be established and maintained at that position; the positions are determined by functional demand, not by where information is lost. Three prospective falsifiable predictions follow (SI S1.6): deleting a class-III maintenance mechanism must abolish absolute-scale sensitivity while sparing relative-scale sensitivity (P-1); a synthetic circuit with no active internal standard must fail to transmit stable absolute scale (P-2); and positions that functionally require absolute position but show no reported maintenance machinery constitute an exhaustive search list whose hits, or exhaustively documented absence, adjudicate the framework (P-3).

The four system analyses below instantiate the framework at four organisational levels. In each case the analysis (i) derives the degeneracy group exactly, (ii) computes the survival set or its Fisher-information image, and (iii) checks which class the system's working statistics belong to, and in two cases goes further, testing prospectively whether the framework's predictions survive contact with data it has never seen.

### The static sensor is blind to its own absolute scale

The phosphorylation–dephosphorylation cycle, the covalent-modification motif at the heart of kinase signalling, converts an input activity ratio into a steady-state modified fraction. With six physical parameters (two catalytic rates, the Michaelis constants, two enzyme totals, substrate total), its steady-state input–output curve is the standard object fitted to dose–response data. The first step is a demonstration: four parameter sets differing in every individual parameter by up to five-fold produce steady-state curves that are not similar but identical to the last digit of double-precision arithmetic (max $|\Delta u| = 0$ over the full titration range). No experiment on this curve, at any precision, distinguishes them. What the curve can return is not six parameters but exactly what the four sets share: two dimensionless combinations: the activity ratio $\xi$ and the saturation ratio $\kappa = K_m/S_T$.

This is Theorem 1 of the study: the steady-state curve depends on the six physical parameters only through these two groups; the degeneracy group is four-dimensional, and all absolute molecular scales are structurally unidentifiable from steady-state fitting. The asymmetric-Michaelis extension raises the identifiable set to the triple $(\xi, \kappa_1, \kappa_2)$ with the closed-form local slope $n_{\mathrm{H}} = 4/(4 - 1/(\kappa_1+\tfrac{1}{2}) - 1/(\kappa_2+\tfrac{1}{2}))$, and leaves the group four-dimensional; embedding the cycle in a positive-feedback switch (critical gain $g^{*} = 0.443$ at $\kappa = 0.05$) moves the operating point but not the group; the same four actions leave every branch of the bistable response invariant to machine precision. Theorem 2 quantifies the barrier at realistic noise: the shape tier ($\ln c$, $\ln \kappa$) is well conditioned (condition number $\chi \approx 30$) while the absolute-scale block carries an exact zero eigenvalue (effective $\chi \sim 10^{300}$; the quoted value is the deep-saturation extrapolation, with the archived code returning up to $\chi \sim 10^{33}$ at the deepest $\kappa$ values in double precision). No increase in sample size or precision crosses this barrier; only an orthogonal measurement axis, substrate or enzyme titration, does. Theorem 3 supplies the reporting protocol for the one informative boundary: the local Hill slope obeys the closed form $n_{\mathrm{H}} = 1 + 1/(2\kappa)$, so $\kappa$ is recoverable from curve shape, but below $\kappa \lesssim 0.01$ the confidence width diverges and the honest statement is a scan-based one-sided bound $\kappa \leq \kappa_{\max}$, never a spurious point estimate (Fig. 2; Extended Data Fig. 2).

The analysis also locates where the cycle should sit on its own $\kappa$-axis if it serves as a sensor. Steady-state mutual information between a log-normal input ensemble and the noisy output has an interior optimum ($\kappa^{*} \approx 0.05$–$0.1$ for narrow input distributions, 2.85 bits; $\kappa^{*} \approx 1$–$3$ for wide ones, 3.04 bits), while the celebrated zero-order ultrasensitive limit^3^ is information-poor for all input widths (1.1–1.4 bits): per snapshot, the infinite-gain switch is a decision element, not a measurement element. The static sensor's contribution to the main line is thus already double-edged: its absolute scale is invisible in principle, and its information-optimal operating point is an interior, scale-relative coordinate, not an absolute one.

### The dynamic encoder carries its information in time, not amplitude

The p53 tumour suppressor responds to DNA double-strand breaks with discrete, stereotyped pulses whose amplitude, duration and first-pulse timing are approximately dose-independent while their number increases with dose.^4–7^ The analysis of this encoder proceeds at four levels, from circuit topology down to the molecular layer; Table 1 collects the four layers with their frozen verdict values (Fig. 3; Extended Data Figs. 3–5, 14, 15).

**Table 1 | The p53 encoder analysis in four layers.** Values quoted verbatim from the frozen verdict cards; unit-level clauses, conventions and scoring scripts: SI S3, S4, S13, S14.

| Layer | Question | Finding | Key values |
|---|---|---|---|
| Topology selection (Theorem 4) | Do pulse observables select the circuit class? | Yes: three output-side signatures, onset continuity, hysteresis width and the range sensitivity $\alpha = \partial\ln A/\partial\ln D$, uniquely assign universality class at ${\sim}15\%$ amplitude precision; published p53 statistics sit deep in the excitable region, excluding pure negative feedback as sufficient architecture | NF circuits destabilise generically through a Hopf onset (Thomas–Soulé theorem^8^; Lemma 1): $A^2 \propto (s - s_H)$ at $R^2 = 0.9926$, hysteresis 0.275 vs normal-form 0.25; excitable onset is all-or-none (Lemma 2): 6.5% variation over 4.6-fold drive, Morris–Lecar control $\alpha = -0.258$; data: $|\alpha| \lesssim 0.03$ over the 100-fold $\gamma$ range; $\alpha_{\mathrm{NCS}} = -0.001 \pm 0.113$ vs $\alpha_{\mathrm{UV}} = +0.981 \pm 0.114$ (bootstrap $\times 20000$), 6.1 combined standard deviations apart (8.6 under the pre-registered one-sided convention, SI S3.3) |
| Channel census | Which observables carry the parameter information? | Timing channels dominate and amplitude is near-empty; population gain dispersion then kills the residual amplitude channel while timing stays pinned | 94.2% of Fisher information on timing (width 60.2%, count 18.2%, period 11.5%, first pulse 4.4%) vs 5.5% on amplitude; mechanism-level recheck on the 23-parameter circuit^10^ preserves the ordering (45.1% nominal; SI S13.3); at gain dispersion 0.6, $\mathrm{CV}(A)$ rises $0.002 \to 0.639$ while $\mathrm{CV}(T)$ stays pinned at $0.016$ |
| Noise structure | Where does fluctuation information live? | Noise amplitude $\sigma$ is an exact Fisher null direction of first-moment statistics; two dispersion signatures locate an independent stochastic sensing stage upstream of the oscillator core | $\lambda_7/\lambda_1 \leq 4.7 \times 10^{-19}$; $\mathrm{Var}(t_1)/\mathrm{Var}(T)$ bounded at 0.03–0.12 for any bare oscillator vs literature 5.8, reproduced as 5.7 at delay dispersion $\mu_d \approx T$; measured $\mathrm{IPI}_{\mathrm{CV}} = 0.30$ requires fast intra-cellular stochasticity (slow heterogeneity pins it at 0.083) |
| Counting channel | What reaches the nucleus? | Dose is transmitted near-losslessly as pulse count; the molecular layer closes underneath the observables | $I(D;N) = 3.07$ of 3.17 bits (96.8%); period from the dominant Jacobian eigenvalue 5.95 h vs simulated 5.48 h; logistic-substitution closed form 5.60 h vs measured 5.5 h; counting-law slope inverts to the repair timescale $\tau_r \approx 7.7$ h (SI S13.4) |

**The encoder’s design logic.** Read off the analysis, it is the main line’s second instance: dose information is placed in an alphabet (counts, intervals, timing)^11,12^ that is invariant under the multiplicative scale group of the intracellular environment, because the amplitude alphabet is both unidentifiable and information-free. The amplitude channel dies twice, geometrically at the single-cell layer and statistically at the population layer, while the exact zero direction of the Fisher spectrum equals the degeneracy-group orbit tangent and the structural parameters attain the Cramér–Rao bound: the three predictions of the conditional digital-limit theorem are verified numerically on this circuit. The same protein switches class with the damage type, excitable dose-counting pulses after double-strand breaks but a graded, non-excitable single pulse after UV^9^. The period-from-constants direction has prior art in simplified cores^13–15^; the six-species phase-share attribution and its recorded trail are contributed here (SI S14).

### Two roads to one digital limit

The static and dynamic analyses each isolate a regime in which measurement becomes calibration-free, reached by opposite means: in the covalent cycle by chemical saturation ($\kappa \to 0$), in the pulse encoder by circuit topology ($\alpha \to 0$). We prove constructively for both circuits that when the analogue scale factor cancels from the likelihood, the inferable information degenerates to counting statistics: the Fisher matrix retains a single nonzero eigenvalue along the partition-counting direction, curvature along the absolute-scale direction vanishes, the residual stochasticity is combinatorial rather than parametric, and only dimensionless ratios survive (Fig. 4). We call this shared structure the digital limit; its unconditional generalisation (a stratified classification theorem, a universal functional theorem, a noise-robustness bound) is given in SI S5 (Extended Data Fig. 6).

Within its stated domain (single channel, single readout, memoryless, passive sensing), counting is the unique calibration-free channel that costs no extra resource. Fold-change, Weber-law and ratiometric schemes are established calibration-free analogue strategies^16^, but all pay extra resources (temporal memory, multi-point sampling, a second channel) and cancel mainly the common-mode gain direction; a fold-change of counts is still counting. The perspective reorganises two celebrated results: the zero-order ultrasensitive switch is the digital limit of the static channel, an information-poor (${\approx}1$ bit per snapshot) decision element; the p53 pulse train is the digital limit of the dynamic channel, a measurement element encoding dose in counts. Between the two limits lies the analogue regime: information-rich (the interior optimum $\kappa^{*}$) but calibration-hungry, usable only where the degenerate scales are pinned externally. The resulting conjecture, if a channel's likelihood-relevant partition statistics are invariant under the degeneracy group of the intracellular environment, they must be of counting type, is stated precisely enough to be falsified, with its known non-counterexamples and its measurement-theoretic boundary in SI S5. Going digital has costs: near-critical maintenance is energetically continuous, and counting at low copy numbers is Poisson-limited; analogue grading wins below roughly 3–4 resolvable levels, counting above. This trade-off rationalises the observed distribution of regimes: p53, facing catastrophic events that demand high-fidelity decisions, counts; metabolic sensing that requires continuous fine-tuning stays analogue.

### Blind adjudications bound the claim empirically

A framework that explains everything explains nothing, so we pre-registered its extrapolations and let public data rule on them. Seven blind adjudications (P2–P7, plus one archived framework self-correction) were executed on six independent published single-cell datasets, each with decision lines frozen before unblinding, scripts with fixed seeds, and verdict cards archived verbatim, hits and falsifications alike. Two adjudications (P6, P7) asked the framework's sharpest question directly, fold-change versus absolute-scale attribution, on sequential-stimulus data. The ledger is part of the result.

Table 2 collects the ledger at adjudication level (Fig. 5; Extended Data Fig. 9).

**Table 2 | The blind-adjudication ledger (adjudication level).** Seven adjudications on six independent published single-cell datasets; decision lines frozen before unblinding, fixed seeds, verdict cards archived verbatim, hits and falsifications alike. Aggregate record: two clean hits, five intermediates, six falsifications, one dataset-level limitation. Clause-level unit values, frozen lines and scoring scripts: SI S7; grid view: Extended Data Fig. 9.

| Adjudication | Dataset | Clause (frozen before unblinding) | Outcome | Verdict |
|---|---|---|---|---|
| P2 | Msn2 trajectories, nine promoters $\times$ four doses^17^ | Event-time channel is insensitive to cell-size stratification | Line reached in 6/7 promoters ($k = 3$; 7/7 at $k = 2$, 6/7 at $k = 4$); $|\Delta| \leq 0.01$ | Hit |
| P2 | same | Size heterogeneity is priced into the molecule-number channel | Stratification penalty $\Delta \in [-0.049, +0.029]$; 19/21 bootstrap CIs cross zero | Falsified |
| P3 | NF-$\kappa$B microfluidic spatial gradient^18^ | Self-referenced event time decodes spatial ordinality | AUC 0.907 and 0.819 across distance terciles | Hit |
| P3 | same | Stimulus-duration axis decodes dose | AUC 0.591/0.521 against the 0.60 line | Falsified |
| P3 | same | High-dose ordinality ($30{\leftrightarrow}100$ ng/ml) | AUC 0.518; the paid peak-amplitude statistic still decodes at 0.761 | Falsified |
| P4 | ERK GPCR-KTR dose series, 7,393 + 7,508 cells, two ligands^19^ | Free event-time dose decoding | All eight AUCs below the 0.60 line | Falsified |
| P4 | same | Stratification invariance | $|\Delta| \geq 0.03$ in 5/10 units | Falsified |
| P5 | GPCR$\to$Ca²⁺ single cells, 195 valid^20^ | Dead-zone clause | 2/4 units outside the [0.40, 0.60) band | Falsified |
| P5 | same | Counting-channel decoding | Intermediate; population response rate and mean count rise monotonically over the 100-fold ladder ($0.015 \to 0.918$; $0.02 \to 4.36$ counts) | Intermediate |
| P6 | NF-$\kappa$B sequential stimulation, 18 units^21^ | Fold-change versus absolute attribution, $|\Delta R^2| > 0.15$ | $\Delta R^2 = -0.0014$, CI $[-0.252, 0.153]$ crosses zero | Intermediate (power-limited) |
| P7 | Chemotaxis background ladder, 25/30 units^22^ | Same attribution on the amplitude coordinate | $\Delta R^2 = +0.118$ / $-0.091$ (erratum re-adjudication); both CIs cross zero | Intermediate (power-limited) |



Six falsifications do the bounding work: the cell-side redemption of invariance pricing does not run through the size-gain axis (P2), and the sufficiency of the free statistics is bounded by unsaturated ordinal axes (P3) and by system morphology (P4, P5), while population-level counting can be intact where per-cell decoding is weak (P5). Two adjudications (P6, P7) addressed the framework's sharpest dividing question, fold-change versus absolute-scale attribution, on sequential-stimulus designs, and both returned power-limited intermediates that we report rather than resolve post hoc. Across the full ledger the count stands at two clean hits, five intermediates, six falsifications and one dataset-level limitation; every clause, frozen line and unit value is in SI S7.

One earlier falsification is on the record by design. A previous version of this programme predicted, by mechanistic extrapolation from the static gain grading of the PdPC, that ERK pulse amplitude should grade with dose; source-figure discrimination on public data returned $\alpha \approx 0$ and falsified it. In the full version history, the discrimination theorem's classifications on data have never been wrong; the one failure was a mechanistic extrapolation across readout levels, and the record is kept because the falsifiability of the framework is a fact, not a rhetorical posture. Two adjudications of a different kind are archived alongside: a proposed mechanistic patch for the NF-κB model–experiment gap (switch-ified IκBα degradation) was tested computationally and ruled wrong (parameter-unreachability test), and the saturated-transient hypothesis for the same gap was quantitatively rejected. The framework adjudicates hypotheses about systems, not only classifications of systems (Extended Data Fig. 7).

### The analogue account ruptures at the receptor layer

The third instance is a stress test that failed, the framework's most instructive result. Bacterial chemotaxis is the canonical precision-sensing system^23^, and its best public single-cell dataset (ref. 22: FRET dose–responses across seven background levels, five foreground steps each, 57–268 cells per background level) admits two statistics that must, under the standard MWC-plus-perfect-adaptation model^24–26^, be generated by one parameter set: the amplitude table $R(B, F)$ (population medians) and the per-cell midpoint distribution $K_{1/2}(B)$. They cannot be. Along the full Pareto front of joint fits, no point achieves amplitude $R^2 > 0.9$ together with midpoint error below $0.1$ dex: amplitude-favouring fits require cooperativity $N \approx 12$–$39$ and a non-physical saturation $a_{\max} = 1.4$–$2.7$; midpoint-favouring fits require $N \approx 2$–$4$ and destroy the amplitude table. Eight minimal repair paths (two subpopulations, imperfect adaptation, protocol mismatch, readout nonlinearity, ligand depletion, background-dependent gain, heterogeneous censored populations, asymmetric cooperativity) were excluded one by one under a uniform acceptance criterion, with all failures double-recorded including three of our own pipeline errors caught and repaired mid-course. The rupture is not a population-sampling artefact: rebuilding both statistics from strictly the same estimable cells leaves the Pareto front unmoved. It is a statement about the standard model class itself. (Fig. 6)

What survives the rupture is the point. The measured per-cell midpoints obey a textbook Weber line, $K_{1/2} = 1.17 \cdot (1.95 + B)$ µM, with maximum deviation 0.10 dex across seven backgrounds spanning four orders of magnitude of ligand (0.01–100 µM), a pure ratio law, invariant under any rescaling of absolute sensitivity. The scale-free statistic is exact; the absolute-scale parameterisation is inconsistent. The receptor layer, the place where the analogue description should be strongest, is where the analysis finds it broken, and finds the surviving, exact description to be the scale-free one.

The dataset's mixed-ligand arm then provided something rarer than a retrospective fit: a prospective test. From the rupture analysis we derived blind predictions for two unseen mixed-ligand arms (one MeAsp-foreground arm on a mixed MeAsp + L-Asp background, one L-Asp-foreground arm on a MeAsp-dominated background) froze them with tolerances before touching the corresponding data, and scored them per clause (SI S8.6). Predictions made on the scale-free coordinate (the per-cell midpoint on the total effective ligand axis) hit in every scored reading (V1, V2, and the $\mathrm{V2}_{\mathrm{alt}}$ ambiguity control: deviations +0.014, +0.038 and +0.031 dex against a frozen $\pm 0.15$ dex tolerance; the point prediction for the foreground-axis midpoint read 0.86 µM against measured 1.24–1.34 µM). Predictions made on the absolute-amplitude coordinate missed in both primary arms (maximum deviations 0.105 and 0.126 against a frozen 0.08 tolerance in the response-amplitude clauses; the amplitude-coordinate structural test likewise missed). Blind prediction succeeds exactly on the statistics the framework says should survive, and fails exactly on the coordinate the analysis found broken. For the main line this is the negative-image instance with a prospective edge: not only does the absolute-scale account of sensing fail on its own terms; its failures and its survivors are both predictable in advance, from the coordinate system alone (Extended Data Fig. 10).

### The field's analogue currency fails its own consistency checks

The fourth instance steps one level up, from single systems to the quantitative currency in which an entire field reports signalling. GPCR pharmacology reports ligand action in absolute-scale parameters, affinity (p$K_{\mathrm{i}}$/p$K_{\mathrm{A}}$), potency (p$\mathrm{EC}_{50}$), and derived bias factors ($\Delta\log(\tau/K_A)$), and its operational model (Black–Leff)^27,28^ supplies fit-free algebraic identities among them, so that internal consistency can be checked directly. We checked three anchor studies and the field database (Table 3; Fig. 7; Extended Data Figs. 11, 12).

**Table 3 | The pharmacological meta-analysis at a glance.** Values quoted verbatim from the frozen verdict cards (SI S9); the repair figures refer to the minimal coordinate-aware extension (code 80).

| Object | Finding | Size | Significance and note |
|---|---|---|---|
| D2R (anchor study) | Functionally inferred affinity of the high-affinity partial agonists departs from their binding affinity; the fracture sits precisely at the ligands on which the study's "significant bias" reports load-bear | $\Delta = -3.5$ dex | 21–23$\sigma$ (bifeprunox); the non-equilibrium drift component is quantitatively explained by the authors' own measured binding kinetics (zero free parameters), while the absolute discrepancy survives |
| AT1R (anchor study) | The same ligand's functional affinity differs between the Gq and $\beta$-arrestin arms, a cross-arm fracture with binding affinity sitting exactly between the two readings | $\Delta = -2.05$ dex | 5.0$\sigma$; replicated in a mutant |
| μOR (anchor study) | One assay column (GIRK) systematically fractured across all ligands, located uniquely by two independent identities; the column on which the paper's therapeutic-window conclusion rests | 0.68–0.97 dex | 3.9–8.4$\sigma$; a minimal coordinate-aware extension resolves the fracture outright ($\Delta\log\mathrm{ML} = +128.2$; residual median $0.811 \to 0.107$ dex) |
| Biased Signaling Atlas (meta level) | Archived ligand pairs violate the model's cell-level identity; the database's inclusion choices silently dropped the fractured μOR column, changing the evidence base of its conclusion | 21.5% of 17,987 pairs beyond 0.3 dex; 4.2% beyond 1 dex; 15 of 66 multi-ligand papers (23%) exceed 0.3 dex median deviation | The Atlas's own arithmetic pipeline is exact (3,931/3,931 identity closures) |

We are explicit about what this instance does and does not claim. It does not claim that biased signalling is unreal; it does not adjudicate which side of any individual discrepancy is correct. It claims that the absolute-scale parameters in which the field reports communication strength fail internal consistency checks at a substantial rate, concentrated at exactly the load-bearing positions of quantitative conclusions, while the study's own consistency checks, which are identity-based and therefore scale-free, are exact wherever they apply. Three robustness and repair layers back the claim (SI S9). The method's own false-positive rate was measured on synthetic null batteries before any verdict was issued (code 72); the headline rates survive common-mode error calibration, group-level aggregation and column-level jackknife (code 73); and the analysis extends from diagnosis to repair: a constraint projection maps each fractured record to its nearest identity-consistent point with an attached fracture probability (code 77), and for the μOR column a minimal coordinate-aware extension of the operational model, admitting one additional arm coordinate, resolves the fracture outright (Table 3; code 80). The fractures are therefore structural, not noise, and repairable at a price the model itself names. The currency fails the test; the test's currency does not.

### Carrier translation across three analysed cascades

The four analyses above treat each system at a single level. We next asked what happens to a carrier when signalling is read across a cascade, exploiting two datasets that each span two signalling levels at once (Table 4; SI S17, S17b).

**Table 4 | Carrier translation across three analysed cascades.** Values quoted verbatim from the frozen verdict cards (SI S17, S17b).

| Cascade | Where the absolute scale dies | What survives the relay | Key values |
|---|---|---|---|
| M3R–ACh $\to$ Ca²⁺ (353 single cells, seven doses, five repeats per cell)^20^ | Across-cell gain distribution; per-experiment $\mathrm{EC}_{50}$ scattered over a 32-fold range | Within-cell repeatability and per-cell rank order; at low dose, the digital event of responding at all | $\mathrm{EC}_{50}$ median 0.63 µM, CV 1.16 across 26 experiments; within-cell repeatability 3.4–3.7% vs across-cell CV 28–36%; Spearman $\rho = 0.58$–$0.90$ ($0.85$–$0.90$ above the lowest two doses); responding fraction $1.00 \to 0.24$ over the 100-fold ladder |
| p53 cascade | At the sensing stage itself: the noise amplitude is an exact Fisher null direction | Dose reaches the nucleus exclusively as event timing and pulse count | Pulse amplitude varies 0.9% over a 12-fold drive range |
| Bacterial chemotaxis chain | Receptor-cluster layer (the parameter rupture of Fig. 6) | The adaptation relay restores the pre-stimulus baseline; the single paid reference sits exactly at the relay that absorbs the rupture | Baseline restored within 2.2% of response amplitude after 5 s washout, history-independent over a 90-fold dose range (2,217 cells, 77,595 pulse events) |

Three cascades, three molecular hardwares, one rule: absolute amplitude does not survive the first stage, and where a paid reference exists, it is positioned at the relay that has to absorb the fracture. A cell cannot read its own absolute amplitude, because it does not know where it sits in the gain distribution; what remains reportable downstream is self-referenced statistics and, at the low-dose end, the digital event of responding at all.

Relation to prior art. Robust perfect adaptation as a design principle was established for bacterial chemotaxis^24^, fold-change detection as a mechanism of scale-invariant sensing was demonstrated in gene regulation and Wnt signalling^16,72^, and cell-to-cell variability of GPCR dose–response parameters was quantified in the source dataset of the first cascade itself^20^. What the cascade analysis adds is the quantitative connection between these observations and the degeneracy framework: the across-experiment $\mathrm{EC}_{50}$ scatter and the within-cell repeatability are shown to be two coordinates of one gain distribution rather than two separate phenomena; the fracture is shown to be confined to that distribution; and the single paid reference is shown to sit at the relay that must absorb it. The rule stated above is thus a constraint on where references have to sit in any cascade whose upstream layer is uncalibrated, not a redescription of scale invariance.

### One pattern, four systems, one ledger

Fig. 1 summarises the convergence. Four systems, four organisational levels, four independent pre-registered analysis programmes, one stratification: absolute-scale quantities are unidentifiable (PdPC degeneracy group, p53 amplitude channel), internally inconsistent (chemotaxis amplitude table), or unreliable at the field level (GPCR affinity fractures), while scale-free statistics are identifiable ($\kappa$ and $n_{\mathrm{H}}$; the p53 timing channels), exact (the Weber line), prospectively predictive (the mixed-ligand midpoint hits), and information-bearing (96.8% counting-channel efficiency). The digital-limit analysis shows the two constructive instances are the same instance: saturation and topology are two roads to the cancellation of the analogue scale factor from the likelihood. The blind-adjudication ledger (two clean hits, five intermediates, six falsifications and one dataset-level limitation across fourteen pre-registered clauses in seven adjudications) marks the empirical boundary of the claim with the same discipline that produced it. The degeneracy group that makes an assay uncalibratable is, inside the cell, the same mathematical object that makes amplitude untransmittable, and the statistics that survive it are the same: counts, ratios, ordinals, event times.


### Hijacked, not interrupted: the encoding chain under oncogenic mutation

The analyses above characterise healthy encoding. We finally asked what the framework says about its breakdown, testing a chain of linked predictions in cancer cells with pre-registered verdict lines at every link (P11–P16, rows of Table 5; Extended Data Fig. 16; SI S19). The chain asserts that oncogenic lesions do not interrupt signalling but hijack its encoding economics: a mutation that prices itself into the cheap direction should leave measurable traces at each cascade level. Link one is a falsification that sharpened the claim: slow repair was predicted to yield more p53 pulses at fixed dose, and the twelve-cell-line record^77^ reversed it (Spearman $\rho = +0.857$, one-sided $p = 0.0068$, seven oscillatory lines), because repair failure switches the encoder from the digital pulse regime to an analogue sustained regime: a language change rather than a silence, and the counting law survives only inside the oscillatory regime.

Link two reads the switch downstream: under matched total exposure (arm ratio 1.65 from archived TP53 integrals, corrected), 20.4% of induced genes remain more than twofold elevated toward the sustained regime, with arrest and terminal-fate genes (CDKN1A, MDM2, TP53I3, BTG2) amplified up to $+2.6$ bits while repair genes (DDB2, XPC) track the exposure integral ($\le 0.5$ bit): two decoder classes, counting-type and level-type, coexist in the same nucleus^79^. Link three measures the Ras-to-ERK channel in isogenic single-allele lines^78^: ligand-dose mutual information falls from 0.43 bits (KRAS$^{\mathrm{WT}}$) to 0.09 (KRAS$^{\mathrm{Q61R}}$) and 0.01 (BRAF$^{\mathrm{V600E}}$), graded with allele strength, while the across-cell gain distribution is unchanged ($|\Delta\mathrm{CV}| \le 0.014$): compression acts on channel capacity, not dispersion. Link four asks where resistance mutations walk: clinically emergent ABL and EGFR variants preserve catalytic activity (fold ratios 2.94, $p = 0.033$; 2.2, $p = 0.002$), with the registered L718X subset splitting the resistance manifold into activity-preserving and binding-disrupting branches^80^. The outer boundary is empirical: at bulk-tissue resolution the receptor–pERK coupling is absent even in wild-type tumours ($|r| \le 0.18$, four cancer types, ${\sim}1{,}000$ tumours), so the chain does not project onto bulk RPPA and its tissue extrapolation requires single-cell-resolved measurement. Cancer, in the record assembled here, is not a broken telephone line: it is a line speaking a cheaper language, with downstream decoders reading the old price list.


## Discussion

The thesis of this paper is narrow and falsifiable: cellular communication is implemented in statistics invariant under the multiplicative-scale degeneracy group of the intracellular environment, because nothing else is reportable across it. Four analyses support it from four directions, two constructive (what can be identified), two destructive (what fails on its own terms), and a pre-registered adjudication ledger records where its extrapolations failed.

**What this study adds.** The information-theoretic literature has established how much information flows through biochemical circuits and what architectures optimise the flow^29–31^; the degeneracy literature has documented that behaviour underdetermines parameters^32–34^. The analysis addresses the complementary epistemic question, what the flow can reveal about the circuit, and turns "parameters are poorly constrained" into a computable object: exactly which combinations are constrained, and how well, with reporting protocols attached. To be explicit: nondimensionalisation, degeneracy groups, Fisher stratification, scan-based one-sided bounds and orthogonal experimental design each have decades of precedent; we claim no methodological novelty, but novelty of application and output: exact degeneracy groups with reporting protocols for two circuits, a formal exclusion from observable signatures to topology classes (upgrading "tried and failed" to "impossible in principle"), a precisely stated falsifiable digital-limit conjecture, and two field-level inconsistency maps whose failure modes were not previously quantified.

**Consistency across systems.** The signature–topology correspondence is consistent with established physiology beyond the p53 case. In Xenopus oocyte maturation the MAPK cascade is a near-digital static switch (Hill coefficient $\gtrsim 35$, positive feedback^35,36^; $\kappa \lesssim 0.015$ by the $n_{\mathrm{H}} = 1 + 1/(2\kappa)$ relation); in mammalian single cells the same ERK pathway produces frequency-modulated fixed-amplitude pulses^37^, the population-graded response arising from the time integral of pulse frequency. The NF-$\kappa$B core circuit lacks a positive loop^38,39^, so by Lemma 1 it cannot occupy the digital-pure regime, and the data agree at finer resolution^40^: inter-peak intervals do not encode (~75–95 min, constant), first-peak amplitude grades only weakly ($\alpha = 0.151$), and oscillation count increases with dose ($\approx +0.6$ per decade). At the tissue level, radioresistant tissues show oscillatory p53 (digital, measurement mode) whereas radiosensitive tissues show sustained p53 (switch, decision mode): the encoding architecture matches physiology.

**The gauge principle.** Because observables depend only on dimensionless combinations, the functionally complete description of a circuit lives in dimensionless coordinates; absolute molecular scales are a redundant gauge. Three practices follow. Experimentally, perturbations and cross-system comparisons are meaningful only when stated on dimensionless groups: comparing absolute rate constants across cell types conflates gauge with physics. Evolutionarily, mutation walks in absolute coordinates while selection reads dimensionless ones; the degeneracy group is the neutral-drift space available to a circuit, which formalises how ion-channel degeneracy supports resilience^33,34^. Functionally, the principle rationalises why biological sensing is so often ratio-based (fold-change detection^16^, morphogen-gradient scaling by tissue size, the kinase/phosphatase activity ratio): a cell has no external calibration standard, so evolution can only have built its sensors on dimensionless quantities. Counting is the limiting case in which the observable is dimensionless by construction.

**The free-first principle.** The record assembled here supports a stronger, cost-effectiveness reading: wherever a calibration-free statistic suffices for the function, it is the one used, and payment for an actively maintained reference occurs exactly where function must reference absolute scale, nowhere else (SI S1.6b). Table 6 collects the per-system carriers, and the pattern holds in every system examined without exception: p53 counts dose while amplitude is dose-independent and functionally silent^9,11^; Msn2 keeps frequency; ERK digitises; chemotaxis reads fold change over two orders of magnitude of absolute ligand level while paying for exactly one reference^62^; Wnt/$\beta$-catenin holds a quantised plateau invariant to ±5% across 0–50 mM LiCl^72^; and the embryo reads relative position $x/L$ off an absolute anchor supplied free by physics ($\alpha = -0.18$, 95% CI $[-0.48, 0.11]$ on 581 profiles, perfect scaling excluded at ${\sim}7\sigma$)^71,73,74^. A single counterexample, a paid reference maintained at a position whose function needs only relative quantities, falsifies the principle. The markers themselves are validated by positive control: in yeast osmotic homeostasis, where function unambiguously references an absolute quantity, all three class-III markers are observed, and causal inhibition of the maintenance mechanism abolishes absolute-scale recovery (SI S1.6c, Table 6).

**Reverse predictions.** The cellular line of this programme has a machine-side sibling: artificial perception networks evolved under controlled perturbation regimes, analysed under the same economics (calibration-freeness is purchased, not assumed). That line produced three verdict-grade laws, each stating that selection prices invariance, and each exporting a falsifiable cellular prediction (Table 7). All three are frozen here in advance of testing; failure adjudicates equally, and if cellular circuits do not obey the cross-domain economics of invariance, the digital-limit principle is a cellular special case.

**An experimental protocol for the unused signature.** The hysteresis signature of the discrimination scheme has never been measured in the p53 literature: all published dose series are one-directional scans. We propose the direct test: scan DNA-damage dose upward and downward in single cells and compare onset and extinction points. The framework predicts zero hysteresis if the core generator is purely excitable; any nonzero width forces reclassification into the subcritical column and downgrades "excluding pure negative feedback" to "excluding supercritical NF". The same experiment adjudicates whether Lemma 1's subcritical escape hatch is realised by nature.

**Boundaries and limitations.** Seven considerations bound the present claims (Table 8); each is registered with the same discipline as the confirmations, and each names the measurement that would lift it.

**The prediction ledger.** The framework's claims were frozen before the corresponding data were inspected, and Table 5 is the complete win-loss record of the prediction campaign reported here and in SI S15–S19. Every entry carries its frozen form, its pre-registered verdict line, and its outcome; falsifications are listed with the same prominence as confirmations, because each falsification contracted the claim to what survived.

**Table 5 | The prediction ledger.** The complete win-loss record of the pre-registered prediction campaign (main text and SI S15–S19); falsifications are listed with the same prominence as confirmations.

| Frozen prediction | Verdict line (pre-registered) | Outcome | Status |
|---|---|---|---|
| Dose is encoded in p53 pulse count, not amplitude | Timing share of Fisher information dominates; amplitude dose-invariant | 94.2% timing vs 5.5% amplitude at representative level (calibration-bounded, E10); amplitude dose-invariant in source data^9^ | Confirmed, calibrated wording |
| Bcd length constant is an absolute anchor under natural length variation (P9-1) | Exponent $\alpha \le 0.7$ with CI upper bound below $1$; perfect scaling requires point $\ge 0.9$ | $\alpha = -0.18$, CI $[-0.48, 0.11]$, 581 embryos; perfect scaling excluded at ${\sim}7\sigma$ | Confirmed |
| Gap-gene readout tracks relative position, supplied free by nuclear counting (S15) | Scaling coefficient $S \approx 1$ mid-embryo; dose-monotone deviation | Mid-embryo median $S = 1.175$ at $1\times$ with 53–59% of boundaries CI-compatible with $S = 1$; CI-weighted mean $S$ falls monotonically with bcd dose ($1.07, 0.96, 0.84$) | Confirmed |
| Wnt signalling reads fold change with a quantised breakdown point (S16) | Plateau within $\pm 20\%$; breakdown localised | Plateau 0.96–1.04 over 0–50 mM LiCl; breakdown onset 50–58 mM | Confirmed |
| Chemotaxis mixed-ligand arm (prospective, blind) | Predicted shifts within frozen band before data release | $+0.014, +0.038, +0.031$ dex, all inside the frozen band | Confirmed |
| Gain-distribution confinement generalises across pulse systems (P10) | Confinement ratio $\ge 3$ in ERK and NF-$\kappa$B | Ratios 0.59–2.17; line not met in either system | Falsified; claim contracted to the M3R-to-Ca cascade |
| Amplitude stratification (Msn2); duration-axis decoding (NF-$\kappa$B); free event-time dose decoding (ERK-KTR) | Per-arm frozen statistical lines (SI S7) | All three lines failed on their respective data | Falsified; removed from the claim set |
| Size-gain operationalisation of cheapest-sufficient-invariance (P2) | Frozen gain-size relation | Relation absent | Falsified; clause contracted |
| Zero hysteresis in p53 dose response (excitable core) | Up-down dose scan: onset equals extinction | Experiment never performed in the literature; protocol frozen here | Awaiting test |
| Slow DNA repair yields more p53 pulses at fixed dose (P11-1) | Negative correlation, $|\rho| \ge 0.7$ | $\rho = +0.857$, $p = 0.0068$: direction reversed; repair failure switches the encoder from pulses to a sustained regime | Falsified; counting law contracted to the oscillatory regime |
| Oncogenic Ras/Raf mutations compress ligand-to-ERK channel capacity (P13) | $\Delta I < -0.2$ bits, paired bootstrap CI excludes zero | G12D $-0.228$, Q61R $-0.335$, BRAF-V600E $-0.413$ bits; G12V $-0.050$ (CI includes zero); graded with allele strength | Confirmed for strong alleles; contracted to an allele-strength-graded claim |
| Downstream decoders misread the p53 regime switch (P14) | output difference above $2\times$ at matched exposure | 20.4% of induced genes differ by more than $2\times$ after AUC correction; arrest genes amplified by the sustained regime, repair genes shape-insensitive | Confirmed with registered AUC correction |
| Oncogenic Ras/Raf alleles widen the across-cell ERK gain distribution (P12a) | CV increase beyond the frozen band | Baseline CV 0.339/0.353/0.353 across the isogenic lines; $|\Delta\mathrm{CV}| \le 0.014$ | Falsified; informative null, compression acts on capacity not dispersion |
| Resistance mutations preserve catalytic activity (P12b) | Median ratio above $2$, one-sided $p < 0.05$ | ABL: 2.94, $p = 0.033$; EGFR DMS: 2.2, $p = 0.002$; L718X subset reverses the rule | Confirmed; manifold split into two branches registered |
| Tissue-level pERK decouples from receptor abundance in mutant tumours (P15) | Mutant slope ratio below $0.8$ with interaction $p < 0.05$ | Coupling absent even in wild-type tumours across four cancer types | Falsified at bulk resolution; extrapolation boundary registered |
| KRAS-G13D channel capacity intermediate between G12V and G12D (P16) | $\Delta I \in [-0.25, -0.05]$ bits, ordering $\mathrm{G12V} > \mathrm{G13D} > \mathrm{G12D}$ | No comparable public dataset as of 2026-09-28 | Awaiting test |

Two literature-level consistencies were added post hoc and are labelled as such (they are not frozen predictions): the cross-species constancy of $\lambda/L$ within 2% under a fixed nuclear count^75^, and the twelve-cell-line invariance of the p53 period near 5 h with all variation confined to envelope properties^77^. Both are consistent with the framework and neither counts toward its score.

**Spin-off discoveries.** The analyses paid for themselves in findings that are independent of the main theorem; Table 9 lists the five sharpest, and SI S18 gives the complete ledger of eleven. Each is stated with its evidence anchor, and each is falsifiable on its own terms.

**Outlook.** The onset-bifurcation classification of encoding has a thirty-year-old precedent in neuronal excitability (Type I/SNIC versus Type II/Hopf)^41^, with the analogue signature living in the frequency rather than the amplitude observable, suggesting a general duality, encoding signature as a function of (bifurcation type $\times$ observable), unifying molecular and electrical pulse codes. A first analysis at textbook Morris–Lecar working points is consistent: at Type-I onset the frequency channel grades as a square-root law while the amplitude channel saturates, and under synaptic-gain drift spanning $0.3\times$–$3\times$ an event-time readout degrades by 0.000% while an amplitude readout drifts by more than 70% on average: gain scale dies, event timing survives; the amplitude leg of the duality fails at subcritical Type-II onset, bounding the duality's domain by bifurcation criticality and the chosen observable (Extended Data Fig. 8). If the digital-limit conjecture holds, it is a design principle: synthetic circuits that must transmit dose under unresolvable parameter uncertainty should be built to count.

The same law has now been found on both sides of the measurement boundary: assay platforms count when they cannot calibrate, and cells, the original uncalibratable instruments, do the same. A companion study extends the analysis to intercellular and neural communication and to artificial networks, where the price of paid amplitude can be measured in decodable bits. The cellular layer, established here, is the anchor: four systems, one currency.

## Methods

### Analysis procedure

Every system passed through the same five-step analysis (dimensionless master equation; exact degeneracy group; Fisher stratification at experimentally realistic noise; one-sided reporting protocols for plateau directions; orthogonal-resolution identification). All verdicts were pre-registered before scoring; failures, including our own pipeline errors, were double-recorded alongside successes and are archived verbatim. Numerical continuation^42^ used custom continuation validated against analytic normal forms (Hopf points to $\pm 10^{-4}$ in the control parameter; $A^2$-linearity $R^2$ reported per branch). Every number quoted in the main text is transcribed from a frozen verdict card, not from session memory (Extended Data Fig. 13).

### PdPC analysis

The symmetric and asymmetric Goldbeter–Koshland steady states were solved by brentq ($\mathrm{xtol} = 10^{-14}$) and nondimensionalised to $(\xi, \kappa)$ and $(\xi, \kappa_1, \kappa_2)$; degeneracy-group actions (rate rescaling; ($K_m$, $S_T$) co-scaling; mixed) were verified pointwise on a 200-point logarithmic titration grid, with invariance holding to machine precision (max $|\Delta u| = 0$). The Fisher information matrix was computed from analytic response Jacobians at observation noise $\sigma = 0.03$; deep-saturation confidence widths for $\ln \kappa$ by Monte Carlo refitting ($n = 400$ synthetic curves per $\kappa$). Mutual-information optima used log-normal input ensembles ($\sigma_\xi = 0.3$ narrow, 1.5 wide), Gaussian output noise $\sigma = 0.03$, adaptive grids converging to $\pm 0.02$ bits. The positive-feedback switch (fixed-point equation $u = F(\xi(1+g\cdot u))$, critical gain $g^{*} = 0.443$ at $\kappa = 0.05$) was analysed identically, with all four group actions verified in both the one-root and three-root regions.

### p53 analysis

Topology classification used a canonical three-variable NF model (analytic equilibrium branch; 6000-point stability scan with the $3 \times 3$ analytic Jacobian; unique destabilisation $s_H = 16.4449$ via complex pair $\pm 0.924i$; limit-cycle amplitudes by solve_ivp at $\mathrm{rtol} = 10^{-9}$, $A^2 \propto (s - s_H)$ at $R^2 = 0.9926$) and a time-scale-separated excitable representative (FitzHugh–Nagumo form; 6.5% amplitude variation over a 4.6-fold drive range), with a subcritical normal-form control (quintic-stabilised; measured hysteresis width 0.275 versus theory 0.25) and a Morris–Lecar Type-II control ($\alpha = -0.258$). The counting law $N \approx (\tau_r/T) \cdot \ln(D_0/D_c)$ was fitted over the sustained-counting window with period self-consistency. The channel census propagated the Fisher information of six pulse statistics ($t_1, T, A, N, w, \mathrm{IPI}$) onto seven dynamical parameters, with three numerical configurations ($M = 400$ central, $M = 4000$ forward, $M = 4000$ central differences) confirming the null-direction ratio $\lambda_7/\lambda_1 \leq 4.7 \times 10^{-19}$ and eigenvector $\sigma$-weight 0.996–0.999. The stochastic population model (200 cells per dose; Euler–Maruyama $dt = 0.001$; dynamical noise $\sigma = 0.10$; lognormal gain dispersion and onset jitter; time unit calibrated by median $\mathrm{IPI} = 5.5$ h) generated the gain-dispersion kill and both dispersion signatures. Experimental placement used source-digitised dose series (Lahav et al.^4^ $\gamma$-irradiation 0.1–10 Gy; Mönke et al.^10^ NCS 25–400 $\mathrm{ng\,ml^{-1}}$; Batchelor et al.^9^ Fig. 1G/H colour-segmentation digitisation with visual axis calibration, bootstrap $\times 20000$: $\alpha_{\mathrm{NCS}} = -0.001 \pm 0.113$, $\alpha_{\mathrm{UV}} = +0.981 \pm 0.114$; duration exponents $-0.21$ / $+0.79$; $\tau_r \approx 7.7$ h).

### Digital-limit analysis

The conditional theorem (two-circuit constructive proof of likelihood cancellation onto partition-counting statistics) and its unconditional generalisation (stratified classification over degeneracy groups; universal functional theorem $T'$; noise-robustness bound Theorem N) are stated and proved in SI S5, with the non-tautology analysis (within-bin amplitude moments; the coefficient of variation; continuous ratio statistics; the boundary against representational measurement theory) and the explicit uniqueness domain (single channel, single readout, memoryless, passive sensing).

### Framework formalism and self-checks

The framework of Results section 1 (layer model, degeneracy group, survival set, three-leg criterion, information-loss identity, three reference classes, predictions P-1/P-2/P-3) is stated with full definitions in SI S1. Its development lineage is part of the recorded trail: seven internal bugs and sixteen external-review items across versions v2.0–v4.2 were recorded and repaired, including the two substantive corrections of the final version (the survival-set definition moved from raw output distributions to transmitted-alphabet distributions, without which the set is generically empty; and the layer-boundary rule restored). Framework self-checks were run as simulations before any system analysis: a p53 cascade three-leg verification (survival set, symbolic reachability, and bandwidth legs checked jointly, with the information-loss identity verified numerically) and a chemotaxis receptor-layer self-check; both are archived with their scripts.

### Blind adjudications P2–P7

Each adjudication had a pre-registration document, a frozen decision table executed verbatim by the scoring script, fixed seeds (SEED = 20260814 for P2/P3, 20260815 for P4–P7), bootstrap $\times 2000$ confidence intervals, and an archived verdict card. P2 (Msn2; Hansen & Zechner^17^, Zenodo 10.5281/zenodo.2755026): nine promoters $\times$ four doses; size-tercile stratification penalty of the molecule-number channel (21 units) versus the calibration-free event-time channel (7 promoters; k = 2/4 robustness arms). P3 (NF-κB spatial gradient; Son et al.^18^, Zenodo 10.5281/zenodo.6858118): per-cell self-referenced free statistics only (first-crossing event time above the cell's own pre-stimulus $\mu_0 + 3\sigma_0$; above-threshold duration); spatial ordinality (near/mid/far terciles), dose ordinality ($10{\leftrightarrow}30$, $30{\leftrightarrow}100$ $\mathrm{ng\,ml^{-1}}$), and stimulus-duration ($15{\leftrightarrow}30$, $30{\leftrightarrow}60$ min) axes; decision lines derived from the null noise floor before unblinding. P4 (ERK re-adjudication; Chavez-Abiega et al.^19^ GPCR-KTR; histamine 7,393 cells and UK14304 7,508 cells, six doses each): free event-time dose decoding ($\tau(3)$ AUC), saturation-region ordering, and stratification invariance. P5 (boundary predictor; Keshelava et al.^20^ GPCR→Ca²⁺; 195 valid cells from 11/27 experiments under frozen inclusion rules, power-reduction qualifier recorded): dead-zone clause, live-zone counting clause, and a paid-peak control arm; an erratum (47b) re-adjudication under the corrected inclusion rule is archived alongside. P6 (fold-change attribution; Wang et al.^21^ NF-κB sequential stimulation; 18 valid units): $\Delta R^2$ between fold-change and absolute-dose readouts, pre-registered line $\pm 0.15$, power-limitation qualifier recorded on the card. P7 (fold-change/Weber attribution; Moore et al.^22^ background ladder): $\Delta R^2$ between fold-ratio and total-concentration coordinates on the amplitude table, line $\pm 0.15$, with 25 condition units in the first pass (code 49) and 30 in the erratum re-adjudication under the corrected unit definition (code 49b); both passes archived. The mixed-ligand prospective arm (code 70; pre-registration P1/P2) froze point predictions and tolerances for two unseen mixed-ligand conditions before the corresponding data were touched, then scored per clause: response-amplitude clauses (tolerance 0.08 in da), amplitude-coordinate structural clause ($R^2$ ordering), and per-cell midpoint clauses on the total effective axis (tolerance $\pm 0.15$ dex). The ERK self-correction used source-figure digitisation of Albeck et al.^37^ (Fig. 3A, MCF-10A; Fig. S4A, 184A1) with pulse counting by peak detection (prominence > 0.15 EKAR). Key adjudications (P2, P3, P4) were independently re-executed in a clean cloud environment from the frozen scripts.

### Chemotaxis analysis

Data source: the Moore et al.^22^ Dryad archive (doi:10.5061/dryad.nvx0k6dzz, CC0). Amplitude tables and per-cell $K_{1/2}(B)$ (log-axis interpolation, right-censoring at $\mathrm{da}_{\max} < 0.5$) were jointly fitted under the MWC-plus-perfect-adaptation class along a Pareto weight sweep; eight repair paths were tested under the uniform acceptance criterion ($R^2 > 0.9$, $K_{1/2}$ error below $0.1$ dex, $a_{\max} \approx 1$); the same-source verdict rebuilt both statistics from strictly the same estimable cells; three internal pipeline errors were caught, repaired and double-recorded. The fold-change/Weber attribution adjudication (P7; codes 49/49b) and the mixed-ligand prospective arm (code 70; pre-registration P1/P2) are described under Blind adjudications; a further prospective cohort (P8) is registered and awaits the necessary dataset.

### GPCR meta-analysis

Black–Leff operational-model identities were checked cell-wise on D2R (Klein Herenbrink et al.^43^), AT1R (Wingler et al.^44^), μOR (Gillis et al.^45^, with the Stahl & Bohn reanalysis^46^ recorded), and the full Biased Signaling Atlas (9,041 rows, 214 papers; 17,987 ligand pairs, full-library identity pipeline code 71). Non-equilibrium occupancy corrections used the authors' own measured binding kinetics with zero free parameters. Significance is reported in dex with full error propagation and tracer-robustness checks. The method's false-positive rate was calibrated on synthetic null batteries before any verdict (code 72); robustness used common-mode error calibration, group-level aggregation and column-level jackknife (code 73); repair used constraint projection with fracture probabilities (code 77) and, for μOR, a coordinate-aware minimal extension of the operational model scored by model comparison (code 80). An admission-ticket rule rejected one candidate system (PTH1R) before scoring; the rejection is double-recorded as a boundary teaching case. The Atlas-wide inconsistency map is provided as an interactive supplement (`06_supplement_interactive/app/index.html`).

### Cancer-chain tests (P11–P16)

Every link of the chain carried its own pre-registration card with frozen verdict lines before data access (SI S19). P11 digitised twelve-cell-line p53 statistics^77,82^ (digitisation error bands registered per figure). P13 re-analysed the isogenic MEF single-cell EKAR3 archive^78^ (Source Data Fig. 2, 8 lines, 6 ligands, 4 doses including ligand-free): the per-cell output was locked as the 0–30 min post-stimulus peak minus the pre-stimulus median; mutual information used equiprobable 8-bin discretisation per line with Miller–Madow correction and 2,000 cell-level paired bootstraps, and was re-verified with five bin counts and a continuous KSG estimator (300 bootstraps). P14 re-analysed the matched-exposure oscillatory-versus-rising p53 archive^79^ (RNA-seq TPM and TMT proteomics from the authors' public repository); the arm imbalance in total p53 exposure (ratio 1.65 by our own TP53 integral check) was corrected on a registered log-linear assumption, and genes with infinite fold ratios (zero TPM in one arm) contribute symmetrically below 1% per side. P12b used the ABL catalytic-efficiency table and the EGFR saturation-mutagenesis Ba/F3 screens^80^ (GEO GSE305057 codon counts); activity was proxied by DMSO-arm variant frequency normalised to synonymous variants, and a self-built enrichment classification was discarded after registration when low-count ratio inflation was detected (amendment A5: author-validated variant lists used instead). P15 used cBioPortal TCGA Pan-Cancer Atlas RPPA (phospho-ERK = mean of MAPK1/3 PT202/Y204; total EGFR) and MC3 mutation calls across COADREAD, LUAD, SKCM and PAAD^81^, with the registered asymmetry that a bulk-level null does not overturn the single-cell result.

### Analysis discipline

Three disciplines made the analysis binding rather than rhetorical. First, pre-registration: decision lines, acceptance criteria and exclusion rules for every adjudication were frozen in writing before the corresponding data were unblinded; scoring scripts execute the frozen tables verbatim. Second, frozen verdict cards: every conclusion in this paper lives on a timestamped card (fixed seeds; bootstrap counts; SHA-256 fingerprints in the project ledger) and is quoted from the card, never from memory; frozen verdicts never move, and corrections take the form of appended errata (e.g., adjudication 47b) rather than edits. Third, double-recording: failures are archived with the same completeness as successes: falsified clauses (six in the present ledger), rejected hypotheses (two in the NF-κB line), and our own pipeline errors (three in the chemotaxis line, one in P5), each with cause and repair. Key adjudications were re-executed independently in a clean cloud environment from the frozen scripts. Version history functions as the record trail: the one mechanistic-extrapolation failure of the programme (the ERK amplitude prediction) is preserved in the record, with its cause and its consequence for the claim's domain.

### Data and code availability

No new experimental data were generated. Public datasets re-analysed: Batchelor et al.^9^ Fig. 1G/H (source-figure digitisation); Lahav et al.^4^ and Mönke et al.^10^ (reported statistics); Tay et al.^40^ (reported statistics and Supp. Fig. 5); Albeck et al.^37^ Figs. 3A/S4A (source-figure digitisation); Hansen & Zechner^17^ (Zenodo 10.5281/zenodo.2755026); Son et al.^18^ (Zenodo 10.5281/zenodo.6858118); Wang et al.^21^ (NF-κB sequential stimulation); Chavez-Abiega et al.^19^ (GPCR-KTR dose series); Keshelava et al.^20^ (Nat. Commun. 9:876); Moore et al.^22^ (Dryad doi:10.5061/dryad.nvx0k6dzz); Klein Herenbrink et al.^43^; Wingler et al.^44^; Gillis et al.^45^; Biased Signaling Atlas (9,041 rows); Morton de Lachapelle & Bergmann^71^ (supplementary Datasets S1/S2, parsed verbatim); Goentoro & Kirschner^72^ (figure-level digitisation, Fig. 4F); Liu et al.^73^ and Petkova et al. per-embryo profiles via the scale-invariance repository of Nikolić et al.^74^ (github.com/mnikolic7/Scale-invariance-in-early-embryonic-development); Gillies et al.^78^ (PMC7569415, Source Data Figs. 2 and 6); Jimenez-Asins et al.^79^ (github.com/albajimenezasins/Proteomics_MSB_2022); Wang et al.^80^ (GEO GSE305057 codon counts); Stewart-Ornstein & Lahav^77^ and Finzel^82^ (figure digitisation with registered error bands); TCGA Pan-Cancer Atlas RPPA and MC3 mutation calls via cBioPortal^81^. All analysis scripts, pre-registration documents, frozen verdict cards, errata and adjudication logs are available at github.com/lihuanlz/cell-communication-audit (archived at Zenodo, doi:10.5281/zenodo.23008824); every figure panel regenerates from a named script against a named archived table. The Atlas-wide inconsistency map ships as an interactive HTML supplement (`06_supplement_interactive/app/index.html`).

## References

1. Vogelstein, B. & Kinzler, K. W. Digital PCR. *Proc. Natl Acad. Sci. USA* **96**, 9236–9241 (1999).

2. Guo, L. et al. Scale degeneracy and its resolution: the route to calibration-free absolute quantification. Preprint at ChemRxiv https://doi.org/10.26434/chemrxiv-2024-19rj6/v10 (2026).

3. Goldbeter, A. & Koshland, D. E. An amplified sensitivity arising from covalent modification in biological systems. *Proc. Natl Acad. Sci. USA* **78**, 6840–6844 (1981).

4. Lahav, G. et al. Dynamics of the p53–Mdm2 feedback loop in individual cells. *Nat. Genet.* **36**, 147–150 (2004).

5. Geva-Zatorsky, N. et al. Oscillations and variability in the p53 system. *Mol. Syst. Biol.* **2**, 2006.0033 (2006).

6. Batchelor, E., Mock, C. S., Bhan, I., Loewer, A. & Lahav, G. Recurrent initiation: a mechanism for triggering p53 pulses in response to DNA damage. *Mol. Cell* **30**, 277–289 (2008).

7. Loewer, A., Karanam, K., Mock, C. & Lahav, G. The p53 response in single cells is linearly correlated to the number of DNA breaks without a distinct threshold. *BMC Biol.* **11**, 114 (2013).

8. Soulé, C. Graphic requirements for multistationarity. *ComPlexUs* **1**, 123–133 (2003).

9. Batchelor, E., Loewer, A., Mock, C. & Lahav, G. Stimulus-dependent dynamics of p53 in single cells. *Mol. Syst. Biol.* **7**, 488 (2011).

10. Mönke, G. et al. Excitability in the p53 network mediates robust signaling with tunable activation thresholds in single cells. *Sci. Rep.* **7**, 46571 (2017).

11. Purvis, J. E. et al. p53 dynamics control cell fate. *Science* **336**, 1440–1444 (2012).

12. Purvis, J. E. & Lahav, G. Encoding and decoding cellular information through signaling dynamics. *Cell* **152**, 945–956 (2013).

13. Monk, N. A. M. Oscillatory expression of Hes1, p53, and NF-κB driven by transcriptional time delays. *Curr. Biol.* **13**, 1409–1413 (2003).

14. Wang, C., Liu, H. & Zhou, J. Contribution of time delays to p53 oscillation in DNA damage response. *IET Syst. Biol.* **13**, 180–185 (2019).

15. Belgacem, I. Sustained limit cycles in the logistic two-gene genetic oscillator: a delay-driven Hopf bifurcation. Preprint at arXiv https://arxiv.org/abs/2605.23722 (2026).

16. Goentoro, L., Shoval, O., Kirschner, M. W. & Alon, U. The incoherent feedforward loop can provide fold-change detection in gene regulation. *Mol. Cell* **36**, 894–899 (2009).

17. Hansen, A. S. & Zechner, C. Promoters adopt distinct dynamic manifestations depending on transcription factor context. *Mol. Syst. Biol.* **17**, e9821 (2021).

18. Son, M. et al. Spatiotemporal NF-κB dynamics encodes the position, amplitude, and duration of local immune inputs. *Sci. Adv.* **8**, eabn6240 (2022).

19. Chavez-Abiega, S., Grönloh, M. L. B., Gadella, T. W. J., Bruggeman, F. J. & Goedhart, J. Single-cell imaging of ERK and Akt activation dynamics and heterogeneity induced by G-protein-coupled receptors. *J. Cell Sci.* **135**, jcs259685 (2022).

20. Keshelava, A. et al. High capacity in G protein-coupled receptor signaling. *Nat. Commun.* **9**, 876 (2018).

21. Wang, A. G., Son, M., Kenna, E., Thom, N. & Tay, S. NF-κB memory coordinates transcriptional responses to dynamic inflammatory stimuli. *Cell Rep.* **40**, 111159 (2022).

22. Moore, J. P., Kamino, K., Kottou, R., Shimizu, T. S. & Emonet, T. Signal integration and adaptive sensory diversity tuning in Escherichia coli chemotaxis. *Cell Syst.* **15**, 628–638 (2024).

23. Micali, G. & Endres, R. G. Maximal information transmission is compatible with ultrasensitive biological pathways. *Sci. Rep.* **9**, 16898 (2019).

24. Barkai, N. & Leibler, S. Robustness in simple biochemical networks. *Nature* **387**, 913–917 (1997).

25. Yi, T.-M., Huang, Y., Simon, M. I. & Doyle, J. Robust perfect adaptation in bacterial chemotaxis through integral feedback control. *Proc. Natl Acad. Sci. USA* **97**, 4649–4653 (2000).

26. Ma, W., Trusina, A., El-Samad, H., Lim, W. A. & Tang, C. Defining network topologies that can achieve biochemical adaptation. *Cell* **138**, 760–773 (2009).

27. Black, J. W. & Leff, P. Operational models of pharmacological agonism. *Proc. R. Soc. Lond. B* **220**, 141–162 (1983).

28. Kenakin, T. Biased receptor signaling in drug discovery. *Pharmacol. Rev.* **71**, 267–315 (2019).

29. Tkačik, G. & ten Wolde, P. R. Information processing in biochemical networks. *Annu. Rev. Biophys.* **54**, 249–274 (2025).

30. Cheong, R., Rhee, A., Wang, C. J., Nemenman, I. & Levchenko, A. Information transduction capacity of noisy biochemical signaling networks. *Science* **334**, 354–358 (2011).

31. Selimkhanov, J. et al. Accurate information transmission through dynamic biochemical signaling networks. *Science* **346**, 1370–1373 (2014).

32. Gutenkunst, R. N. et al. Universally sloppy parameter sensitivities in systems biology models. *PLoS Comput. Biol.* **3**, e189 (2007).

33. Prinz, A. A., Bucher, D. & Marder, E. Similar network activity from disparate circuit parameters. *Nat. Neurosci.* **7**, 1345–1352 (2004).

34. Goaillard, J.-M. & Marder, E. Ion channel degeneracy, variability, and covariation in neuron and circuit resilience. *Annu. Rev. Neurosci.* **44**, 335–357 (2021).

35. Ferrell, J. E. & Machleder, E. M. The biochemical basis of an all-or-none cell fate switch in Xenopus oocytes. *Science* **280**, 895–898 (1998).

36. Xiong, W. & Ferrell, J. E. A positive-feedback-based bistable ‘memory module’ that governs a cell fate decision. *Nature* **426**, 460–465 (2003).

37. Albeck, J. G., Mills, G. B. & Brugge, J. S. Frequency-modulated pulses of ERK activity transmit quantitative proliferation signals. *Mol. Cell* **49**, 249–261 (2013).

38. Hoffmann, A., Levchenko, A., Scott, M. L. & Baltimore, D. The IκB–NF-κB signaling module: temporal control and selective gene activation. *Science* **298**, 1241–1245 (2002).

39. Nelson, D. E. et al. Oscillations in NF-κB signaling control the dynamics of gene expression. *Science* **306**, 704–708 (2004).

40. Tay, S. et al. Single-cell NF-κB dynamics reveal digital activation and analogue information processing. *Nature* **466**, 267–271 (2010).

41. Rinzel, J. & Ermentrout, G. B. Analysis of neural excitability and oscillations. In *Methods in Neuronal Modeling* (eds Koch, C. & Segev, I.) 135–169 (MIT Press, 1989).

42. Kuznetsov, Y. A. *Elements of Applied Bifurcation Theory* 3rd edn (Springer, 2004).

43. Klein Herenbrink, C. et al. The role of kinetic context in apparent biased agonism at GPCRs. *Nat. Commun.* **7**, 10842 (2016).

44. Wingler, L. M. et al. Angiotensin and biased analogs induce structurally distinct active conformations within a GPCR. *Science* **367**, 888–892 (2020).

45. Gillis, A. et al. Low intrinsic efficacy for G protein activation can explain the improved side effect profiles of new opioid agonists. *Sci. Signal.* **13**, eaaz3140 (2020).

46. Stahl, E. L. & Bohn, L. M. Low intrinsic efficacy alone cannot explain the improved side effect profiles of new opioid agonists. *Biochemistry* **61**, 1923–1935 (2022).

47. Shoval, O., Alon, U. & Sontag, E. D. Symmetry invariance for adapting biological systems. *SIAM J. Appl. Dyn. Syst.* **10**, 857–886 (2011).

48. Bar-Or, R. L. et al. Generation of oscillations by the p53–Mdm2 feedback loop: a theoretical and experimental study. *Proc. Natl Acad. Sci. USA* **97**, 11250–11255 (2000).

49. Ciliberto, A., Novák, B. & Tyson, J. J. Steady states and oscillations in the p53/Mdm2 network. *Cell Cycle* **4**, 488–493 (2005).

50. Gouzé, J.-L. Positive and negative circuits in dynamical systems. *J. Biol. Syst.* **6**, 11–15 (1998).

51. Bahadur, R. R. Sufficiency and statistical decision functions. *Ann. Math. Stat.* **25**, 423–462 (1954).

52. Krantz, D. H., Luce, R. D., Suppes, P. & Tversky, A. *Foundations of Measurement*, Vol. I (Academic Press, 1971).

53. Chen, H.-L., Doty, D., Reeves, W. & Soloveichik, D. Rate-independent computation in continuous chemical reaction networks. *J. ACM* **70**, Article 22 (2023).

54. Hamerly, R., Bandyopadhyay, S. & Englund, D. Asymptotically fault-tolerant programmable photonics. *Nat. Commun.* **13**, 6831 (2022).

55. Krishna, S., Jensen, M. H. & Sneppen, K. Minimal model of spiky oscillations in NF-κB signaling. *Proc. Natl Acad. Sci. USA* **103**, 10840–10845 (2006).

56. Sachdev, S., Grönloh, M. L. B., Gadella, T. W. J., Bruggeman, F. J. & Goedhart, J. Highly biased agonism for GPCR ligands via nanobody tethering. *Nat. Commun.* **15**, 4687 (2024).

57. Adler, M. & Alon, U. Fold-change detection in biological systems. *Curr. Opin. Syst. Biol.* **8**, 81–89 (2018).

58. Tostevin, F. & ten Wolde, P. R. Mutual information between input and output trajectories of biochemical networks. *Phys. Rev. Lett.* **102**, 218101 (2009).

59. Govern, C. C. & ten Wolde, P. R. Optimal resource allocation in cellular sensing systems. *Proc. Natl Acad. Sci. USA* **111**, 17486–17491 (2014).

60. Bowsher, C. G. & Swain, P. S. Environmental sensing, information transfer, and cellular decision-making. *Curr. Opin. Biotechnol.* **28**, 149–155 (2014).

61. Walczak, A. M., Tkačik, G. & Bialek, W. Optimizing information flow in small genetic networks. II. Feed-forward interactions. *Phys. Rev. E* **81**, 041905 (2010).

62. Lazova, M. D., Ahmed, T., Bellomo, D., Stocker, R. & Shimizu, T. S. Response rescaling in bacterial chemotaxis. *Proc. Natl Acad. Sci. USA* **108**, 13870–13875 (2011).

63. Pantazis, Y. & Katsoulakis, M. A. A relative entropy rate method for path space sensitivity analysis of stationary complex stochastic dynamics. *J. Chem. Phys.* **138**, 054115 (2013).

64. Hao, N. & O'Shea, E. K. Signal-dependent dynamics of transcription factor translocation controls gene expression. *Nat. Struct. Mol. Biol.* **19**, 31–39 (2012).

65. Cai, L., Dalal, C. K. & Elowitz, M. B. Frequency-modulated nuclear localization bursts coordinate gene regulation. *Nature* **455**, 485–490 (2008).

66. Venkatachalapathy, H., Dallon, S., Yang, Z., Azarin, S. M., Sarkar, C. A. & Batchelor, E. Pulsed stimuli enable p53 phase resetting to synchronize single cells and modulate cell fate. *Mol. Syst. Biol.* **21**, 390–412 (2025).

67. Chickarmane, V., Ray, A., Sauro, H. M. & Nadim, A. A model for p53 dynamics triggered by DNA damage. *SIAM J. Appl. Dyn. Syst.* **6**, 61–78 (2007).

68. Novák, B. & Tyson, J. J. Design principles of biochemical oscillators. *Nat. Rev. Mol. Cell Biol.* **9**, 981–991 (2008).

69. Shimizu, T. S., Tu, Y. & Berg, H. C. A modular gradient-sensing network for chemotaxis in Escherichia coli revealed by responses to time-varying stimuli. *Mol. Syst. Biol.* **6**, 382 (2010).

70. Shoval, O. et al. Fold-change detection and scalar symmetry of sensory input fields. *Proc. Natl Acad. Sci. USA* **107**, 15995–16000 (2010).

71. Morton de Lachapelle, A. & Bergmann, S. Precision and scaling in morphogen gradient read-out. *Mol. Syst. Biol.* **6**, 351 (2010).

72. Goentoro, L. & Kirschner, M. W. Evidence that fold-change, and not absolute level, of β-catenin dictates Wnt signaling. *Mol. Cell* **36**, 872–881 (2009).


73. Liu, F., Morrison, A. H. & Gregor, T. Dynamic interpretation of maternal inputs by the Drosophila segmentation gene network. *Proc. Natl Acad. Sci. USA* **110**, 6724–6729 (2013).

74. Nikolić, M. et al. Scale invariance in early embryonic development. Preprint at arXiv:2312.17684 https://arxiv.org/abs/2312.17684 (2023).

75. Gregor, T., Bialek, W., de Ruyter van Steveninck, R. R., Tank, D. W. & Wieschaus, E. F. Diffusion and scaling during early embryonic pattern formation. *Proc. Natl Acad. Sci. USA* **102**, 18403–18407 (2005).

76. Cheung, D., Miles, C., Kreitman, M. & Ma, J. Adaptation of the length scale and amplitude of the Bicoid gradient profile to achieve robust patterning in abnormally large Drosophila melanogaster embryos. *Development* **141**, 124–135 (2014).

77. Stewart-Ornstein, J. & Lahav, G. p53 dynamics in response to DNA damage vary across cell lines and are shaped by efficiency of DNA repair and activity of the kinase ATM. *Sci. Signal.* **10**, eaah6671 (2017).

78. Gillies, T. E. et al. Oncogenic mutant RAS signaling activity is rescaled by the ERK/MAPK pathway. *Mol. Syst. Biol.* **16**, e9518 (2020).

79. Jimenez-Asins, A. et al. Time-series transcriptomics and proteomics reveal alternative modes to decode p53 oscillations. *Mol. Syst. Biol.* **18**, e10588 (2022).

80. Wang, Y. et al. Deep mutational scanning reveals EGFR mutations conferring resistance to the 4th-generation EGFR tyrosine kinase inhibitor BLU-945. *npj Precis. Oncol.* **9**, 294 (2025).

81. Sanchez-Vega, F. et al. Oncogenic signaling pathways in The Cancer Genome Atlas. *Cell* **173**, 321–337 (2018).

82. Finzel Pérez, A. Upstream control and downstream responses of p53 are involved in its tumor suppression functions upon genotoxic stress. PhD thesis, Freie Universität Berlin (2016).

83. Muzzey, D., Gómez-Uribe, C. A., Mettetal, J. T. & van Oudenaarden, A. A systems-level analysis of perfect adaptation in yeast osmoregulation. *Cell* **138**, 160–171 (2009).

84. Lemière, J. & Chang, F. Quantifying turgor pressure in budding and fission yeasts based upon osmotic properties. *Mol. Biol. Cell* **34**, ar133 (2023).

85. Silberberg, J. M., Kéer, S., Böhm, P. J. N., Jordan, K., Wiénberg, M., Grass, J. & Hänelt, I. KdpD is a tandem serine histidine kinase that controls K$^+$ pump KdpFABC transcriptionally and post-translationally. *Nat. Commun.* **15**, 3223 (2024).

---

**Acknowledgements.** This work re-analyses public single-cell and pharmacological data, and we thank the teams who produced, curated and shared these records. We thank E. Batchelor and colleagues, G. Lahav and colleagues, G. Mönke and colleagues, S. Tay and colleagues, and J. Albeck and colleagues for the p53 and NF-$\kappa$B live-cell recordings and reported statistics on which the pulse-encoder analysis rests; A. S. Hansen and C. Zechner for the Msn2 promoter-by-dose dataset (Zenodo 10.5281/zenodo.2755026); M. Son and colleagues for the NF-$\kappa$B spatial-gradient dataset (Zenodo 10.5281/zenodo.6858118); A. G. Wang, S. Tay and colleagues for the NF-$\kappa$B sequential-stimulation trajectories; S. Chavez-Abiega and colleagues for the GPCR-KTR dose series; A. Keshelava, V. L. Katanaev and colleagues for the GPCR-to-Ca²⁺ single-cell archive (27 experiments, per-cell peaks across seven doses with five repetitions), complete enough to support both the boundary-predictor adjudication P5 and the two-stage cascade analysis; J. P. Moore, K. Kamino, T. S. Shimizu, T. Emonet and colleagues for the chemotaxis FRET dose-response archive, released CC0 with foreground, background and mixed-ligand arms complete enough to support prospective blind tests (Dryad doi:10.5061/dryad.nvx0k6dzz); C. Klein Herenbrink and colleagues, L. M. Wingler and colleagues, and A. Gillis and colleagues for the D2R, AT1R and $\mu$OR primary pharmacology datasets, and E. L. Stahl and L. M. Bohn for their independent reanalysis; and the Biased Signaling Atlas team for maintaining an open, machine-readable pharmacology database whose archived provenance made a field-level identity check possible at all; and A. Morton de Lachapelle and S. Bergmann for the per-boundary precision and scaling tables of the Drosophila gap-gene read-out, released as complete supplementary datasets; and L. Goentoro and M. W. Kirschner for the Wnt fold-change study whose figure-level record supported quantitative digitisation. We also thank the maintainers of Dryad and Zenodo for the public infrastructure on which pre-registered reanalysis depends. We thank T. Gregor, W. Bialek and colleagues for archiving the per-embryo Bicoid and gap-gene profiles (Liu et al. 2013; Petkova et al. 2019) in the open scale-invariance repository that made adjudication P9 possible.

We also thank T. E. Gillies, M. Pargett, J. G. Albeck, F. McCormick and colleagues for the isogenic Ras/Raf single-cell EKAR3 archive with per-condition source data complete enough for information-theoretic reanalysis; A. Jimenez-Asins, G. Lahav and colleagues for releasing the matched-exposure p53 transcriptome and proteome repository; Y. Wang, J. T. Poirier and colleagues for depositing the EGFR saturation-mutagenesis codon counts (GEO GSE305057); and the cBioPortal and TCGA Pan-Cancer Atlas teams for maintaining open, harmonised tumour-level molecular data.

**Competing interests.** The author declares no competing interests. This preprint reports an independent research project of the author; it received no dedicated funding and is not associated with any product or service of the affiliation.

---

## Figure captions

**Fig. 1 | One discipline, four systems, one ledger.** a, The five-step analysis discipline and the degeneracy-group concept (four parameter sets, one curve; machine-precision invariance). b, The synthesis matrix: rows = four systems; columns = absolute-scale identifiable / internally consistent / information-bearing versus scale-free identifiable / consistent / information-bearing, with the archived verdict for each cell. c, The blind-adjudication ledger (P2–P7, fourteen clauses): hits, intermediates and falsifications, each with dataset and frozen decision line.

**Fig. 2 | The static sensor.** a, PdPC steady-state curves across $\kappa$. b, Degeneracy-group invariance (four sets, one curve, max $|\Delta u| = 0$). c, Fisher eigenspectrum, shape versus scale tiers ($\chi \approx 30$ versus ${\sim}10^{300}$). d, Hill-slope closed form $n_{\mathrm{H}} = 1 + 1/(2\kappa)$ versus numerics; the one-sided-bound regime. e, Mutual information $I(\xi; u)$ versus $\kappa$ with interior optima; the zero-order limit as decision element.

**Fig. 3 | The dynamic encoder.** a, NF bifurcation diagram with $\sqrt{}$-onset (continuation-grade; $A^2$-linearity $R^2 = 0.9926$). b, Subcritical hysteresis loop, the third signature (0.275 versus 0.25). c, Discrimination diagram with model-side and experimental-side placements (p53-DSB $\alpha = -0.001 \pm 0.113$; p53-UV $\alpha = +0.981 \pm 0.114$; 6.1 combined $\sigma$; ERK, NF-κB, Msn2-type placements). d, Channel census on the excitable representative under the stated calibration: 94.2% timing versus 5.5% amplitude (mechanism-level recheck preserves the ordering; SI S13.3); amplitude two-layer death ($\mathrm{CV}(A)$ $0.002 \to 0.639$; $\mathrm{CV}(T)$ pinned $0.016$). e, The digital limit: chemical saturation ($\kappa \to 0$) and circuit topology ($\alpha \to 0$) as two roads to the same likelihood cancellation; counting law $N \approx (\tau_r/T) \cdot \ln(D_0/D_c)$ with $\tau_r \approx 7.7$ h on the digitised scale (absolute dose scale convention-dependent; shape corroborated at mechanism level, SI S13.3).

**Fig. 4 | Two roads to one digital limit.** Schematic of the likelihood-cancellation argument. Left: the static covalent-modification sensor reaches calibration-freeness by chemical saturation ($\kappa \to 0$); the dynamic pulse encoder reaches it by circuit topology ($\alpha \to 0$). Both routes cancel the analogue scale factor from the likelihood: the Fisher matrix collapses onto a single partition-counting direction, residual noise becomes combinatorial, and only dimensionless ratios survive. The analogue regime between the two limits is information-rich but calibration-hungry.

**Fig. 5 | Blind adjudications bound the claim.** a, P2 (Msn2): stratification penalty of the molecule-number channel (21/21 units centred on zero) versus event-time insensitivity ($|\Delta| \leq 0.01$). b, P3 (NF-κB gradient): spatial ordinality (AUC 0.907/0.819), dose saturation failure, duration-axis falsification. c, P4 (ERK–KTR): event-time decoding AUCs below line across two ligands; stratification-invariance falsification. d, P5 (GPCR→Ca²⁺): dead-zone violation and weak per-cell decoding against monotonic population-level counting. e, P6/P7 fold-change attribution: power-limited intermediates on NF-κB sequential stimulation ($\Delta R^2 = -0.0014$) and the chemotaxis background ladder ($\Delta R^2 = +0.118$ / $-0.091$ across the erratum pair). Frozen decision lines marked on every panel.

**Fig. 6 | The receptor-layer rupture and the prospective test.** a, Pareto front of joint fits (no point in the acceptance region; $N \approx 2$–$4$ versus $N \approx 12$–$39$). b, The exact Weber line $K_{1/2} = 1.17 \cdot (1.95 + B)$ µM across seven backgrounds (max deviation 0.10 dex). c, Eight repair paths, all excluded (summary grid). d, Same-source verdict: the front is unmoved when both statistics come from strictly the same cells. e, Prospective mixed-ligand arm: blind predictions on the scale-free midpoint coordinate hit (+0.014/+0.038/+0.031 dex within the frozen $\pm 0.15$ dex band); blind predictions on the amplitude coordinate miss (0.105/0.126 against the frozen 0.08 band).

**Fig. 7 | The pharmacological meta-analysis.** a, D2R fracture map (bifeprunox $\Delta = -3.5$ dex; load-bearing coincidence). b, AT1R cross-arm fracture ($\Delta = -2.05$ dex, 5.0$\sigma$; binding affinity between the arms). c, μOR GIRK-column systematic fracture (two identities locate the same column) and its repair by the coordinate-aware minimal extension ($\Delta\log\mathrm{ML} = +128.2$; residual $0.811 \to 0.107$ dex). d, Atlas-wide identity violations (21.5% > 0.3 dex; 4.2% > 1 dex), the constraint-projection repair map, and the silently dropped column.

**Table 6 | The currency per system (the free-first pattern).**

| System | Carrier statistic (survival set) | Class | Paid reference | Evidence |
|---|---|---|---|---|
| p53 (DSB response) | dimensionless timing $t_1/\mathrm{IPI} \approx 5.7$ reproduced (literature 5.8); pulse count $N$ | I | none detected; noise scale $\sigma$ is a Fisher null direction | S3, S4, S14 |
| NF-κB | fold change; pulse count increases with dose | I | none | S6.1, S7.5 |
| ERK | digitised event presence and duration | I | none | S6.3, S7.3 |
| Msn2 | event frequency | I | none | S6.5, S7.1 |
| Bacterial chemotaxis | Weber fold change $\Delta L/L$ | I (readout) | methylation-maintained effective $K_d$ (class III) | S8, S7.6 |
| GPCR pharmacology (field level) | none portable across laboratories | negative result | parameters drift unchecked | S9 |
| PdPC static sensor | ratio and threshold $g^* = 0.443$ | I | none | S2 |
| Bicoid to gap genes (Drosophila embryo) | relative position $x/L$ read off Bcd threshold crossings | I (readout) | scaling paid by nuclear-density counting ($\propto N/L^2$, a geometric counter); anterior hyper-scaling ($S \approx 1.7$) and posterior hypo-scaling ($S \approx 0.89$) locate the terminal anchors | S1.6b, S15 |
| Wnt/β-catenin | fold change over background (FCD); plateau 0–50 mM LiCl within ±5%, breakdown at 50–58 mM | I | synthesis–degradation ratio of the destruction complex (no maintained reference) | S16 |
| Yeast osmotic homeostasis (positive control) | absolute volume and turgor | III | internal glycerol set point held by the Hog1–glycerol–Fps1 integrator loop; all three markers fire; Hog1 inhibition abolishes recovery (0.644 to 0.070) | S1.6c |

Class I: calibration-free statistics (ratios, counts, event times); class III: actively maintained paid references. The right-hand pattern, free carriers everywhere and payment concentrated at functionally necessary anchors, is the cost-effectiveness corollary of SI S1.6b.

**Table 7 | Reverse predictions from the machine line.** Three verdict-grade laws from the companion artificial-perception programme, each exporting a falsifiable cellular prediction, frozen in advance of testing.

| Law | Cellular prediction | Falsifier | Status |
|---|---|---|---|
| Cheapest sufficient invariance | a pathway's coding regime exactly matches the pricing of its perturbation environment | systematic over-insurance | frozen; the P2 adjudication ruled out the size-gain operationalisation and the clause has contracted accordingly |
| Estimator ordering | sensory adaptors pay content-noise prices while homeostatic setpoints pay reference-drift prices | the two sensitivities do not decouple | frozen, awaiting test |
| Two-factor code switching | coding regime is insensitive to noise severity but switches under noise-type swap at matched severity | no switch under a direct synthetic-biology assay | frozen; cheapest of the three to test |

**Table 8 | Boundaries and limitations.** Seven considerations bound the present claims.

| # | Boundary | What it bounds |
|---|---|---|
| 1 | Analogue information is not excluded in principle; where cells report graded quantities, the payment theorem predicts the calibration cost is paid explicitly, and morphogen biology reads as a catalogue of such payments | the scope of the scale-free reading |
| 2 | Information statements concern steady-state, single-readout channels | the digital-limit principle bounds what is robustly transferable per unit of calibration effort, not what is transferable in principle |
| 3 | Lemma 2 assumes time-scale separation with fast-subsystem folds; Lemma 1's Hopf-onset argument is genericity-based | the reach of the topology-classification lemmas |
| 4 | The chemotaxis rupture is a statement about a model class, not the bacterium | structural versus statistical repair is a registered open problem with a prospective arm awaiting the necessary dataset |
| 5 | The GPCR meta-analysis quantifies inconsistency, not misconduct | its deliverable is an inconsistency map plus a reconciliation checklist that any extracted bias factor should pass before use |
| 6 | The p53 molecular-layer decomposition carries three registered boundaries (SI S14.5): the eigenvalue period overshoots simulation by 8.5%, the simulated width sits below the archived band, and the Wip1 width response conflicts between two literature sources | the molecular-layer closure; the first-pulse and interpulse closures remain anchored on the same measured bands |
| 7 | Gain-distribution confinement does not generalise beyond the M3R-to-Ca cascade (pre-registered ratios 0.59–2.17 against the frozen line of 3 on ERK and NF-$\kappa$B; SI S17c) | amplitude in those systems dies by within-cell pulse-to-pulse stochasticity; the carrier-level conclusion (frequency and count codes) is unaffected |


**Table 9 | Spin-off discoveries (five sharpest; full ledger in SI S18).**

| Spin-off | Statement | Evidence |
|---|---|---|
| Repair time from the counting-law slope | The slope of the p53 pulse-counting law inverts the DSB repair timescale ($\tau_r \approx 7.7$ h), turning a fit parameter into a measurand; no $\gamma$H2AX assay needed | S14 (Lemma 3 by-product) |
| A computable topology diagnostic | The dispersion ratio $t_1/\mathrm{IPI}$, literature value 5.8 (reproduced at 5.7), is impossible for a bare oscillator (bound 0.03–0.12) and certifies an independent stochastic sensing stage upstream | S4, S14 |
| Rupture is confined to the gain distribution | In the M3R-to-Ca²⁺ cascade, across-experiment $\mathrm{EC}_{50}$ scatter (32-fold) decomposes into across-cell gain dispersion with 3.5% within-cell repeatability; the cell is precise, the population is not | S17 |
| The paid reference sits at the absorbing relay | In chemotaxis, the single paid reference (methylation-maintained affinity) occupies exactly the middle stage that resets the upstream parameter rupture (baseline restoration to 2.2%, history-independent) | S17b, S8 |
| A field database fails identity checks | 21.5% of Biased Signaling Atlas records violate the operational model's own cell-level identity by more than $0.3$ dex, and its inclusion choices silently dropped the fractured $\mu$OR column | S9 |


*Extended Data figures Extended Data Figs. 1–16 (per-system analysis details; asymmetric-Michaelis and positive-feedback-switch panels; dispersion signatures and noise-colour plane; counting-channel capacity; digital-limit theorem diagrams; NF-κB dual-regime analysis and the two rejected patches; Morris–Lecar duality instance table; full blind-adjudication unit tables; chemotaxis repair-path details and per-background Weber residuals; GPCR per-study fracture maps and tracer-robustness; analysis artefact chain and double-recorded failure ledger; p53 molecular-layer decomposition; closed-form cross-validation; cancer-chain verdict ledger).*
