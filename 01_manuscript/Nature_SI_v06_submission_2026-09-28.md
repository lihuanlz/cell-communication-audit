# Supplementary Information

**Scale-free statistics are the currency of cellular communication**

Huan Li

Independent Researcher, Shanghai 201899, P. R. China

Corresponding author: HL@liangtabio.com (H.L)



Companion to the main text. Every number in this document is transcribed from the archived verdict cards and script outputs listed in S0; nothing is quoted from memory.

---

## S0. Archive map and reproducibility contract

### S0.1 Deposit map

The deposit snapshot (github.com/lihuanlz/cell-communication-audit; archived at Zenodo, doi:10.5281/zenodo.23008824) is organised by function:

- `01_manuscript/`: manuscript and SI, markdown sources of record and docx renders.
- `02_figures/`: every main and Extended Data figure (PNG and SVG).
- `03_data_csv/`: every numeric table behind the figures and verdicts (CSV/JSON), including `p53_molecular_layer/` (S14) and `tumour_chain_P11_P16/` (S19).
- `04_code/`: every audit script cited in this paper, flat with zero-padded filenames (`code04_...` to `code93_...`), plus `CODE_MAP.csv` (Table S0, S0.4) and the raw chemotaxis FRET `.mat` files; subdirectories hold the Extended Data figure scripts (`ed_scripts/`), the 2026-09-27 pre-registered batch (`preregistered_audits_2026_09_27/`, S15-S17), the tumour-chain scripts (`tumour_chain_P11_P16/`, S19) and frozen assets (`assets/`).
- `06_supplement_interactive/app/`: the interactive Atlas inconsistency map.
- `07_audit_cards_tumour_chain/`: frozen verdict cards P11-P16.
- `08_audit_trail/`: the audit-trail documents cited in this SI (verdict cards, adjudication reports, reconciliation notes, pre-registration texts), stored under their historical workspace-alias paths so that every path quoted in this SI resolves verbatim inside the deposit.
- Root: `README.md`, `05_validation_records_index.md`, and `_archive_alias_map_2026-09-26.csv`, which records the historical workspace folder names used during the project; the deposit itself uses the layout above.

### S0.2 Pre-registration and verdict discipline

Every empirical claim in the adjudication programme carries one of two archival statuses. *Sealed* means the decision table, thresholds, and scoring script were frozen before the corresponding data values were read. *Registered* means the analysis plan was written before execution but the data had been seen during pipeline development. Frozen verdicts never move: when an error in an executed clause was discovered after freezing, the correction was appended as an erratum that re-adjudicates under the corrected rule, with the original verdict preserved alongside (errata 44b for P2, 47b for P5, 49b for P7). Successes and failures are double-recorded in the same ledger; the full failure inventory is given in S10.2.

### S0.3 Script-number reconciliation (closed 2026-09-23)

The absorbed NC lineage (SI v1.5) promised reproducibility via "scripts 4, 6, 9, 10, 11, 12, 13, 15, 16, 17". Scripts 16 and 17 have no same-named physical files under `08_audit_trail/01_cell_line_1/`: the NC numbering was a logical numbering of that manuscript lineage, and the physical artefacts live under different numbers in `08_audit_trail/03_cell_line_3/`. The promise was not empty, but the numbering requires reconciliation. The reconciliation is frozen in the table below; the deposit-side mapping of record is `04_code/CODE_MAP.csv` (S0.4). Its content:

| NC logical number | Content | Physical file |
|---|---|---|
| code 4 | PdPC identifiability audit (Theorems 1–3 numerical cores) | `04_code/code04_pdpc_identifiability_audit.py` |
| code 6 | NF-κB negative-feedback oscillator audit (soft mode) | `04_code/code06_nfkB_negative_feedback_oscillator_audit.py` |
| code 9 | NF-κB spiky-mode dose audit (with erratum log) | `04_code/code09_nfkB_spike_mode_audit.py` |
| code 10 | asymmetric-K_m degeneracy group analysis | `04_code/code10_asymmetric_Km_degeneracy_group.py` |
| code 11 | discrimination diagram v0.9 (experimental placements) | `04_code/code11_discrimination_map_v09.py` |
| code 12 | NF-κB class-transfer parameter unreachability test (rejected patch 1) | `04_code/code12_nfkB_GLM_prediction_test_switchlike_IkB.py` |
| code 13 | p53 Fisher stratification (audit step-3 closure) | `04_code/code13_p53_fisher_stratification_S7.py` |
| code 15 | NF-κB saturated-transient hypothesis test (rejected patch 2) | `04_code/code15_nfkB_saturation_transient_test.py` |
| code 16 | p53 channel-level identifiability census (94.2%/5.5%, λ₇/λ₁ bound, 96.8%) | `04_code/code51r3_p53_fisher_audit.py` (r1/r2 retained alongside) |
| code 17 | Track C2 source-digitisation bootstrap re-check | `04_code/code52_TrackC2_discrimination_diagram.py` |

All scripts above, the adjudication scripts 44/44b/45/46/47/47b, the p53 census family (51r1-r4, 52), the chemotaxis arm (51chemo, 53-62, 70), the GPCR codes (63-73, 77, 80) and codes 82-93 (S13.3/S14) live flat in `04_code/` of the deposit. One naming collision is registered: `code52_protocol_artifact_test_raw_FRET_coordinate_analysis.py` and `code52_TrackC2_discrimination_diagram.py` are two different scripts from two lines sharing the number 52; the deposit uses full filenames throughout, and `CODE_MAP.csv` disambiguates every row.

### S0.4 Table S0

Table S0 is `04_code/CODE_MAP.csv` in the deposit: one row per archived script or figure generator, giving the package filename, the original workspace path, the NC logical number where the absorbed lineage had one, the SI section where its result is reported, and a one-line description.

---

---

## S1. The mathematical framework

This section states the framework behind Results section 1, condensed from the archived framework document (v4.2, 2026-08-15). Status at the outset: the information-loss identity and mutual-information monotonicity (S1.7) are theorems given the model of S1.1; the survival set (S1.3) and the three-leg criterion (S1.5) are operational definitions; no tightness, reachability, or capacity claim is made.

### S1.1 Layer model and the degeneracy group

The input is a stochastic process $S(t)$ from an input ensemble $\mathcal{S}$ with distribution $p_S$. An input statistic $X = f(S)$ (a difference ratio, sign, count, event time, amplitude) is defined independently of layer parameters. Layer $i$ is a stochastic map: given its upstream input and parameters $\theta_i \in \Theta_i$, its output symbol is distributed as $Y_i \sim P(Y_i \mid \text{upstream}, \theta_i)$, with biochemical noise absorbed into $P$. A cascade is a Markov chain $L_1 \to \cdots \to L_n$, $Y_0 = S$, in which $Y_i$ depends only on $Y_{i-1}$ and $\theta_i$. The degeneracy group of layer $i$ is the group of parameter transformations leaving the output distribution invariant:

$$G_i = \left\{ g : \Theta_i \to \Theta_i \;\middle|\; P(Y_i \mid Y_{i-1}, \theta) = P(Y_i \mid Y_{i-1}, g\cdot\theta),\ \forall\, Y_{i-1}, Y_i \right\}. \tag{S1.1}$$

The identifiable content of the layer is exactly the set of orbit invariants of $G_i$: two parameter values on one orbit cannot be distinguished from the output at any measurement precision. For a fixed statistic $X$, the conditional distribution is $P(Y_i \mid X, \theta_i) = \int P(Y_i \mid S, \theta_i)\, dP(S \mid X)$, with $dP(S \mid X)$ induced by the input ensemble.

### S1.2 The layer-boundary rule

Layers are set by information structure, not biochemical topology. The rule (restored in v4.2 from v3.3 §2.2, and required for pathways with feedback): **layer boundaries are degeneracy equivalence classes**. Parameters unidentifiable from the same output belong to one layer; a feedback parameter is assigned to the layer from whose output it is unidentifiable. Cross-layer feedback (adaptation canonically) is thereby absorbed into the corresponding layer, and the Markov property of S1.1 holds under this assignment. Without this rule the framework cannot be applied to real pathways containing feedback.

### S1.3 The survival set

Let $a_i$ be the **transmitted readout functional** of layer $i$, the alphabet map the downstream layer consumes (leg 2, S1.5); let $\Theta_i^{\mathrm{work}}$ be the physiological working range and $\varepsilon$ calibrated to the system noise level. The survival set is

$$T_i = \left\{ X \;:\; \sup_{\theta \in \Theta_i^{\mathrm{work}},\, g \in G_i}\ \mathrm{TV}\Big( P\big(a_i(Y_i) \mid X, \theta\big),\ P\big(a_i(Y_i) \mid X, g\cdot\theta\big) \Big) < \varepsilon \right\}. \tag{S1.2}$$

Two points carry the weight of a recorded correction (Bug 7, S1.8). First, the total-variation condition acts on the distribution of $a_i(Y_i)$, not on the raw output distribution: acting on $P(Y_i \mid \cdot)$, any unidentifiable gain (the generic case) makes every output distribution sensitive to $g$, the survival set is identically empty, and the definition contradicts the interface algebra of S1.4, under which difference ratios and signs survive unconditionally. Second, the gain-degeneracy note: $G_i$ in (S1.1) is defined on conditional distributions given the input, and strictly a pure gain $k$ changes the conditional distribution, so it is not in $G_i$. The true source of gain degeneracy is **joint confounding**: with both channel gain $k$ and input ensemble $p_S$ unknown, $(k, p_X)$ and $(ck, p_{X/c})$ induce identical output distributions, and the degeneracy lives in the "input ensemble unknown" direction. The alphabet-level TV condition of (S1.2) is the operational projection of this joint degeneracy; the two-level division of labour ($G_i$ for parameter degeneracy under known input, (S1.2) for the observable consequences of joint degeneracy) is a stated convention.

### S1.4 Interface algebra

Which statistics survive depends on the function class of the interface $\varphi_i$. The generic biological interface is affine, $y = kx + b$, with $b$ a basal term (spontaneous activation, basal expression, leak) and the gain $k$ generically unidentifiable. Fold-change does not survive: in $(kx_1 + b)/(kx_0 + b)$ neither $k$ nor $b$ cancels, and fold survives only at $b = 0$. Difference ratios survive: $\Delta y = k\,\Delta x$, so $(y_3 - y_2)/(y_2 - y_1) = (x_3 - x_2)/(x_2 - x_1)$ with $k, b$ both cancelling; the sign survives for $k > 0$; counts and event times survive because threshold-crossing order and timing are affine-invariant. The unconditional survivors are therefore difference ratios, signs, counts, and event times. Amplitude survives if and only if the gain is pinned from outside; an increment survives where difference ratios survive and the gain is pinned to the precision the downstream consumer requires, $|\Delta k/k| \le \delta$. Since basal terms are nearly ubiquitous, the generic survival set is $\{\text{sign, increment, difference ratio, count, event time}\}$; the fold-change literature (the Goentoro/Shoval line^16,47^) holds where $b \approx 0$ or the interface is actively normalised, a special rather than generic case.

### S1.5 The three-leg criterion

$X$ crosses interface $\varphi_i$ if and only if three conditions hold simultaneously. Leg 1, **algebraic survival**: $X \in T_i$ in the sense of (S1.2); failure mode, scale leak. Leg 2, **symbolic reachability**: $I(X; Y_i) > 0$; failure mode, the statistic is invariant but uninformative. Leg 3, **temporal bandwidth**: the statistic's timescale $\tau_X$ lies inside the interface bandwidth $\mathrm{BW}(\varphi_i)$; failure mode, erasure by filtering. The leg-2 alphabet table maps interface type to transmitted alphabet: linear/affine interfaces transmit continuous concentration; a switch transmits $\{\text{above}, \text{below}\}$ (a quantiser); a pulse generator transmits $\{\text{pulse}\} \times \{\text{timing}\}$. The leg-3 bandwidth is defined in the time domain: modulate the input at timescale $\tau$ (sinusoidal or square wave) and take the set of $\tau$ whose output response stands above the noise floor by a factor of at least 1; this is directly measurable in perturbation–response experiments, and the frequency-domain equivalent holds only where the interface is linearisable. **Adaptor algebra identity:** adaptation is the device that actively subtracts the baseline and pulls the affine interface back onto the increment channel; an adaptor outputs increments (which still require pinned gain), not difference ratios directly. This is the algebraic identity of adaptation and explains its evolutionary prevalence: what an adaptor purchases is the survival of the difference family.

### S1.6 Three reference classes and the prospective predictions

The naive dichotomy, free scale-free statistics versus paid absolute quantities, has a hole: molecular constants ($K_d$, $K_m$) are absolute in units yet free of cellular maintenance, so a signature "wherever absolute scale is referenced there must be active maintenance" would flag every $K_m$-driven reaction as a false positive. The classification axis is therefore not the object type but **who sets and maintains the reference, on what timescale**. Class I, calibration-free statistics (difference ratios, signs, counts, event times): no reference needed. Class II, structural constants ($K_d$, $K_m$): set by natural selection, maintained by physics (molecular structure); on physiological timescales the cell neither pays nor can adjust them; they are dimensional converters mapping an external absolute concentration $[L]$ to the internal dimensionless occupancy $[L]/(K_d + [L])$, and trigger no maintenance signature. Class III, paid references (set points, internal standards, and physiologically regulated effective constants such as the methylation-adjusted effective $K_d$ of the chemotaxis receptor or the FliM-remodelled effective CheY-P threshold of the motor): set by the cell, maintained by active feedback, drifting if unmaintained, paid for continuously.

A site is class III when at least two of three operational markers hold, with reference-value specificity required: (1) negative feedback acts directly on the reference set point itself; (2) direct perturbation of the reference value is followed by recovery to the same set point (recovery of a downstream response does not count); (3) maintaining the reference consumes extra energy, and inhibiting the energy supply makes the reference drift. The **payment position** is \{points where function references absolute scale\} $\cap$ class III, determined by functional demand, not by where information is lost: a terminal decision referencing only relative quantities pays nothing; one referencing absolute position must establish and maintain a local reference there. The class-III signature is a forward search instruction, not an executable reverse falsification: an exhaustive search finding no maintenance triggers re-examination rather than a verdict against the framework, because "not found" and "absent" cannot be separated in practice; class-II references lie outside its jurisdiction.

Three prospective falsifiable predictions follow. For steady state $S_0$ and perturbation $\Delta S$, define absolute-quantity sensitivity $dY/dS_0$ and relative-quantity sensitivity $dY/d(S/S_0)$; kill thresholds are fixed per system before data are opened (form fixed, no numerical value pre-committed). **P-1 (deletion):** deleting a class-III maintenance mechanism predicts $dY/dS_0 \to 0$ with $dY/d(S/S_0)$ preserved; retained absolute sensitivity with no substitute maintenance falsifies the framework. **P-2 (synthetic circuit):** a circuit with no active internal standard (no class III) but passive structural constants (class II) predicts $dY/dS_0 \approx 0$ with nonzero relative sensitivity; stable absolute-scale transmission falsifies the framework. **P-3 (search target):** generate the list of positions that should be class III but have no reported maintenance (for example, cells reading absolute position in a morphogen field); a hit confirms, an exhaustively documented absence triggers re-examination.

### S1.6b Corollary: the free-first (cost-effectiveness) principle

**Statement.** Wherever a class-I statistic suffices for the function, it is the one used; class-III payment occurs exactly at the positions where function references absolute scale, and nowhere else. Formally this is the contrapositive sharpening of the payment-position rule above: the payment set is not a subset of the demand set but equal to it. **Falsifiable form:** a single position whose function references only relative quantities yet which hosts an active, energy-consuming reference-maintenance mechanism falsifies the corollary; this direction is complementary to P-3, which searches for unpaid positions that should pay, while the corollary excludes paid positions that need not pay.

**Worked example: bacterial chemotaxis, both ends measured.** The free end is the readout: single-cell FRET measurements show the pathway responds to the fold change of attractant, with output dynamics invariant to multiplying the absolute ligand level by a constant over about two orders of magnitude^62^, implementing the fold-change detection predicted from the circuit's nonlinear integral feedback^69,70^. The paid end is exactly one reference: the receptor's effective affinity is held by the CheR/CheB methylation feedback, which satisfies all three class-III operational markers (the feedback acts on the reference set point itself; perturbation of methylation level is followed by recovery to the same set point; the adaptation consumes cellular resources and the reference drifts when adaptation is disabled). One system therefore exhibits the full corollary in a single view: everything that can be free is free (ratio readout), and payment is concentrated at exactly one functionally necessary anchor.

**Second example: p53 dose encoding by count, not amplitude.** For double-strand-break induction, pulse amplitude and duration are fixed independently of dose while pulse number grows with dose^9^, and pulse dynamics, not amplitude level, controls downstream fate^11^. Our channel census independently finds 94.2% of Fisher information in the timing channel against 5.5% in amplitude, with the noise scale $\sigma$ an exact null direction (S4): the dimensionless carriers are used, the calibratable-amplitude channel is functionally silent, and the cell does not even pay to stabilise its own noise level. A twelve-cell-line live-cell survey^77^ (3825 single-cell trajectories, five irradiation doses) corroborates the division of labour from the population side: wherever p53 oscillates, the period is pinned near 5 h in every line, while dose and cell-line identity modulate only envelope properties (pulse width, amplitude, sustainability), and pharmacological shift of DNA-repair efficiency or ATM activity converts dynamics between pulsatile and sustained without moving the period. What varies across contexts lives in the amplitude channel; the time anchor does not move.

**Third example: PdPC does not overbuild.** The information-optimal operating point of the static sensor is interior, not at the zero-order ultrasensitive limit (S2.6): the system pays for exactly as much cooperativity as the information return justifies, and no more.

**Fourth example: morphogen scaling is paid for by a geometric counter.** Gap-gene boundaries in the Drosophila blastoderm are read at relative positions $x/L$, and our independent recomputation from the per-boundary scaling coefficients of Ref.^71^ (S15) shows that the scaling itself reduces to nuclear-density counting, a class-I geometric object: because Bcd degradation tracks nuclear density ($\propto N/L^2$ as nuclei divide), embryo length normalisation comes free with development and no dedicated length-measuring reference is maintained. The two registered deviations, anterior hyper-scaling ($S \approx 1.7$ at 1$\times$ bcd) and mild posterior hypo-scaling ($S \approx 0.89$), locate where terminal anchoring adds what the free mechanism cannot supply, consistent with the corollary's demand that payment appear exactly at functionally necessary anchors.

### S1.7 The information-loss identity and the weakest-layer theorem

The end-to-end survival set is defined through the composite map, guaranteeing split-independence: $S_{\mathrm{cascade}} = \{ X : X \text{ passes the three legs for } f = \varphi_n \circ \cdots \circ \varphi_1 \}$.

**Theorem (information-loss identity).** By the cascade Markov property $S \to Y_1 \to \cdots \to Y_n$, define the loss at layer $i$,

$$L_i = I(S; Y_{i-1}) - I(S; Y_i), \quad i = 1, \dots, n, \quad Y_0 = S. \tag{S1.3}$$

Then

$$H(S) = I(S; Y_n) + \sum_{i=1}^{n} L_i, \tag{S1.4}$$

and mutual information is monotone, $I(S; Y_1) \ge \cdots \ge I(S; Y_n) \ge 0$. The **weakest layer** $i^{*} = \arg\max_i L_i$ causes the largest loss of input information. For any $X = f(S)$, $I(X; Y_n) \le I(S; Y_n)$, and the data-processing inequality gives the per-statistic chain $I(X; Y_1) \ge \cdots \ge I(X; Y_n)$: end-to-end transferable content is bounded by the weakest-layer loss, which the identity locates without any capacity concept.

**Intersection bound (loose, labelled as such):**

$$S_{\mathrm{cascade}} \subseteq T_1 \cap \varphi_1^{-1}(T_2) \cap \cdots \cap T_n. \tag{S1.5}$$

Passing layer by layer is necessary but not sufficient end to end; the gap comes from loss of leg-2 sensitivity under composition (a downstream layer discards a component of the upstream symbol) and from alphabet re-quantisation. No tight bound is provided; the intersection's role is to locate the cause of death.

### S1.8 Development lineage

The lineage is part of the audit trail (v3.3 and v4.2 appendices). Versions v2.1/v3.0/v3.1 repaired six internal bugs; v3.2/v3.3 answered two external-review rounds of sixteen items; v4.1 was a condensed rewrite upgrading the information-loss identity and removing the capacity concept; v4.2 repaired two further items. The first is Bug 7, self-introduced in v4.1: leg 1's condition acted on the raw output distribution, under which the survival set is generically empty and contradicts the interface algebra; the fix moved the TV condition to the readout-functional distribution (S1.3). The second is the restoration of the layer-boundary rule (S1.2), which the v4.1 condensation had dropped. The full genealogy of internal bugs 1–6 and the sixteen external items is archived in the v3.3 appendix.

### S1.9 Framework self-check simulations

Before any system audit, the three legs and the identity were checked numerically on minimal stochastic models, failures double-recorded.

**p53 cascade (code 50, two-arm v0.2).** Three layers (damage sensing; p53–Mdm2 pulse train; promoter readout), affine observation $y = k\cdot\mathrm{p53} + b$, per-cell unidentifiable $k \sim \mathrm{LogNormal}(0, \sigma_k)$, $b = 0.3$; arm A digital (dose in pulse count, amplitude $1.0 \pm 0.15$), arm B analogue (amplitude growing linearly with dose), so that arm-B amplitude carries real information and its death strongly tests leg 1. Sweeping $\sigma_k$ from 0 to 0.6 (4000 cells per point): arm-B amplitude information collapses 1.437 $\to$ 0.158 (scale leak confirmed); fold-change decays in both arms (arm B 0.376 $\to$ 0.177), refuting the auditor's expectation that ratios cancel the gain, because with $b \neq 0$ the gain does not cancel, agreeing with the algebra row "fold survives only at $b = 0$" (double-recorded, neither side's numbers altered); counts stay flat (0.37–0.51), immune to gain degeneracy. Leg 3, sweeping a downstream integrator: transmitted count information reads 0.79, 1.11, 0.34, 0.14 at $\tau = 0.5, 1, 3, 8$ h, collapsing as $\tau$ approaches the 5.5 h pulse spacing; bandwidth mismatch kills a surviving statistic. The identity closes bit-exactly on arm A: $H(D) = 2.322$, $I(D; N) = 1.729$, $I(D; Y_{\mathrm{count}}) = 0.435$, $L_1 = 0.593$, $L_2 = 1.294$, and $0.435 + 0.593 + 1.294 = 2.322$; the weakest layer is the readout layer (information dies at event detection, not encoding). Limitations are recorded: toy model; mutual information from binned interpolation (absolute values biased, trends reliable); agreement with data is a consistency check, not a test; the tests are P-1/P-2/P-3.

**Chemotaxis receptor layer (v0.3, codes 51–53 with later rounds).** Checked against the P7 erratum unit table (30 background-by-foreground units) as calibration reference. Three findings matter for the framework. First, the naive MWC occupancy-increment kernel is falsified: with $\theta(L) = L/(K_d + L)$ the predicted response is multiplicative, $R^2(\log F) = 0.000$ against the empirical 0.652, predicting 0.003 against an empirical 0.424 at $B = 100$, $F = 40$; occupancy increment is not the receptor-layer response kernel. Second, raw amplitude leaks under gain degeneracy ($R^2(\log F)$ falls 0.687 $\to$ 0.192 as $\sigma_k$ goes 0 $\to$ 0.6), while per-cell normalisation survives by trivial algebra alone: the external audit's point that analyst-side normalisation is the analyst cancelling $k$ on the cell's behalf was accepted and recorded as Bug 8, downgrading that arm from mechanistic evidence to an algebraic control. Third, the mechanistic core (MWC activity states with Barkai–Leibler exact adaptation) is natively fold-type: the cell-autonomous readout (increment from the cell's own adapted baseline, no analyst normalisation) survives gain degeneracy and adaptation noise with only mild degradation (worst-case $R^2(\log r)$ 0.576 $\to$ 0.356), the first mechanistic demonstration of class-III paid pinning, but what feedback preserves is the layer's native (fold) coordinate, not the empirical increment structure. The empirical kernel instead has the form $R = A \cdot \ln(1 + F/K_c)/(1 + (B/B_s)^p)$, $K_c \approx 0.46$–$0.57$ µM, $B_s \approx 63$–$70$ µM, $p \approx 0.8$–$0.9$ ($R^2 = 0.967$; the earlier bootstrapped phenomenological fit gave $R^2 = 0.942$), an absolute-increment form anchored at a structural constant that the standard free-energy-difference family does not generate. The recorded framework-level correction: "increment survival = difference-ratio survival + pinned $k$" was rewritten as "feedback pins baseline and gain, preserving the layer's native statistic", and whichever layer performs the fold-to-increment coordinate transformation carries an additional paid reference structure. All of these are self-consistency checks sharing one calibration table; the independent test is the registered P8 arm. The later rounds of this file (the amplitude-versus-$K_{1/2}$ joint inconsistency and its eight excluded repair paths) belong to SI S8, not to this self-check.

---

---

## S2. PdPC static-sensor analysis

The phosphorylation–dephosphorylation cycle (PdPC; Goldbeter–Koshland cycle) converts an input activity ratio into a steady-state modified fraction. With substrate total $S_T$, kinase $E_1$ (maximum rate $V_1$), phosphatase $E_2$ ($V_2$), and modified fraction $u$, the symmetric-Michaelis steady state obeys

$$\xi = \frac{u}{1-u}\cdot\frac{\kappa + 1 - u}{\kappa + u}, \qquad \xi \equiv \frac{V_1}{V_2},\quad \kappa \equiv \frac{K_m}{S_T}. \tag{S2.1}$$

The $\kappa \to 0$ limit is the zero-order switch of Goldbeter and Koshland; $\kappa \to \infty$ degenerates to a graded Michaelian response $\xi \approx u/(1-u)$.

### S2.1 Theorem 1 (degeneracy group): six physical parameters, two identifiable combinations

**Theorem 1.** The steady-state dose–response curve of Eq. (S2.1) depends on the six physical parameters $(k_1, E_{1T}, k_2, E_{2T}, K_m, S_T)$ only through the two dimensionless combinations

$$\xi = \frac{k_1 E_{1T}}{k_2 E_{2T}}\cdot x \;\;(\text{titration variable } x), \qquad \kappa = \frac{K_m}{S_T}.$$

Equivalently, the likelihood-invariant degeneracy group is four-dimensional, generated by (a) the joint product rescalings of $k_1 E_{1T}$ and of $k_2 E_{2T}$ ($\xi$ invariant), (b) $(K_m, S_T) \to (sK_m, sS_T)$ ($\kappa$ invariant), and (c) an independent enzyme-total direction. All absolute molecular scales are therefore structurally unidentifiable from steady-state fitting.

*Proof sketch.* The steady state is a rational function of $(\xi, \kappa)$ alone, so any parameter transformation preserving both combinations leaves the full curve invariant; the four generators above are independent and preserve them, and any transformation moving either combination moves the curve. $\square$

*Numerical verification.* Against a reference curve, rescaling $(K_m, S_T)$ by $\times 5$ gives $\max|\Delta u| = 0.0$ over the full titration range; arbitrary joint $(k, E)$ rescalings give $\max|\Delta u| = 1.7 \times 10^{-18}$ (machine precision).

**Asymmetric Michaelis constants ($\kappa_1 \neq \kappa_2$).** With $\kappa_i = K_{m,i}/S_T$ the steady state is

$$\xi = \frac{u}{1-u}\cdot\frac{\kappa_1 + 1 - u}{\kappa_2 + u}. \tag{S2.2}$$

The parameter space rises to seven dimensions and the identifiable set to the triple $(\xi, \kappa_1, \kappa_2)$, so the group dimension is $7 - 3 = 4$, unchanged from the symmetric case: asymmetry enlarges the identifiable set, not the group. On a 200-point logarithmic titration grid ($\kappa_2 = 0.01$, $K_{m,1} = 0.03$ reference), the four independent group actions are exact ($\max|\Delta u| = 0$, $2.2 \times 10^{-16}$, $0$, $0$ respectively), while a counterexample action moving $K_{m,1}$ alone shifts the curve by $1.45 \times 10^{-1}$, as expected for a non-group direction.

**Closed-form Hill coefficients.** Evaluating the logit slope of Eq. (S2.1) at $u = 1/2$ gives the exact closed forms

$$n_H = 1 + \frac{1}{2\kappa} \tag{S2.3}$$

(symmetry) and

$$n_H = \frac{4}{4 - 1/(\kappa_1 + \tfrac12) - 1/(\kappa_2 + \tfrac12)}, \tag{S2.4}$$

which agree with numerical logit slopes to machine precision (relative deviations $\leq 1.1 \times 10^{-14}$ across $r = \kappa_1/\kappa_2 \in [0.03, 30]$). Asymmetry has a sign effect: at $\kappa_2 = 0.01$, $r = 0.03$ raises $n_H$ to $98.973$ (versus $51.0$ for the symmetric reading), because saturation on the modifying-enzyme side contributes more to ultrasensitivity; $r = 30$ lowers it to $5.068$.

### S2.2 Theorem 2 (Fisher barrier): shape identifiable, scale not

**Theorem 2.** Under the additive Gaussian measurement model $y_{ij} = u(\xi_j; \theta) + \varepsilon_{ij}$, $\varepsilon_{ij} \sim \mathcal{N}(0, \sigma^2)$, with Fisher matrix $F = \sigma^{-2}\sum_j \nabla_\theta u(\xi_j)\,\nabla_\theta u(\xi_j)^{\mathsf T}$ conditioned on the titration design (20 points over 3 decades, $\sigma = 0.03$, Western-blot-grade noise):

(i) the shape parameterisation $\theta = (\ln c, \ln \kappa)$ is well conditioned, $\chi \approx 19$–$34$ for $\kappa \in [0.005, 1]$;

(ii) the absolute parameterisation $\theta = (\ln c, \ln K_m, \ln S_T)$ carries an exact zero eigenvalue along the degeneracy-group orbit ($\chi \sim 10^{300}$ in the deep-saturation extrapolation).

*Proof sketch.* (i) is a direct computation of the $2 \times 2$ shape block. (ii) follows from Theorem 1: the group action leaves every $u(\xi_j)$ invariant, so $F$ annihilates the orbit tangent; the barrier is structural, not a practical difficulty, and no increase in sample size or precision crosses it. $\square$

*Extrapolation note.* In double-precision arithmetic $\lambda_{\min}$ hits the floating-point noise floor near $10^{-14}$ (including small negative values), so the archived code returns effective condition numbers up to $\chi \sim 10^{33}$ at the deepest $\kappa$ values; $\sim 10^{300}$ is explicitly labelled as the deep-saturation extrapolation. Both figures support the identical conclusion, because the barrier criterion is the cross-magnitude stratification between the shape block and the absolute block. A referee reproducing the code should see finite saturated values up to ${\sim}10^{33}$ at the deepest $\kappa$, not $10^{300}$; this is floating-point saturation, not a contradiction. The barrier does not depend on the symmetry assumption: for the asymmetric model (S2.1) the shape block $(\ln c, \ln \kappa_1, \ln \kappa_2)$ stays at $\chi = 37.9$–$207.4$ across $r \in [0.03, 30]$, while the absolute block retains a numerical-zero eigenvalue at every asymmetry ratio.

From Eq. (S2.3), $\kappa$ is recoverable from curve shape exactly: $\kappa = 1/\big(2(n_H - 1)\big)$, with the threshold pinned at $\xi = 1$ and the closed-form slope $d\xi/du\big|_{1/2} = 4 - 2/(\kappa + \tfrac12)$.

### S2.3 Theorem 3 (one-sided protocol for the deep zero-order regime)

**Theorem 3.** For the fixed design above, the 95% confidence width on $\kappa$ inflates as $\kappa$ decreases, and below $\kappa \lesssim 0.01$ only a one-sided upper bound is reportable:

| $\kappa$ | 0.001 | 0.005 | 0.01 | 0.05 | 0.1 | 0.5 | 1.0 |
|---|---|---|---|---|---|---|---|
| 95% CI inflation | $\times 2920$ | $\times 5.3$ | $\times 2.4$ | $\times 1.34$ | $\times 1.29$ | $\times 1.33$ | $\times 1.45$ |
| $n_H$ | 501 | 101 | 51 | 11 | 6 | 2 | 1.5 |

*Protocol construction.* The estimate is replaced by a scan over the curve family: the reported quantity is $\kappa \le \kappa_{\max}$, where $\kappa_{\max}$ is the smallest value such that all curves with $\kappa \le \kappa_{\max}$ are likelihood-equivalent to the best fit within noise over the entire scan grid. Because the bound is a property of the degenerate curve family rather than of any point estimate, it is invariant to the optimizer and its initialisation: the optimizer-invariance criterion that distinguishes a scan-based bound from a spurious point estimate.

*Worked example ($\kappa = 10^{-3}$).* The CI inflation reaches $\times 2920$ and $n_H = 501$; every switch in the deep zero-order region looks the same, so the honest statement is the upper bound, never a point estimate. This is exactly the one-sided reporting protocol already in use for the RIA ($\kappa \ge 1$) and ELISA ($\kappa_3 \ge 18.8$) assay classes.

### S2.4 Orthogonal resolution: substrate and enzyme titration

Substrate titration moves $\kappa = K_m/S_T$ point by point along the $S_T$ axis, while enzyme titration holds $\kappa$ fixed and moves $\xi$; the two titration directions are orthogonal in the $(\xi, \kappa)$ plane. Combining the two titration classes (implementable in an in vitro reconstituted system) converts absolute $K_m$ and $S_T$ from unidentifiable in principle to identifiable in principle; the residual practical barrier is audited by the same Fisher/profile-likelihood machinery of S2.2–S2.3. This is the separation-uniqueness analogue of the assay-side sample–platform orthogonality.

### S2.5 The positive-feedback switch: $g^* = 0.443$

Embedding the cycle in positive feedback, $u = F\big(\xi(1 + g\,u)\big)$ with $F$ the inverse Goldbeter–Koshland steady state (equivalently $\xi(u) = M(u)/(1 + g\,u)$, solved by brentq at xtol of $10^{-14}$), folds the gain–dose curve beyond a critical feedback gain. At the main case $\kappa = 0.05$ (baseline $n_H = 11$ from Eq. (S2.3)) the fold condition locates

$$g^* = 0.443. \tag{S2.5}$$

Beyond $g^*$ the apparent Hill coefficient rises from the baseline $11$ to $\approx 480$ at $g = 10$, and the up/down fold points separate by a hysteresis width of $0.44$–$0.58$ dex at $g = 10$–$20$. The degeneracy group survives the fold: all four group actions of Theorem 1, applied pointwise in both the one-root and the three-root regions, leave every branch of the bistable response invariant with $\max|\Delta u| = 0$. Classification is thus a property of the operating regime, not of molecular identity: one molecular topology, two regimes, one degeneracy group.

### S2.6 Information optima: the zero-order limit is not optimal

Treating the cycle as a channel with log-normal input $\ln \xi$ (median 1, spread $\sigma_\xi$) and noisy output $u$ ($\sigma = 0.03$), the steady-state mutual information $I(\xi; u)$ has an interior optimum in $\kappa$:

| $\kappa$ | 0.002 | 0.01 | 0.05 | 0.1 | 0.3 | 1 | 3 | 10 |
|---|---|---|---|---|---|---|---|---|
| narrow input, $\sigma_\xi = 0.3$ | 1.42 | 2.12 | 2.84 | **2.85** | 2.44 | 1.87 | 1.58 | 1.46 |
| wide input, $\sigma_\xi = 1.5$ | 1.11 | 1.41 | 2.07 | 2.42 | 2.85 | **3.03** | 3.04 | 3.02 |

(bits). The optimum tracks the input dynamic range: $\kappa^* \approx 0.05$–$0.1$ ($n_H \approx 6$–$11$) for narrow input ensembles and $\kappa^* \approx 1$–$3$ ($n_H \approx 1.5$–$2$) for wide ones. The zero-order ultrasensitive limit is information-poor for all input widths (1.1–1.4 bits): per snapshot, the infinite-gain switch is a decision element, not a measurement element. The variational reading is that a covalent cycle's saturation is matched to its input distribution rather than pushed to the zero-order limit.

*Source archive:* `08_audit_trail/01_cell_line_1/manuscript/doc_audit_and_revision_reports/PdPC_identifiability_audit_memo.md`; NC SI v1.5 S1, S3 (`08_audit_trail/01_cell_line_1/manuscript/supplementary_information_SI_v02.md`); `NC_manuscript_v11_fulltext_EN.md` Methods (positive-feedback switch; asymmetric Michaelis, script 10).

---

---

## S3. p53 topology selection

Pulse-generator circuits are partitioned by interaction graph into the NF class (pure negative feedback: Bar–Or^48^, Monk^13^, Geva–Zatorsky Models I–III^5^) and the EXC class (a positive-feedback loop or equivalent fast autocatalysis, with a resting state coexisting with finite-amplitude excursions; Ciliberto-type NF+PF^49^, Mönke et al. FHN reduction^10^). The observables are the pulse amplitude $A(D)$ and the per-damage-event pulse count $N(D)$, with discrimination indices

$$\alpha = \frac{\partial \ln A}{\partial \ln D} \;\;(\text{analogue index}), \qquad \nu = \frac{\partial N}{\partial \ln D} \;\;(\text{digital index}). \tag{S3.1}$$

### S3.1 Lemmas 1–3 and Theorem 4

**Lemma 1 (NF onset is Hopf, $A \propto \sqrt{D - D_c}$).** An NF-class circuit has at most one equilibrium $e^*(D)$ (Thomas' first conjecture, proved by Soulé^8^ and Gouzé^50^: multistationarity requires a positive loop in the interaction graph). *Proof sketch.* (i) Uniqueness is the Thomas–Soulé contrapositive above. (ii) Generic codimension-1 destabilisation of a unique equilibrium is either fold-type (extra equilibria, contradicting (i); transcritical/pitchfork are non-generic) or a complex pair crossing the imaginary axis (Hopf). (iii) The supercritical Hopf normal form gives limit-cycle radius $r \approx \sqrt{\mu/\ell_1}$, $\mu \propto (D - D_c)$, $\ell_1$ the first Lyapunov coefficient, hence

$$A(D) \propto \sqrt{D - D_c}, \qquad \frac{dA}{dD} \sim \frac{1}{\sqrt{D - D_c}} \to \infty \;\; (D \to D_c^+). \tag{S3.2}$$

$\square$ **Corollary 1a.** Every NF circuit possesses a dose interval with $\partial A/\partial D \neq 0$, and its analogue sensitivity diverges at onset. *Exception note (subcritical Hopf).* Strongly nonlinear NF circuits can produce finite-amplitude onset via a limit-cycle fold, but this carries observable hysteresis, amplitude still varies with $D$ above the fold, and full-range $\partial A/\partial D \approx 0$ is unattainable; the digital signature cannot be rescued (Theorem 4, second signature).

**Lemma 2 (EXC amplitude is drive-independent, exact in the singular limit).** For the FitzHugh–Nagumo representative $\dot v = f(v) - w$, $f(v) = v - v^3/3$, $\dot w = \varepsilon\,(v + a - bw - I(D,t))$, with drive entering only the slow equation, the fast nullcline does not move with $D$; in the singular limit $\varepsilon \to 0$ the excursion is fixed by knee geometry (take-off $v_L = -1 \to v_R = 2$; return at $v = -2$), so $A \in \{3, 4\}$ by baseline convention, independent of $I(D,t)$ exactly; finite separation introduces $O(\varepsilon)$ corrections. *General form.* For a fast–slow system $\dot x = f(x,w)$, $\dot w = \varepsilon g(x,w,D)$ whose fast subsystem has attracting branches terminating at folds $w_\pm$: if the drive enters only through the slow equation, then in the singular limit

$$A = x_{\text{fast-branch}}(w_-) - x_{\text{fast-branch}}(w_+), \quad \text{independent of } D \text{ exactly}. \tag{S3.3}$$

If the drive also enters the fast equation, the fast nullcline moves with $D$ and amplitude acquires an $O(1)$ but weak ($|\alpha| \ll 1$) drive dependence: the leaky-digital regime. *Necessity of positive feedback* follows from Lemma 1(i) by the converse use of Thomas' theorem: excursion–rest coexistence requires a separatrix structure, which requires a positive loop. *General-dimension note.* For fast–slow systems of dimension three or higher, Fenichel theory carries the statement: an excursion is slow flow to a fold surface $F$, a fast jump off it, and landing on the other attracting branch; amplitude is the distance between the jump endpoints, the drive shifting only the arrival time at $F$. Premises: two-time-scale structure, fold destabilisation, transversal arrival at $F$ (canards excepted; non-fold fast destabilisations lie outside the domain).

**Lemma 3 (counting law).** With damage decaying by repair, $D(t) = D_0 e^{-t/\tau_r}$, and repair slower than the pulse period ($\tau_r \gg T$), pulses fire at period $T$ within the supra-threshold window, giving logarithmic counting with dynamic range $\sim \ln D_0$ rather than $D_0$:

$$N(D_0) \approx \frac{t_{\text{above}}}{T} = \frac{\tau_r}{T}\ln\frac{D_0}{D_c}. \tag{S3.4}$$

**Theorem 4 (topology selection).** Given $(A(D), N(D))$ measured across a dose range at $\le 15\%$ amplitude precision: (1) if $\alpha \lesssim 0.1$ and $\nu > 0$ (amplitude invariant, count growing), the circuit must contain a positive-feedback loop or excitable structure (EXC class); (2) if $A^2$ is linear in $D$ with $\alpha \gtrsim 0.5$ ($\sqrt{}$-onset, amplitude drifting with dose throughout), behaviour is consistent with the pure-NF class; (3) the classes are separated by an order of magnitude in the $(\alpha, \nu)$ plane (NF: $\alpha \sim 0.5$–$\infty$; EXC: $\alpha \sim 0.01$–$0.1$), so the discrimination is robust to measurement noise. *Proof structure.* (1) by Lemma 1 (NF cannot be digital) and Lemma 2 (EXC must be digital) in a squeeze; (2) by the Lemma 1 scaling law; (3) by the numerics of S3.2. The operative discrimination is the three-signature triple (onset continuity; hysteresis width under up/down scans; range $\alpha$); no single signature suffices, and the subcritical escape hatch is caught independently by the second and third.

### S3.2 Continuation-grade numerics

The NF representative (three-ODE Goodwin-type, analytic equilibrium branch) has its unique destabilisation at $p^* = 1.477$, $s_H = 16.4449$, with crossing pair $-0.000 \pm 0.924i$ (Hopf, corroborating Lemma 1(i)–(ii)); the cycle above is a stable small-amplitude branch (supercritical) with $A^2$ linear in $(s - s_H)$ at $R^2 = 0.9926$ over $s - s_H \in [0.05, 3.80]$. The delayed-NF cross-check ($\tau = 8$) gives $R^2 = 0.971$, $s_c = 0.0406$, and a $3.6\times$ amplitude change over a $2.5\times$ drive range. The subcritical normal form $\dot z = (\mu + i)z + |z|^2 z - |z|^4 z$ shows a finite-amplitude onset ($r = 1.00$) that mimics a digital signature but is caught by both remaining signatures: hysteresis width $0.275$ against the normal-form prediction $0.25$, and amplitude still varying above onset ($r: 1.00 \to 1.13$ over $\mu: 0 \to 0.35$). The EXC representative (FHN, slow drive) holds amplitude to $6.5\%$ over a $4.6$-fold drive range ($\alpha = 0.014$, against $360\%$ for NF at matched range). The Morris–Lecar Type-II control (drive entering the fast equation) gives $A = 83.7 \to 68.4$ mV over a $2.1$-fold drive range, $\alpha_{ML} = -0.258$: the leaky-digital signature, an order of magnitude below the NF analogue signature, confirming Lemma 2's drive-position classification on a second system. The counting law verifies at $R^2 = 0.920$ with $D_c = 0.261$, and its slope is a free period estimator: back-computed $T = 2.14$ au against the measured $2.30$ au.

### S3.3 The Batchelor et al.^9^ built-in control and the sigma reconciliation

Batchelor et al.^9^ provide a control built into a single study: the same protein in the same cell line (MCF7) shows excitable, fixed-amplitude, dose-counting pulses after double-strand breaks but a graded, non-excitable single pulse after UV. Digitisation of Fig. 1G/H (colour-channel bar segmentation, visual axis calibration; error bars $\pm 14$–$30$ AU) gives NCS (100/200/400 ng ml⁻¹) amplitudes $137.9/193.1/137.9$ AU, non-monotonic, a $1.40\times$ spread, so amplitude carries no monotone dose information; UV (2/4/6/8/10 J m⁻²) amplitudes $84.8/103.3/181.5/366.3/331.5$ AU, monotone $4.3\times$ over 2–8 J m⁻² ($\sim 14$ error bars). Durations concord (NCS exponent $-0.21$; UV $127.7 \to 461.0$ min, exponent $+0.79$). Bootstrap resampling ($\times 20000$ within the bar-top extraction noise model) gives

$$\alpha_{NCS} = -0.001 \pm 0.113 \;\; (95\%\ \mathrm{CI}\ [-0.23, 0.22]), \qquad \alpha_{UV} = +0.981 \pm 0.114 \;\; (95\%\ \mathrm{CI}\ [0.78, 1.22]). \tag{S3.5}$$

**Sigma reconciliation (frozen cards versus main text).** The archived frozen verdict cards quote the two-arm separation as $\approx 8.6$ joint standard deviations. That figure uses the cards' pre-registered single-sided denominator: the slope difference $\Delta\alpha = \alpha_{UV} - \alpha_{NCS} = 0.982$ divided by the larger single-arm bootstrap error, $0.982/0.114 \approx 8.6$. The main text instead uses the standard independent-error combination $\sqrt{0.113^2 + 0.114^2} \approx 0.1605$, giving $0.982/0.1605 \approx 6.1$ combined standard deviations. Both statements are correct under their stated denominator conventions: a referee reproducing the single-sided convention lands on $8.6\sigma$; a referee using the independent-error combination lands on $6.1\sigma$. The frozen cards are retained verbatim per archive discipline; this paragraph is the reconciliation note, and no conclusion depends on the choice (the EXC/analogue regions are separated by an order of magnitude in $\alpha$).

**By-product: the counting-law slope inverts the repair timescale.** Lemma 3's slope gives $\tau_r \approx \nu \cdot T = 1.4 \times 5.5\ \mathrm{h} \approx 7.7\ \mathrm{h}$, consistent with the fast component of DSB repair (hour scale); the slope is thereby promoted from fit parameter to measurand. *Script-number reconciliation.* NC SI v1.5's logical "script 17" is the archived `code52_TrackC2` (`04_code/code52_TrackC2_discrimination_diagram.py`; data `code52_TrackC2_results.json`), per the frozen mapping table of SI S0.

### S3.4 TrackC2: empirical validation on the discrimination diagram

TrackC2 (script `code52`, verdict card `verdict_card_TrackC2_discrimination_diagram_2026-09-18.md`, closed 2026-09-18 at text-plus-digitised level) places the measured $(\alpha, \nu)$ signatures on the discrimination diagram:

| data | $\alpha = \partial\ln A/\partial\ln D$ | $\nu = \partial N/\partial\ln D$ | region | source tier |
|---|---|---|---|---|
| NCS (DSB) | $-0.001 \pm 0.113$ | $2.0 \pm 0.5$ | EXC | $\alpha$: Fig. 1G digitisation, bootstrap $\times 20000$; $\nu$: Mönke et al.^10^ text |
| $\gamma$-IR (DSB) | $0.05 \pm 0.10$ | $1.4 \pm 0.4$ | EXC | Batchelor et al.^9^ text + Lahav et al.^4^ $N: 1 \to \sim 6$, 0.1–10 Gy |
| UV | $+0.981 \pm 0.114$ | $0.0$ (+0.2) | analogue/NF | $\alpha$: Fig. 1H digitisation, bootstrap $\times 20000$; $\nu$: the paper's "single pulse" |
| model anchor NF | 1.13 | $\approx 0$ | n.a. | archived TrackB value |
| model anchor EXC | 0.04 | $\tau_r/T$ (parameter) | n.a. | script 51 / memorandum |

Both DSB arms land deep in the EXC region and the UV arm in the analogue region, co-varying with the study's own excitability classification: the signature–topology correspondence of Theorem 4 holds on real data. Verbatim anchors: "DSBs trigger a series of p53 pulses with fixed amplitude and duration, independent of the damage dose, whereas the number of pulses increases with higher damage"; the UV arm is a "single pulse that increases in amplitude and duration in proportion to the UV dose." Honest boundary: the closure is text-level plus digitised-level, not pixel-level re-digitisation; independent error bars for $\nu$ require WebPlotDigitizer extraction of Lahav et al.^4^ and Mönke et al.^10^ (registered upgrade; tightens uncertainties, moves no placement).

### S3.5 Hysteresis-protocol design (two-way dose scan)

The second signature (hysteresis width) has never been measured directly for p53: all published dose series are one-way scans. The protocol converts the strongest attack surface of Theorem 4 (the subcritical escape hatch) into a falsifiable experiment: EXC predicts zero hysteresis; any nonzero hysteresis forces reclassification to the subcritical column, with "excluding pure negative feedback" downgraded to "excluding supercritical NF" (the conditional-strength clause). *System.* MCF7 p53-CFP / Mdm2-YFP single cells (the Lahav et al.^4^ readout, for comparability with all existing placements); live-cell time-lapse, 20 min frames, $\ge 16$ h. *Dose ladder (two-way scan).* NCS $25 \to 50 \to 100 \to 200 \to 400 \to 200 \to 100 \to 50 \to 25$ ng ml⁻¹, dwelling 4–6 h per step (2–3 pulse periods), medium exchange on the descending limb; or an equivalent stepped $\gamma$-irradiation design ($0.5 \to 1 \to 2.5 \to 5 \to 10 \to 5 \to 2.5 \to 1 \to 0.5$ Gy, cumulative), with a one-way ascending control arm on the same cell population. *Frozen readout criteria.* Onset = first-pulse dose on the up-scan; extinction = last-pulse dose on the down-scan; width $w = |\mathrm{onset} - \mathrm{extinction}|$/active window; width below $5\%$ is EXC-consistent (the prediction), width $\ge 5\%$ is a subcritical signature triggering reclassification and the downgrade clause. *Power.* At $\le 15\%$ amplitude precision, $\alpha$ is estimated to $\approx \sigma_A/\Delta\ln D \approx 0.075$, below the NF/EXC boundary gap of $0.2$; $n = 30$–$50$ cells per arm suffice, and the ladder density (3–4 points per decade) bounds the width resolution to half a step. *Registered risks.* (i) NCS is an irreversible damaging agent, so descending-limb "extinction" may be damage-accumulation death, countered by pulsed dosing with washout and censoring of dying cells; (ii) stepped $\gamma$-irradiation confounds cumulative with instantaneous dose: the two arms are mutual controls, both dose conventions reported; (iii) cell-cycle-position confounding: re-test on synchronised controls.

*Source archive:* `08_audit_trail/01_cell_line_1/manuscript/doc_audit_and_revision_reports/topology_selection_theorem_pulse_observables_circuit_topology.md`; `TrackB_theorem_rigorisation_report.md`; `TrackC_experimental_data_placement_report.md`; `08_audit_trail/02_cell_line_2/results/reconciliation_p53_Batchelor2011_two_arm_structure_v01.md`; `08_audit_trail/03_cell_line_3/results/verdict_card_TrackC2_discrimination_diagram_2026-09-18.md`; NC SI v1.5 S5.1, S10, S11, S12.

---

---

## S4. p53 channel census and dispersion signatures

### S4.1 Fisher propagation protocol

The audit object is the EXC-class representative, a FitzHugh–Nagumo reduction with exponential damage repair,

$$\dot v = s\,(v - v^3/3 - w) + \sigma\,\xi(t), \qquad \dot w = s\,\varepsilon\,(v + a - b\,w - I(t)), \qquad I(t) = I_0 e^{-t/\tau_r}, \tag{S4.1}$$

with locked parameter set $s = 20$, $\tau_r = 20$, $a = 0.7$, $b = 0.8$, $\varepsilon = 0.08$, $I_0 = 0.6$, $\sigma = 0.06$. FHN is a universality-class representative, not a p53 molecular model; all conclusions carry the class-level identifiability reading. Six pulse statistics (first-pulse time $t_1$, mean period $T$, amplitude $A$, count $N$, pulse width $w$, and inter-pulse-interval dispersion $\mathrm{IPI}_{sd}$) are propagated onto the seven dynamical parameters $(a, b, \varepsilon, s, \tau_r, I_0, \sigma)$. Monte Carlo ensembles ($M = 400$ or $4000$ trajectories) give the statistic covariance $\Sigma_g$; per-parameter log-perturbations ($\pm 2\%$ and $\pm 5\%$ double runs; $\pm 8\%$ for the discrete $N$ channel) give the Jacobian by finite differences; the Fisher matrix is $F = J^{\mathsf T} \Sigma_g^{+} J$ (pseudoinverse, truncation $10^{-12}$). Three numerical configurations were run for the spectral claims: $M = 400$ central differences, $M = 4000$ forward differences, and $M = 4000$ central differences. *Double-recorded failures.* Round r1 died on an off-by-one in the pulse-width extraction (zero mean and variance), which the pseudoinverse propagated into three spurious "zero directions" (bug fixed and recorded). Round r2 revived $w$ (share 57.7%) but a semi-implicit integrator mis-step (using the updated $v$ in the $w$ equation) drifted $N$ to $8.3$; restoring synchronous integration returned $N = 9.0$. The original pre-registered proposition P2 ($\mathrm{Var}(t_1)/\mathrm{Var}(T) \ge 4$ for the bare model) failed across all regimes (0.10–0.28) and is on record as falsified; its diagnosis (a missing upstream stochastic stage) produced the re-pinned proposition P2′, fixed before the r3 run. Both failures and the revision trail are preserved in the verdict card.

### S4.2 Channel shares and the gain-dispersion kill

The Fisher information of the seven dynamical parameters distributes over the six channels as:

| channel | $w$ width | $N$ count | $T$ period | $A$ amplitude | $t_1$ first pulse | $\mathrm{IPI}_{sd}$ |
|---|---|---|---|---|---|---|
| share | 60.2% | 18.2% | 11.5% | **5.5%** | 4.4% | 0.3% |

Timing-type channels carry 94.2% of the total parameter information; amplitude carries 5.5% (rounding note: the unrounded timing share is 94.196%, reported as 94.2%; the rounded per-channel shares sum to 94.3%). This is the information-theoretic quantification of Lemma 2: amplitude is not merely dose-invariant, it carries almost no information about any dynamical parameter. The second layer is statistical: with per-cell gain dispersion $k_i \sim \mathrm{LogN}(0, \sigma_k)$, growing $\sigma_k$ from $0$ to $0.6$ drives $\mathrm{CV}(A)$ from $0.002$ to $0.639$ while $\mathrm{CV}(T)$ stays pinned at $0.016$. The amplitude channel's two-layer death is thereby closed: geometric at the single-cell layer (knee-geometry lock, 5.5% share) and statistical at the population layer (uncalibratable gain dispersion), while timing channels survive both layers.

### S4.3 Moment-hierarchy separation: noise is a null direction

The softest Fisher direction is the noise amplitude $\sigma$ itself, with eigenvector $\sigma$-weight $0.996$–$0.999$ across the three numerical configurations of S4.1. Under the $M = 4000$ central-difference configuration,

$$\lambda_7/\lambda_1 \;=\; 4.7 \times 10^{-19} \quad (\text{read as the upper bound } \lambda_7/\lambda_1 \le 4.7 \times 10^{-19}), \tag{S4.2}$$

below the floor of double precision. Because this direction is already pressed against the finite-difference resolution floor, only the upper bound is reportable; an exact condition number would require $M \gg 4000$ or analytic sensitivities, and this limitation is registered as such. The content: first-moment (mean-type) pulse statistics are exactly blind to the fluctuation-generating noise; information about $\sigma$ lives exclusively in second-moment, dispersion-type observables.

### S4.4 Analysis step 3 closure: zero direction and Cramér–Rao attainment (script 13)

The dynamic-encoder analogue of the static Fisher barrier closes audit step 3. With counting law $N(D) = \rho \ln(D_0(D)/D_c)$, $D_0 = 0.40 + 0.09 D$, per-dose counts $n \sim \mathrm{Poisson}(N)$ over $M = 200$ cells, amplitude observable $\bar y \sim \mathcal{N}(sA, \sigma^2)$ with $\sigma = 15\%\, sA$, and parameter vector $\theta = (\ln D_c, \ln \rho, \ln A, \ln s)$ (structural parameters $D_c, \rho$; scale nuisance $s$): (i) in the counting-only digital limit the Fisher spectrum is $\lambda = [26935, 1898, 0, 0]$, exactly two nonzero directions; (ii) adding the amplitude observable adds exactly one nonzero direction (the product $sA$), $\lambda = [106667, 26935, 1898, 0]$, and the zero-direction eigenvector $(0,0,1,-1)/\sqrt 2$ coincides with the degeneracy-group orbit tangent at cosine $1.000000$; (iii) Monte Carlo estimation ($R = 400$, Poisson MLE) against the Cramér–Rao bound gives ratios $1.24$ ($\ln \rho$), $1.36$ ($\ln D_c$), and $0.95$ for the estimable combination $\ln(sA)$, while fixing $s$ at a wrong value leaves the likelihood invariant under the compensation $\hat A = sA/s$: the whole group orbit is degenerate. The structural block $(\ln D_c, \ln \rho)$ has condition number $14.2$ (same order as the PdPC shape block, $\sim 30$); the scale direction is an exact zero ($\lambda_{\min}/\lambda_{\max} = 0.00\mathrm{e}{+00}$, not a floating-point approximation). The dynamic encoder's Fisher barrier is isomorphic to the static sensor's: identifiable counting structure against unidentifiable amplitude scale. All three predictions of the conditional digital-limit theorem are verified numerically on this circuit.

### S4.5 The sensing-stage signature

For any bare oscillator the dispersion ratio is bounded at $O(1)$: $\mathrm{Var}(t_1)/\mathrm{Var}(T) = 0.03$–$0.12$, because without a delay stage the first-pulse variance is constrained by single-pulse trigger noise to the IPI-variance scale. The literature value $5.8$ (Lahav et al.^4^: $t_1$ dispersion $\pm 240$ min versus IPI $\pm 100$ min; measured on MCF7 single cells at 5 Gy $\gamma$-irradiation, the 146-cell variability set of Geva-Zatorsky et al.^5^ Fig. 4 and its Lahav et al.^4^ predecessor, with between-cell period dispersion $5.5 \pm 1.5$ h and oscillating fractions $\approx$ 50%/90% at 5/10 Gy) is therefore physically impossible for a bare oscillator. Adding an upstream exponential stochastic delay $\tau_d \sim \mathrm{Exp}(\mu_d)$, standing in for the ATM damage-sensing stage, flips the ratio monotonically across four orders of magnitude: $\mu_d/T = 0.10 \to 64$; $\mu_d/T = 0.73$ (the literature calibration $240/330$ min) $\to 1953$. Under realistic noise calibration ($\sigma^* = 0.25$, giving intra-cellular IPI CV of $0.316$, approximately the literature $0.30$), $\mu_d \approx 1.0 \cdot T$ reproduces the ratio $5.7 \approx$ the measured $5.8$ (a quadratic extrapolation from the calibrated $0.73 \to 3.01$ segment, not a simulated point). Verdict: the $t_1$/IPI dispersion ratio is a sufficient signature of an independent stochastic sensing stage upstream of the oscillator core, and its magnitude inverts the upstream delay dispersion, the third row of the topology discrimination table. The mechanistic prior art (Mönke et al.^10^: stochastic DSB repair plus a Type-I excitable system reproducing pulse-number and IPI heterogeneity) is thereby theorematised into a necessary signature with a quantitative inversion; the present model does not claim a molecular identity for the delay stage.

### S4.6 The noise-colour plane

Slow parameter heterogeneity (per-cell $\pm 5$–$20\%$ variation in $a$, $\varepsilon$, $s$, $\tau_r$) inflates only the population $\mathrm{CV}(T)$ (up to $0.20$) while the intra-cellular IPI$_\mathrm{CV}$ stays pinned at $0.083$; fast white noise $\sigma$ moves both axes together ($\sigma = 0.25$: IPI$_\mathrm{CV} = 0.32$, with $N$ shifting only $+16\%$). The (population $\mathrm{CV}(T)$, intra-cellular IPI$_\mathrm{CV}$) plane thus separates fluctuation sources by time colour, and the measured IPI$_\mathrm{CV} = 0.30$ can arise only from fast intra-cellular stochasticity. The Geva-Zatorsky et al.^5^ reading (slow noise explaining inter-cellular variability) is refined, not contradicted: intra-cellular irregularity must be fast, while population period dispersion may be slow. This is the fourth row of the discrimination table.

### S4.7 Counting-channel capacity

The mutual information between dose and pulse count is

$$I(D; N) = 3.07 \;/\; 3.17\ \text{bits} \;\; (96.8\%) \tag{S4.3}$$

over a four-fold dose range ($I_0: 0.30 \to 1.20$), with the counting law $N \approx 10.5 \ln I_0 + c$ linear throughout. The counting channel is near-lossless: p53 counts not only because counting is degeneracy-immune, but because the digital channel's efficiency ceiling is essentially attained. (Binned mutual-information estimate following the script-50 convention: trends reliable, absolute values mildly biased.)

*Source archive:* `08_audit_trail/03_cell_line_3/results/` scripts `code51` r1–r4 (r3 final; data `code51r3_p53_Fisher_audit_results.json`, `code51r4A_noise_calibration_robustness.json`, `code51r4B_Fisher_M4000.json`); `verdict_card_P53_Fisher_audit_2026-09-18.md`; `recon_card_P53_Fisher_audit_2026-09-18.md`; `closure_card_p53_signalling_line_2026-09-18.md`; NC SI v1.5 S7 (script 13), S9; `08_audit_trail/01_cell_line_1` script 13.

---

## S5. The digital limit: theorems, conjecture, boundaries

This section condenses the theorem-level material behind Results section 4 from the NC-lineage SI (v1.5, section S2) and the associated NC-lineage working note internally numbered S27 (an archived programme document, not a section of the present Supplementary Information; its verification scripts carry the same prefix); open items are labelled as open.

### S5.1 The conditional theorem

**Setting and statement.** Parameters $\theta \in \Theta \subset \mathbb{R}_+^d$ (positive rates, concentrations, Michaelis constants); observables depend only on dimensionless monomial ratios (Buckingham $\pi$); $G$ is the connected likelihood-preserving transformation group. The output is an event stream, each event carrying an amplitude $y \in \mathbb{R}_+$ and a bin $k \in \{1, \dots, K\}$ under scale-free binning. Hypotheses: H1 (multiplicative degeneracy), $G$ acts diagonally, $g \cdot \theta = (c_1\theta_1, \dots, c_d\theta_d)$; H2 (digital limit), conditioned on the bin, $p(y \mid k; \theta)$ depends on $\theta$ only along $G$ orbits and is constant on $\Theta/G$; H3 (event independence), Bernoulli/Poisson-type streams. Conclusions: C1, $\dim G \ge 1$ implies $G$ contains a one-dimensional generalised-scale subgroup; C2, the occupancy vector $n = (n_1, \dots, n_K)$ is sufficient for $\theta$; C3, any partition statistic invariant under $G$ and likelihood-relevant is a function of $n$.

**Proof sketch** (Lemmas DL1–DL3 here; DL4 in S5.3). DL1 (C1): under the log map, the ratio-preserving diagonal group is a linear subspace acting by translation, whose one-parameter subgroups map back to generalised scalings $\theta_i \mapsto \theta_i t^{a_i}$; H1 is defended by Buckingham $\pi$ (both audited circuits have explicitly multiplicative groups), with additive baseline-shift degeneracies explicitly out of scope. DL2: under H2–H3 the occupancy vector is multinomial given $N$, or independent Poisson, with orbit-constant probabilities; verified constructively in both circuits (PdPC at $\kappa \to 0$: independent Bernoulli molecules, hence binomial; p53: Poisson count jitter after scale cancellation). DL3 (C2): Fisher–Neyman applied to $p(n;\theta) = [N!/\prod_k n_k!]\prod_k q_k(\theta)^{n_k}$.

**CV counterexample dissolved.** CV is a scale-invariant partition statistic not of counting type, and the earlier v0.8 statement was open to it; but in the digital limit the amplitude distribution is constant on $\Theta/G$ (H2), so CV is not likelihood-relevant and is excluded by the tightened C3, which quantifies over likelihood-relevant statistics only.

**Constructive verification (p53, code 13).** Counts only: Fisher spectrum $\lambda = [26935, 1898, 0, 0]$, exactly two nonzero structural directions. With the amplitude observable: $\lambda = [106667, 26935, 1898, 0]$, the zero-direction eigenvector $(0,0,1,-1)/\sqrt{2}$ having cosine $1.000000$ with the group orbit tangent. Structural parameters attain the Cramér–Rao bound (Monte Carlo versus CRB ratios 1.24 for $\ln\rho$, 1.36 for $\ln D_c$, 0.95 for $\ln(sA)$; $R = 400$); fixing the scale $s$ at a wrong value leaves the likelihood invariant after $\hat{A} = sA/s$. The PdPC side is the static analogue (SI S2).

### S5.2 The unconditional generalisation

Unconditional uniqueness is restated as a stratified classification over degeneracy groups, a different surviving class per row, the table itself pre-registered as the honest fallback should the strongest row fail. $G_0$, diagonal multiplicative scaling ($y \mapsto ay$): the proven layer of S5.1, surviving class counts. $G_1 = G_0 \ltimes$ additive baseline shifts ($y \mapsto ay + b$): ratio statistics (fold-change, ratiometric) die here, and the log-linearisation of Lemma DL1 stops. $G_2$, all strictly increasing continuous amplitude reparametrisations ($y(t) \mapsto m(y(t))$): complete amplitude uncalibratability formalised, surviving class event-time/ordinal functionals. On the regular class $\mathcal{Y}_0$ (continuous signals with finitely many strictly monotone segments separated by strict extrema), Lemma L1′ (proven) identifies the maximal invariant of $G_2$ as the full ordinal structure $R_y = \{(s,t) : y(s) \le y(t)\}$, not the extremum skeleton of the earlier sketch (recorded correction: the skeleton suffices only on an excursion-disjoint subclass). The finite-sample version is classical (ranks as maximal invariant, Lehmann–Romano); event-time functionals are functionals of $R_y$; the amplitude direction is unrecoverable from $R_y$, the group-theoretic statement that amplitude must be paid for.

**Theorem T′ ($G_2$ layer, noiseless, regular class).** Assume H2′ (event resolvability, $N \cdot w \lesssim T_{\mathrm{win}}$) and H3′ (amplitude asymptotically parameter-free, $P_\theta(y) = Q_\theta(R_y)\, H(A_y \mid R_y)$ with $H$ independent of $\theta$; this formalises the $G_2$ premise rather than adding physics, since parameter information in the amplitude distribution would constitute a calibration channel contradicting it). Then: (i) $R_y$ is sufficient for $\theta$ (H3′ factorisation plus Fisher–Neyman); (ii) **universal functional theorem**: the minimal sufficient statistic of any H3′-compliant subfamily is a functional of $R_y$ (Bahadur^51^); (iii) where the minimal sufficient statistic lands inside $\sigma(R_y)$ is a per-family classification question; (iv) the amplitude axis contributes zero information about $\theta$.

Two recorded corrections accompany T′. The strengthening implicit in the earlier sketch, "$R_y$ itself is minimally sufficient", was falsified on the Gamma renewal family, and hypothesis H4 (model richness) was retired because it served only that clause. The verification table includes its negative rows: for the parametric Gamma$(k,\lambda)$ renewal family, $(N, \sum X_i, \prod X_i)$ is sufficient (exponential family), a strict coarsening of $R_y$ (pin: datasets identical in $(N,\sum X,\prod X)$ but differently ordered have bitwise-identical likelihood ratios; a 1% change in $\sum X$ shifts the ratio by $\approx 1.9\times10^{-2}$); for the free-shape renewal family, the interval multiset is minimal sufficient, iid exchangeability making the inter-event order likelihood-free (permutation leaves the ratio bitwise invariant to a $10^{-16}$ floating-point tail; moving one interval by 0.1% changes it by $7.5\times10^{-4}$). Both failure modes are lossless compressions inside the surviving class: the surviving-class clause is intact and only the strengthening died. A negative control (pulse amplitude $\equiv c\theta$) violates H3′ and kills the conclusion as expected, confirming H3′ is load-bearing. Gamma pins: `numeric_check_S27_Gamma_family_V3.py`.

**Theorem N (noise robustness, pairwise-comparison distortion).** Let $y_i = x_i + \varepsilon_i$ with $\varepsilon_i$ coordinate-wise independent, and $D_n$ the fraction of pairs whose order flips. N1 (universal concentration): changing one $\varepsilon_i$ flips at most $n-1$ pairs, so bounded differences $c_i = 2/n$, $\sum c_i^2 = 4/n$, and McDiarmid give

$$P\big(|D_n - E D_n| \ge t\big) \le 2\exp(-n t^2/2), \tag{S5.1}$$

needing only coordinate independence (no Gaussianity, no margin assumption): concentration is free, the tool textbook-level, no novelty claimed. N2 (expectation is a margin functional): for $\varepsilon_i \sim N(0,\sigma^2)$ iid,

$$E[D_n] = \int \Phi\big(-\delta/(\sqrt{2}\sigma)\big)\, d\mu_n(\delta), \qquad E[D_n] \le \tfrac{1}{2}\,\mu_n(\delta < \sqrt{2}\sigma z) + \Phi(-z)\ \ (z > 0), \tag{S5.2}$$

with quantile corollary $E[D_n] \le p/2 + \Phi(-\eta/(\sqrt{2}\sigma))$ when a fraction $\ge 1-p$ of pairs have margin $\ge \eta$: the distortion level is set solely by margin mass within $O(\sigma)$ of zero. N3 (atom wall is an information limit, not a technical gap): an atom of mass $\rho$ at zero margin (plateaus, exact ties) makes $\sim\rho^2$ of pairs tied, distorted under any fixed reference (one half under a random-tie convention), a floor $\Theta(\rho^2)$ independent of method: noise kills order exactly where no order exists.

Numeric pins (`numeric_check_S27_noise_concentration_bound.py`): plateau-free signal ($n = 300$, $\sigma = 0.05$): theory 0.053671 versus empirical 0.053472, empirical std $2.4\times10^{-3}$, consistent with (S5.1); 40% plateau: non-tied pairs 0.0608 versus 0.0610, tied-pair flip rate identically 1.0000, total $(1-\rho^2)\cdot0.0608 + \rho^2 = 0.1840$ versus empirical 0.1843; p53 anchor ($\sigma = 15\%$ of amplitude): pulse-baseline flip probability $\Phi(-4.71) \approx 1.2\times10^{-6}$, so cross-magnitude comparisons are essentially indestructible at experimental noise and the fragile region is exactly the tied pairs.

**Open items (honest list).** (i) Measure-theoretic fine print G-a (dominating measure, measurable selection, factorisation on $\mathcal{Y}_0$), registered as routine but not verified line by line. (ii) The parallel event-time jitter formulation: the continuous-jitter bound follows from the same margin analysis with time margins, a one-sentence generalisation not yet written in. (iii) The noise bound for L1 in general functional norms remains open (the pairwise metric is covered by N1–N3). (iv) The physical-realism boundary of the $G_2$ premise stands ($G_2$ is the limiting narrative of complete amplitude uncalibratability); time-axis reparametrisation is a separate open item. We note explicitly that concentration was at one stage registered as requiring tail control of the margin distribution; N1 dissolved that requirement, and margin dependence survives only in the level statement N2 and the atom wall N3.

### S5.3 Lemma DL4 (likelihood-ratio partitions are level sets of counts)

**Statement and proof.** A statistic $T(\mathrm{data})$ is **likelihood-relevant** if for every pair $(\theta_1, \theta_2)$ the likelihood ratio $\Lambda = p(\mathrm{data};\theta_1)/p(\mathrm{data};\theta_2)$ factorises through $T$; this is the precise form of "carrying parameter information". The multinomial likelihood ratio $\Lambda(n;\theta_1,\theta_2) = \prod_k [q_k(\theta_1)/q_k(\theta_2)]^{n_k}$ factorises exactly through $n$; when $\{q_k\}$ generically separates parameters on $\Theta$, the finest common refinement of the likelihood-ratio partitions is the family of level sets of $n$, so every likelihood-relevant statistic is a function of $n$ ($n$ is minimal sufficient). On the invariance side, $p(n; g\cdot\theta) = p(n;\theta)$ by Lemma DL2, so every counting-type statistic, including count ratios such as the dPCR $\lambda$ ratio, is $G$-invariant. Combining: a partition statistic invariant under $G$ and likelihood-relevant is a function of $n$, hence of counting type. $\square$

### S5.4 The conjecture

**Statement.** Let $G$ be the degeneracy group of a channel's observable statistics, and let "unresolved parameter uncertainty" mean uncertainty along the orbits of $G$. The conjecture: if the likelihood-relevant partition statistics of a channel's output are invariant under $G$, those statistics must be of counting type.

**Non-tautology analysis.** Partitioning does not restrict legitimate statistics to counts: within-bin amplitude moments are equally partition statistics, and the conjecture's entire content is that they must either transform under the scale subgroup or carry no parameter information. Scale-invariant non-informative statistics are not counterexamples: CV is scale-invariant and not of counting type, but in the digital limit it carries no dose information and falls outside the likelihood-relevant statement. Fold-change and ratiometric measurement are likewise not counterexamples: they read ratios of amplitude coordinates, cancel only the diagonal common-mode direction, and remain exposed to non-diagonal degenerate directions such as baseline offsets; they are instances of the ratio principle, of which the digital limit is the extreme case. On the existence side we remain honest: constructive proofs cover two circuits, and a third-system constructive proof is the most direct existence test. A first attempt on a third system (Morris–Lecar model, August 2026) passed two of the three predictions but failed the zero-direction–dose-axis alignment check; the failure was attributed to a Fisher-construction defect and adjudicated as construction-unresolved rather than a counterexample, and the item remains open. The originating script predates the audit-archive convention and is not part of the deposit.

**Uniqueness domain.** Single channel, single readout, memoryless, passive sensing. Multi-time-point ratios (fold-change, Weber's law), spatial-gradient ratios, and dual-channel ratiometric schemes are established calibration-free analogue strategies, but all pay extra resources (temporal memory, multi-point sampling, a second channel, active perturbation), all cancel mainly the common-mode gain direction, and a fold-change of counts is still counting.

**Boundary against representational measurement theory.** The classification of scale types by permissible transformation groups (ordinal = strictly increasing, interval = affine, ratio = multiplicative) was completed at the representational layer by Krantz, Luce, Suppes and Tversky^52^ (Foundations of Measurement, 1971–1990). Our stratification $G_0 \subset G_1 \subset G_2$ is consistent with that classification, but the contribution here sits at the inferential layer (sufficiency filtering, noise bounds, audit instruments): the conjecture constrains which output statistics can carry parameter information under degeneracy and claims nothing about scale-type uniqueness in the representational sense.

### S5.5 Rate-independent chemical computation as the third-level instance

The same limit appears in computation theory: rate-independent chemical computation, computation robust to kinetic-parameter uncertainty, is exactly the class of piecewise rational linear functions (Chen, Doty, Reeves and Soloveichik^53^), and becomes Turing-universal with inhibition. Measurement (the assay-side digital limit, derived from calibration-freeness alone), encoding (the cell-side limit, reached constructively on the two circuits), and computation (rate independence) each possess a digital limit in which parameter ignorance is neutralised, and the three derivations are mutually independent. The audit identifies the shared structure as cancellation of the analogue scale factor from the likelihood; the conjecture of S5.4 claims this equivalence is universal for biological channels under unresolved parameter uncertainty.

### S5.6 The cost of going digital

This subsection condenses the energy-cost memorandum, whose stated status is a perspective-level organising framework, not a theorem: constants depend on modelling assumptions; scaling structure is stable. A channel built from components with static relative mismatch $\varepsilon_0$ (neuromorphic CV $\approx 20\%$; cellular protein-expression variation $\sim 30\%$) faces three distinct error sources: thermal noise (capacitance law $E = kT/(2\varepsilon^2)$, each added bit quadrupling energy), gain mismatch (static and systematic, averaging only across components or by calibration), and threshold mismatch.

**Analogue channel.** Gain mismatch imposes an uncalibrated hard ceiling $D \le 1/\varepsilon_0$ resolvable levels, $C_A \le \log_2(1/\varepsilon_0)$ bits ($\varepsilon_0 = 0.2 \to 2.32$ bits; $0.3 \to 1.74$; $0.05 \to 4.32$), which no thermal-noise energy can break. Below the ceiling one component suffices at cost $kTD^2/2$; above it, redundancy averaging over $M = 4\varepsilon_0^2 D^2$ components at optimal allocation gives $E_A(D) = 2kTD^2$ (mismatch costs a factor of 4). Across a depth-$L$ cascade with independent per-stage errors, $\varepsilon_{\mathrm{tot}} = \sqrt{L}\,\varepsilon_1$ forces $\varepsilon_1 = 1/(\sqrt{L}D)$, so

$$E_A(L,D) = 2kT\,L^2 D^2 \quad (\text{correlated worst case } \propto L^3 D^2), \tag{S5.3}$$

consistent with observed error accumulation in deep analogue optical networks (Hamerly et al.^54^).

**Counting channel.** $N$ all-or-none events with Poisson fluctuation $\sqrt{N}$ give $D = \sqrt{N}$ levels, hence $N = D^2$; per-event reliability is barrier-set, $\varepsilon_{ev} = e^{-\Delta G/kT}$, so $\Delta G = \alpha kT$ per event with $\alpha = \ln(1/\varepsilon_{ev}) \approx 7$ at $10^{-3}$, giving $E_B(D) = \alpha kT D^2$: same $D^2$ scaling, constant $\alpha$ versus 2, so single-stage analogue is cheaper by $\alpha/2 \approx 3.5$ (Sarpeshkar's low-precision analogue win). The mismatch asymmetry is the core: gain mismatch cancels entirely, all-or-none events carrying no amplitude into the count (the dPCR $\beta$ cancellation and p53 $\alpha \to 0$ mechanism); threshold mismatch is logarithmically suppressed, $N = (\tau_r/T)\ln(D_0/D_c)$ giving $\delta N/N = \varepsilon_0/\ln(D_0/D_c)$, a soft ceiling $C_B \le \log_2[\ln(D_0/D_c)/\varepsilon_0]$ ($D_0/D_c = 100$, $\varepsilon_0 = 0.2$: 23 levels $\approx 4.5$ bits versus the 2.32-bit analogue ceiling, stretchable via $D_0/D_c$). Counts regenerate at every cascade stage (von Neumann regeneration), so

$$E_B(L,D) = \alpha kT\,L\,D^2. \tag{S5.4}$$

**Crossover.** With depth $L$ defined as the number of transmission stages through a noisy medium (co-localised stages such as scaffolds do not increment depth), $E_A/E_B = 2L/\alpha$, giving

$$L^{*} = \alpha/2 \approx 3\text{–}4 \text{ stages}, \tag{S5.5}$$

independent of $D$ since $D^2$ cancels; correlated errors move it earlier, $L^{*} \approx \sqrt{\alpha} \approx 2.6$; inside the ceiling region analogue is cheaper still, so (S5.5) is conservative. Shallow systems favour analogue, deep systems favour counting, and the crossover is set by event reliability alone. The biological mapping is consistent: scaffolded PdPC modules at $L \approx 1$ stay graded; freely diffusing kinase cascades (RTK–Ras–MAPK, 3–4 stages) sit at the crossover; myelinated axons at $L \sim 10^2$–$10^3$ implement per-node regeneration; p53 decoding counts. Recorded limitations: the force is in the scaling structure ($L^2$ versus $L$; hard versus soft ceiling), not the coefficients; $\sqrt{L}$ compounding assumes independent errors; Poisson counting assumes independent events (bursts modify the $\sqrt{N}$ exponent); $\alpha$ is an Arrhenius-type barrier model; the comparison covers energy per channel use, not bandwidth, so in high-bandwidth shallow settings the analogue advantage exceeds what is stated (photoreceptors staying graded until the ganglion-cell handoff is consistent). Going digital buys mismatch immunity and depth robustness at a continuing energetic price; the analogue regime remains optimal where transmission is shallow and degenerate scales are pinned by other means.

---


---

## S6. Cross-system consistency and the NF-κB/ERK analysis record

This section collects the model-side and literature-side audits that place the four main-text systems in a common signature frame, and records the two hypothesis-level adjudications (codes 12 and 15) that the framework ran against itself in the NF-κB line. Source archive: `08_audit_trail/01_cell_line_1/code_and_figures/` codes 6/8/9/11/12/15 with their run logs, and `08_audit_trail/01_cell_line_1/code_and_figures/public_data_audit_three_dataset_ruling_report.md` (2026-08-04).

### S6.1 NF-κB oscillator analysis (Krishna et al. model^55^, soft and spiky regimes)

The Krishna et al.^55^ three-variable minimal NF-κB oscillator was audited under the three-signature protocol (onset continuity, hysteresis width, range sensitivity α) in both published parameter regimes.

*Soft regime* (B = 3, δ = 0.005, n = 1; code 6, reproduced as code 12's baseline arm). A single supercritical Hopf onset at C_H = 2.49 × 10⁻³ (code 6 original: 2.62 × 10⁻³); onset amplitude follows A² ∝ (C − C_H) at R² = 0.996; windowed α = +0.729 (code 6: 0.714); hysteresis absent.

*Spiky regime* (B = 954.5, δ = 0.029, n = 1; code 9). The equilibrium branch is unique (six-point root-count check) and bounded by dual Hopf boundaries computed analytically with finite-difference Jacobians: low C_H = 0.00589 (λ = ±4.2951i) and high C_H = 0.16741 (λ = ±17.1915i). A 110-point two-direction continuation gives maximum up/down amplitude difference 2.05 × 10⁻³: hysteresis width zero. Dual-initial-condition probing returns the same attractor in 8/8 cases, correcting the fold/bistability observation in the code 6 header note, recorded as a long-transient artefact (erratum). A fine sweep resolves a canard near-vertical segment at C ≈ 0.00600 (A: 0.061 → 0.175). Windowed α: rising limb +0.549 (window [0.009, 0.032]), canard-inclusive +0.844, falling limb −0.438; period spans 0.37–3.20 (about 8.8-fold). Spikiness: at C = 0.035 the duty above half-maximum is 10.8%, a genuinely spiky waveform.

Three-signature placement for both regimes: continuous onset, zero hysteresis, |α| > 0.3 in every window; both fall in the NF analogue region. The placement is robust to parameter regime, so the gap to the experimental placement of Tay et al.^40^ (α = 0.151, constant first-peak area, dose-encoding pulse count) is a model-experiment discrepancy, not a regime artefact.

### S6.2 The two rejected patches (codes 12 and 15)

Two mechanistic hypotheses for closing the NF-κB model-experiment gap were proposed (by external review) and adjudicated computationally. Both were rejected and are double-recorded.

*Code 12, parameter-unreachability test of the switch-ified IκBα degradation patch.* The hypothesis: replacing the IKK-driven degradation term I/(ε + I) by a switch-like I⁴/(ε⁴ + I⁴) should move the model from NF to EXC (α < 0.3, fold/canard structure). Results: (i) fold/bistable multi-value points: 0/150 grid points in all four parameter combinations; (ii) in the soft regime the oscillation window was abolished entirely (steady state migrates to the saturated degradation branch): oscillation destroyed, not converted; (iii) in the spiky regime the window shifted to (9.11 × 10⁻³, 0.20) but structure was unchanged (relaxation limit cycle, amplitude span 1.95-fold, wide-window α = −0.005, period span 10.8-fold): still an NF-class relaxation oscillator; (iv) excitability probes below the window returned response peaks equal to the perturbation itself (response ratio 3000 = perturbation span 3000): zero amplification, no threshold. Verdict: the patch is falsified; the NF → EXC transition is structural, not a matter of parameter sharpness. A side result: the graded sub-window response also weakens the reading of Tay's damped oscillations as a sub-Hopf NF resting side, which would predict a graded first peak contrary to the stereotyped one observed.

*Code 15, saturated-transient hypothesis.* The hypothesis: transporter saturation (K_T on nuclear import) with saturating stimulus C_max could reconcile the Krishna et al. model with the Tay et al. three numbers^40,55^ (α = 0.151, ν = +0.58 per decade, A2/A1 ∈ [0.5, 0.8]). A deterministic grid of 24 (K_T, C_max) combinations × 5 doses (120 step responses) fits α within ±0.07 in 8 of 24 combinations, but in all of them the count slope stays too high: minimum ν = 0.900 (C_max = 0.008, K_T = 0.03), above the Tay upper bound 0.83. A population-level heterogeneity re-check agrees. Verdict: mechanism acknowledged, quantitatively rejected at both levels.

### S6.3 ERK digitisation and the static-extrapolation falsification

A previous programme version (v0.8) predicted, by mechanistic extrapolation from the static gain grading of the PdPC, that mammalian ERK pulse amplitude should grade with dose in the analogue region (α ≳ 0.5, since the pathway has not crossed g*). Source-figure digitisation of Albeck et al.^37^ (Fig. 3A, MCF-10A, 0/10/50/200 pg ml⁻¹ EGF, single-cell EKAR trajectories; uncertainty ±0.03 EKAR) falsified the prediction: pulse peak heights were constant at −1.03 ± 0.02, −1.04 ± 0.02 and −1.04 ± 0.02 across 10/50/200 pg ml⁻¹ (|α| ≲ 0.05 over >1.3 decades), while pulse count rose at ν ≈ +0.39 per decade over 10–50 pg ml⁻¹ and saturated at 200 pg ml⁻¹; the second cell line (184A1, Fig. S4A) reproduced the pattern. The falsified item is the extrapolation, not a theorem; it is preserved in the version history as the framework self-correction ledger entry (main-text Discussion). The same three-dataset audit confirmed the p53 built-in control quantitatively (α_NCS = 0.00 ± 0.28; α_UV = +0.99 ± 0.32) and partially overturned the NF-κB "frequency does not encode" sentence: peak interval is constant (75–95 min) but oscillation count grows with dose (1.7 → 4.0 over four decades, +0.6 per decade) while first-peak α = 0.151 sits in the leaky-digital band.

### S6.4 Xenopus MAPK and tissue-level p53

Embedding the PdPC loop in positive feedback gives a critical feedback gain g* ≈ 0.443; beyond g* the gain-dose curve folds, the apparent Hill coefficient rises from a baseline of 11 (κ = 0.05, closed form 1 + 1/(2κ)) to about 480 at g = 10, and hysteresis width reaches 0.44–0.58 dex (g = 10–20), while the PdPC degeneracy group remains exact throughout the bistable region. This places the Xenopus oocyte maturation switch in-framework: the MAPK cascade there has n_H ≳ 35 with positive feedback, a near-digital static switch corresponding via n_H = 1 + 1/(2κ) to κ ≲ 0.015. The same ERK molecular identity thus occupies two regimes across organisational levels: a bistable decision element in the oocyte, and fixed-amplitude frequency-modulated pulses in mammalian single cells (S6.3). At tissue level, radiation-resistant gut shows oscillatory p53 (counting mode) whereas radiation-sensitive thymus and spleen show sustained p53 (switch/decision mode); coding architecture matches physiological function.

### S6.5 Msn2-type frequency encoding

The minimal excitable model of Msn2-type frequency encoding (FHN under sustained drive; code 8) shows, inside the activity window [0.32, 1.43] au, amplitude invariance to ±1.0% (α = 0.020) with monotonically increasing frequency, ν_f = 0.128 au⁻¹ (R² = 0.96). Onset is a continuous Hopf followed by a canard explosion (I_expl ≈ 0.334), and a narrow subcritical hysteresis of width 0.0095 (about 0.9% of the window) is detected: the hysteresis signature correctly captures this subcritical exception while the signature triple still places the system on the EXC digital side. This is the model-side counterpart of the P2 adjudication dataset (S7.1).

### S6.6 Morris–Lecar neuronal duality instance table

The encoding-signature duality (signature = f(bifurcation type × observable)) was instantiated at the Morris–Lecar textbook working point. Intrinsic layer G₀: under joint scaling a·(I, ḡ, C) trajectories are pointwise invariant (5/5 waveforms, correlation = 1.0, frequency ratio = 1.0). Transmission layer G₂: under a seven-point synaptic-gain drift sweep of 0.3×–3×, event-time readout degrades by 0.0000% across the whole sweep while amplitude readout changes by a mean of 72.9% (peak 200%): gain scaling dies, event timing survives. Type-I/SNIC: the frequency channel grades by the square-root law while amplitude saturates. Type-II/Hopf (subcritical control): onset is a large limit cycle, neither channel grades, control amplitude slope α = −0.258. The amplitude leg of the duality fails on subcritical Type-II: the duality's domain is bounded by bifurcation criticality and observable, as the discrimination theorem requires.

---

---

## S7. Blind adjudications P2–P7: full record

This section is the verbatim empirical ledger of the seven pre-registered blind adjudications. Every clause had its decision line frozen in writing before unblinding; scoring scripts execute the frozen tables verbatim with fixed seeds (SEED = 20260814 for P2/P3, 20260815 for P4–P7) and bootstrap ×2000 confidence intervals; verdict cards are archived verbatim and never edited, with corrections appended as errata (44b, 47b, 49b). Archive roots: `08_audit_trail/01_cell_line_1/results/P2P3P4_blind_adjudication_local_originals/`, `08_audit_trail/01_cell_line_1/results/Keshelava2018_P5/` and `08_audit_trail/01_cell_line_1/results/Keshelava2018_P5_47berrata_readjudication/`, `08_audit_trail/02_cell_line_2/results/P6_Wang2022_FCD_attribution/`, `08_audit_trail/02_cell_line_2/results/P7_Moore2024_FCD-Weber/` (+ `08_audit_trail/03_cell_line_3/results/P7_Moore2024_FCD-Weber/`).

### S7.0 The ledger

Seven adjudications (P2–P7) on six independent published single-cell datasets carry fourteen pre-registered clauses: P2 is three independent frozen clauses (P2-1 and P2-2 falsified, P2-3 a hit), and P5-3 is scored as a dataset-level limitation rather than a hit/falsification. Totals: 2 clean hits, 5 intermediates, 6 falsifications, 1 dataset-level limitation (14 clauses).

| Adjudication | Clause | Dataset | Frozen decision line | Result value | Verdict |
|---|---|---|---|---|---|
| P2 (code 44/44b) | P2-1: size-tercile stratification penalty of molecule-number channel | Msn2, Hansen & Zechner^17^ (Zenodo 10.5281/zenodo.2755026) | Δ_mol ≥ +0.03 in ≥6/7 promoters | 0/7 promoters pass; 21 units, Δ ∈ [−0.0491, +0.0294], 19/21 bootstrap CIs cross zero | Falsified |
| P2 | P2-2: fluorescence (AU) channel penalty bounded | same | Δ_AU ≤ 0.01 and Δ_AU < Δ_mol in ≥6/7 promoters | 0/7 promoters pass | Falsified |
| P2 | P2-3: event-time channel stratification-insensitive | same | \|Δ_time\| ≤ 0.01 (k = 3) in ≥6/7 promoters | 6/7 promoters pass at k = 3; 7/7 in the k = 2 arm, 6/7 in the k = 4 arm (DCS2 \|Δ\| = 0.0102) | Hit |
| P3 (code 45) | P3-1: dose ordinality by free statistics | NF-κB gradient, Son et al.^18^ (Zenodo 10.5281/zenodo.6858118) | both units τ(3) AUC ≥ 0.60 | 0.6548 (10↔30 ng ml⁻¹, pass) / 0.5184 (30↔100, fail) | Intermediate |
| P3 | P3-2: spatial ordinality (near/mid/far terciles) | same | both units AUC ≥ 0.55 | 0.9070 (mid > near) / 0.8193 (far > mid), both pass | Hit |
| P3 | P3-3: stimulus-duration axis | same | both units AUC ≥ 0.60; all units below line → falsified | 0.5911 (30 > 15 min) / 0.5214 (60 > 30 min), both below | Falsified |
| P4 (code 46) | P4-1: free event-time dose decoding, unsaturated region | ERK-KTR, Chavez-Abiega et al.^19^ (histamine 7,393 cells; UK14304 7,508 cells; 6 doses each) | 4 units τ(3) AUC ≥ 0.60; any reversal (AUC < 0.5) or all below line → falsified | all 8 AUCs below 0.60 (max 0.5819; one at 0.4695 < 0.5) | Falsified |
| P4 | P4-2: saturation-region ordering | same | hit iff AUC_H < 0.60 and D ≥ 0.05; falsified iff AUC_H ≥ 0.60 or D ≤ 0 | AUC_L = 0.5669, AUC_H = 0.5201, D = +0.0468 [+0.032, +0.062] | Intermediate |
| P4 | P4-3: size-gain stratification invariance | same | hit iff ≥8/10 units \|Δ\| ≤ 0.01; falsified iff ≥3/10 units \|Δ\| ≥ 0.03 | 5/10 units \|Δ\| ≥ 0.03 (4/10 units \|Δ\| ≤ 0.01) | Falsified |
| P5 (code 47b, erratum of 47) | P5-1: dead-zone clause (event time τ) | GPCR→Ca²⁺, Keshelava et al.^20^ (Nat. Commun. 9:876); 195 valid cells, 11/27 experiments | hit iff 4/4 τ-AUC ∈ [0.40, 0.60); falsified iff ≥2 units outside [0.40, 0.60] | 2/4 units out of band (0.629, 0.659; in-band 0.568, 0.520) | Falsified |
| P5 | P5-2: live-zone counting clause | same | hit iff 4/4 count-AUC ≥ 0.60; falsified iff 4/4 < 0.60 | 2/4 units at line (0.629, 0.659; below: 0.570, 0.522) | Intermediate |
| P5 | P5-3: paid-peak calibration control | same | ≥3/4 units peak-AUC ≥ 0.60 → calibration pass | only 2/4 units ≥ 0.60 (0.731, 0.744; below: 0.551, 0.523) | Dataset-level limitation |
| P6 (code 48) | P6: fold-change vs absolute-dose attribution | NF-κB sequential stimulation, Wang et al.^21^ | ΔR² ≥ +0.15 hit; ≤ −0.15 falsified; else intermediate (slope-sign constraint) | ΔR² = −0.0014, 95% CI [−0.2518, +0.1533], 18 valid units, power qualifier | Intermediate |
| P7 (code 49, erratum 49b) | P7: fold/Weber vs absolute-total attribution | chemotaxis background ladder, Moore et al.^22^ (Dryad 10.5061/dryad.nvx0k6dzz) | ΔR² ≥ +0.15 hit (fold slope > 0); ≤ −0.15 falsified; else intermediate | first pass: ΔR² = +0.1181 [−0.2956, +0.4677], 25 units; erratum: ΔR² = −0.0909 [−0.4446, +0.2766], 30 units; both CIs cross zero | Intermediate |
| **Totals** | **14 clauses** | **6 datasets** | | | **2 hits / 5 intermediates / 6 falsifications / 1 dataset-level limitation** |

### S7.1 P2: Msn2 amplitude stratification decay (Hansen & Zechner^17^)

Pre-registration `preregistration_P2_amplitude_stratified_decay_prediction_v01.md` (frozen 2026-08-14; seed 20260814; bootstrap 2000; main analysis duration 50 min). The adjudication tested whether cell-size heterogeneity is priced into the molecule-number channel: if invariance must be purchased, stratifying cells into size terciles should degrade dose decoding by molecule-number statistics (P2-1), the fluorescence channel should show a bounded smaller penalty (P2-2), and the calibration-free event-time channel should be stratification-insensitive (P2-3). Units were 7 scored promoters (ALD3, DCS2, DDR2, HXK1, RTN2, SIP18, TKL2) × 3 dose pairs (100–275, 275–690, 690–3000 molecules per cell); two SIP18 mutant promoters were recorded descriptively and never entered scoring.

Results. P2-1: 0/7 promoters reach Δ_mol ≥ +0.03; across the 21 stratification units Δ ∈ [−0.0491, +0.0294], and the appended CI computation (code 44b, which leaves the frozen verdict untouched and whose 63/63-row reproduction check matched the formal run bit-exactly) shows 19/21 bootstrap 95% CIs crossing zero. P2-2: 0/7 promoters pass. P2-3: 6/7 promoters satisfy |Δ_time| ≤ 0.01 at k = 3 (only DCS2 fails, Δ = −0.0102), and the robustness arms pass 7/7 at k = 2 (max |Δ| = 0.0062) and 6/7 at k = 4 (DCS2 at |Δ| = 0.0102 against the frozen 0.01 line); Δ_time is also stable across the four alternative analysis durations (|Δ_mol| ≤ 0.0223 in all duration arms). Verdict: P2-1 falsified, P2-2 falsified, P2-3 hit. The cell-side redemption of invariance pricing does not run through the size-gain axis, and the clause contracted accordingly in the main text.

### S7.2 P3: NF-κB spatial gradient (Son et al.^18^)

Pre-registration `preregistration_P3_free_statistic_sufficiency_NFkB_spatial_gradient_v01.md`; code 45 v0.1.1; seed 20260814; bootstrap 2000. Only per-cell self-referenced free statistics were admitted: first-crossing event time above the cell's own pre-stimulus μ₀ + 3σ₀, and above-threshold duration; any cross-cell calibration was a forbidden operation. Decision lines were derived from the null noise floor before unblinding.

Results. P3-2 (spatial ordinality, near/mid/far terciles): AUC 0.9070 (mid > near) and 0.8193 (far > mid), both above the 0.55 line: clean hit. P3-1 (dose ordinality): 0.6548 for TNF 10↔30 ng ml⁻¹ (pass at 0.60) but 0.5184 for 30↔100 ng ml⁻¹ (fail): intermediate, because the pre-registration scores falsification only when all units fail. P3-3 (stimulus-duration axis): 0.5911 (30 > 15 min) and 0.5214 (60 > 30 min), both below the 0.60 line: falsified by the literal frozen table. The paid peak-amplitude descriptive arm still decodes where the free statistics fail (peak AUC 0.7610 at 30↔100 ng ml⁻¹; 0.7130 at 10↔30), and the count-channel descriptive arm sits at 0.5791/0.5094. Response rates by dose were 0.831/0.950/0.991 (10/30/100 ng ml⁻¹). Overall: P3 falsified (a falsified clause exists), with the spatial hit and the dose-saturation boundary recorded verbatim.

### S7.3 P4: ERK re-adjudication on GPCR-KTR dose series (Chavez-Abiega et al.^19^)

Pre-registration `preregistration_P4_ERK_readjudication_GPCR_KTR_trajectory_level_v01.md` (v0.3, re-anchored to 6 dose levels after schema verification in zero-statistics windows); code 46 v0.1.5; seed 20260815; bootstrap 2000; run 2026-08-14T22:28. Main-analysis ligands: histamine (7,393 cells included of 7,393; zero exclusions under the frozen rules) and UK14304 (7,508/7,508), six doses each. The adjudication asked whether the digital-side placement of ERK (S6.3) survives at trajectory level: free event-time statistics τ(3) should decode dose ordinality in the unsaturated region (P4-1), the saturated region should show ordering loss (P4-2, the cross-system test of the P3 boundary statement), and the P2-contracted stratification clause should hold (P4-3).

Results. P4-1: all eight τ(3) AUCs fall below the 0.60 line (histamine 0.5819, 0.5633, 0.5636, 0.4695; UK14304 0.5561, 0.5663, 0.5377, 0.5096), including one reversal below 0.5: falsified. P4-2: AUC_L = 0.5669, AUC_H = 0.5201, gradient D = +0.0468 with 95% CI [+0.032, +0.062]; D is positive but below the 0.05 hit line: intermediate. P4-3: 5/10 units have |Δ| ≥ 0.03 (histamine −0.1091, +0.1311, −0.1190; UK14304 −0.0393, +0.0536), exceeding the ≥3/10 falsification line; 4/10 units have |Δ| ≤ 0.01: falsified, so the size-gain axis is priced in mammalian ERK though not in yeast Msn2. Descriptive arms (not scored): S1P series (7,891 cells), Akt channel, paid peak/integral amplitudes, inhibitor arms. Overall: P4 falsified. The failure is consistent with the morphology risk registered before unblinding (near-equal times-to-peak in this transient single-peak morphology).

### S7.4 P5: boundary predictor on GPCR→Ca²⁺ transients (Keshelava et al.^20^)

Pre-registration `preregistration_P5_boundary_predictor_first_test_GPCR_calcium_transient_v01.md` (v0.4; v0.1 frozen 2026-08-14, with two zero-statistics-window methodology revisions v0.2/v0.3 recorded, and the v0.4 erratum appendix; seed 20260815; bootstrap 2000). The boundary predictor v1 was tested on HEK293 M3R Fura-2 within-cell dose ladders (7 acetylcholine doses, 5 repeats each, 27 experiments): P5-1 dead-zone clause (event time τ should not decode: AUC ∈ [0.40, 0.60)), P5-2 live-zone clause (response count should decode: AUC ≥ 0.60), P5-3 paid-peak calibration control.

Inclusion record. The frozen timing-quality rules excluded 16 experiments (no valid anchor/chain); 11/27 experiments entered. The original code 47 run (208 cells, 0 further exclusions) returned an as-run verdict of intermediate. A post-adjudication audit against the authors' independent SD3 peak-value truth table then identified a feature-implementation defect: the baseline-referenced peak definition counted the recorded whole-run background drift (+0.25 ratio units) and flow-switching artefacts as response (the D1 response rate, ~6% by the authors' count, was misreported as 87.5%). By the erratum rule, the v0.2.0 as-run verdict was voided-by-erratum, double-recorded and never deleted, and the re-adjudication code 47b under the corrected inclusion rule is the authoritative verdict: 195 valid cells (13 excluded), with the power-reduction qualifier 195 ∈ [150, 200) recorded on the card.

Authoritative (47b) results. P5-1: τ-AUCs 0.629 and 0.659 (D1↔D2, D2↔D3) fall outside the dead-zone band, 0.568 and 0.520 (D5↔D6, D6↔D7) inside: 2/4 units out of band, falsified. P5-2: count-AUCs 0.629, 0.659 at line and 0.570, 0.522 below: 2/4, intermediate. P5-3: peak-AUCs 0.731 and 0.744 at line, 0.551 and 0.523 below: only 2/4, so the dataset-level weak-decoding limitation is attached. Overall: P5 falsified (the boundary predictor dies on this morphology; the dying side is the dead-zone clause). The population descriptive arm, recorded alongside, is unambiguous: response rate rises 0.015 → 0.272 → 0.544 → 0.769 → 0.826 → 0.892 → 0.918 and mean count 0.02 → 0.87 → 2.24 → 2.70 → 3.62 → 4.13 → 4.36 across D1–D7 (100 nM to 10 µM). Population-level counting is intact where per-cell decoding is weak; the clause contracted accordingly.

### S7.5 P6: fold-change attribution on NF-κB sequential stimulation (Wang et al.^21^)

Pre-registration `preregistration_P6_FCD_attribution_adjudication_NFkB_sequential_grid_v01.md` (frozen 2026-08-15, zero decision statistics computed at freezing); code 48 v0.1.0; seed 20260815. Units are the S2 attenuation A_c = med(peak_S2 | condition)/med(peak_S2 naive) for ordered ligand pairs × dose groups (TNF/IL-1β/LPS/PAM; 4 pulses at 120 min intervals); the pre-registration designed 36 units and 18 valid units entered scoring, with the power-limitation qualifier recorded on the card. Competing univariate models A ~ log₁₀(fold) versus A ~ log₁₀(abs1); frozen line ±0.15 on ΔR², with a fold-slope sign constraint.

Result. ΔR² = R²_fold − R²_abs1 = −0.0014 (R²_fold = 0.0063, R²_abs1 = 0.0078; fold slope −0.0086), bootstrap 95% CI [−0.2518, +0.1533] crossing zero: intermediate, power-limited. The self-built feature arm (local-baseline-referenced peaks, pre-registered robustness arm ①) gives Δ = −0.1783 (R²_fold = 0.0608, R²_abs1 = 0.2391, n = 18), which would cross the falsification line; per the pre-registered divergence rule, the author-feature main judgement stands, and the self-built arm is double-recorded here as a non-scoring control and presented alongside in the Discussion. The fold-change-versus-absolute attribution on this dataset is undecidable at the available power and is reported rather than resolved post hoc.

### S7.6 P7: fold-change/Weber attribution on the chemotaxis background ladder (Moore et al.^22^)

Pre-registration `preregistration_P7_FCD_Weber_verdict_chemotaxis_single_cell_FRET_v01.md` (frozen 2026-08-15; schema anomalies registered before unblinding, including the 230831_FOV2 background conflict resolved by measured s values); codes 49 v0.1.0 and 49b v0.2.0 (erratum); seed 20260815. Units are (B, F) condition pairs with B > 0 on the MeAsp background ladder (6 backgrounds × 5 foregrounds = 30 designed units; B = 0 group descriptive only). Response R(B,F) is the per-condition median of per-cell Δa (FRET activity drop on attractant step, 7 repeats per level). Competing models R ~ log₁₀ r (fold r = (B+F)/B) versus R ~ log₁₀ T (total T = B+F); frozen line ±0.15 on ΔR².

First pass (code 49). 25 condition units (the B = 0.01 µM group measured at B = 0 in the affected files entered at their measured background); ΔR² = +0.1181 (R²_fold = 0.2145, R²_abs = 0.0964, fold slope +0.1399), 95% CI [−0.2956, +0.4677]: intermediate. Fold-total separability on this ladder is limited by design, sep_rho = −0.4553.

Erratum re-adjudication (code 49b). A unit-definition/background-assignment correction (the B = 0.01 µM group correctly assigned, per the pre-registered measured-s rule) brought the full 30 units into scoring: ΔR² = −0.0909 (R²_fold = 0.0815, R²_abs = 0.1724, fold slope +0.0549), 95% CI [−0.4446, +0.2766], sep_rho = −0.4450: intermediate again. Both CIs cross zero and the point estimate changes sign between the passes; robustness arms (amplitude-feature ΔR² +0.0174/−0.0768; cell-level ΔR² 0.0) do not rescue a decision. The amplitude table does not adjudicate its own coordinate system, which is precisely the ambiguity the receptor-layer rupture (main-text Results 6) makes structural. The erratum pair 49/49b is archived in full, with the 49b unit table also deposited in the 03 lineage folder.

### S7.7 Cloud re-execution records

Key adjudications were re-executed independently in a clean cloud environment on 2026-08-15 from the frozen scripts, against the same archived datasets: P2 and P3 (`08_audit_trail/01_cell_line_1/results/P2P3_cloud_recheck_2026-08-15/`, codes 44 and 45 with unit tables, descriptive/robustness arms, judgement logs and verdict JSONs) and P4 (`08_audit_trail/01_cell_line_1/results/P4_cloud_recheck_2026-08-15/`, code 46). The cloud verdict cards reproduce the local originals: identical frozen margins and seeds, identical point estimates (e.g. P2 promoter-level Δ values bit-identical; P4-2 AUC_L = 0.5668941796791128 identical to the local card), and identical verdicts (P2-1/P2-2 falsified with P2-3 hit; P3 overall falsified; P4-1/P4-3 falsified, P4-2 intermediate). No cloud-local discrepancy exists in the ledger.

### S7.8 Pre-registration texts and the ERK extrapolation falsification

The pre-registration texts are archived at `08_audit_trail/05_programme_design/preregistration_P2_amplitude_stratified_decay_prediction_v01.md`, `preregistration_P3_free_statistic_sufficiency_NFkB_spatial_gradient_v01.md`, `preregistration_P4_ERK_readjudication_GPCR_KTR_trajectory_level_v01.md`, `preregistration_P5_boundary_predictor_first_test_GPCR_calcium_transient_v01.md`, `preregistration_P6_FCD_attribution_adjudication_NFkB_sequential_grid_v01.md`, and `preregistration_P7_FCD_Weber_verdict_chemotaxis_single_cell_FRET_v01.md`. Each contains a blindness declaration (inventory of schema-level inspection before unblinding), the frozen decision table, guardrail/empty-run rules, and revision records confined to zero-statistics windows. Three adjudications carry errata appended after freezing under the never-move rule: 44b (CI recomputation, no verdict change), 47b (feature-implementation defect, re-adjudication authoritative), 49b (unit-definition correction, both passes intermediate). One pipeline error (P5) and one guardrail-empty-run history (P4 v0.1.2, abolished before running) are double-recorded. Separately, the ERK static-extrapolation falsification (S6.3) is preserved in the version history as the programme's one failure of mechanistic extrapolation across readout levels: the v0.8 prediction (amplitude α ≳ 0.5) was overturned by source-figure discrimination on Albeck et al.^37^ (α ≈ 0), and the record is kept because the falsifiability of the framework is a fact of the archive, not a rhetorical posture.

---

## S8. Chemotaxis receptor-layer rupture (main-text Results 6)

### S8.1 Data and statistics construction

The dataset is Moore et al.^22^ (Dryad doi:10.5061/dryad.nvx0k6dzz, CC0): single-cell FRET dose-responses of the bacterial chemotaxis receptor layer across seven background levels B = 0 / 0.01 / 0.1 / 0.3 / 1 / 10 / 100 µM MeAsp, five foreground step levels per background, 57 to 268 cells per background level. Two statistics were constructed from the same raw trajectories. (i) The amplitude table R(B, F): per-condition population medians of the plateau response, assembled as a 30-row unit table. (ii) The per-cell midpoint distribution K1/2(B): for each cell, the foreground concentration at half-maximal response, interpolated on the log-concentration axis; cells whose maximal observed response da_max < 0.5 were right-censored. Per-background K1/2 medians were 2.03 / 2.90 / 2.17 / 2.58 / 3.49 / 13.92 / 119.56 µM for the seven backgrounds. Under the standard model class (MWC activity core with methylation-based perfect adaptation), both statistics must be generated by a single parameter set; the audit tests whether any such set exists.

### S8.2 Pareto-front protocol and acceptance region

Joint fits of the amplitude table and the K1/2 distribution were scanned along a Pareto front by varying the weight λ between the amplitude R² objective and the K1/2 absolute-error objective (in dex). The acceptance region, fixed before scoring, requires amplitude R² > 0.9, K1/2 error < 0.1 dex, and a physical saturation amplitude amax approximately equal to 1. The endpoints of the front are:

| λ | amplitude R² | K1/2 error | parameters |
|---|---|---|---|
| 0 | 0.980 | 0.622 dex | Ki = 30.7, Ka = 133, N = 39, amax = 1.39 |
| 30 | 0.511 | 0.074 dex | Ki = 4.3, N = 2.0 |

No point on the front enters the acceptance region. Fitted amax values of 1.4 to 2.7, where the physical value is approximately 1, are a symptom of model strain. The rupture is thereby quantified: no standard MWC-plus-perfect-adaptation parameter set is jointly consistent with the amplitude table and the midpoint distribution.

### S8.3 Eight repair paths and the double-recorded pipeline errors

Each repair path below is a minimal extension or perturbation of the standard model, scored under the uniform acceptance criterion of S8.2 (amplitude R² > 0.9 and K1/2 error < 0.1 dex and amax approximately 1).

| # | Hypothesis | Test | Result | Why rejected |
|---|---|---|---|---|
| 1 | Two subpopulations (fast-adapting plus non-adapting; models M2a/M2b) | Joint fit with mixture weight α free | Directional failure; α driven to 0.05–0.08; worse than the single-species fit | The non-adapting subpopulation's zero-ligand-anchored increment shape conflicts with the amplitude table; the failure is directional, not a matter of tuning |
| 2 | Imperfect adaptation (M2c, β < 1) | β added to the adaptation core | Best compromise R² = 0.953 / K1/2 error 0.182 dex at λ = 1 | The improvement derives mainly from Ka → ∞ rather than from β; remains outside the acceptance region |
| 3 | Protocol mismatch (step vs waveform pre-adaptation) | Re-analysis by protocol shape | Step-protocol responses are logF-dominated; waveform FCD plateau appears only at B in [10, 100] µM | Excluded as the primary cause of the rupture |
| 4 | FRET readout-layer nonlinearity (variants V1/V2/V3) | Power-law readout remap | V3 selects a stable q ≈ 2.0 in the expansive direction; best 0.950 / 0.257 | Readout nonlinearity excluded as the sole cause |
| 5 | TCS ligand depletion | Depletion ratio κ ≈ 48 ≫ 1 examined | After solving the implicit free-ligand equation, the Pareto front is digit-for-digit unchanged | Formally excluded |
| 6 | Background-dependent gain retuning (G-Ki / G-N / G-both) | Effective parameters free per background | K1/2 side reaches 0.052 dex but amplitude R² collapses to 0.68; N0 and h hit non-physical boundaries | The model meets the midpoint statistic only by abandoning the amplitude table |
| 7 | Heterogeneous population plus censoring statistics (code 60) | Population spread σa free | Best 0.946 / 0.179; the optimizer actively shrinks σa* to 0.11 | The fitted direction is opposite to the hypothesis |
| 8 | Asymmetric cooperativity (T-state N ≈ 2–4, R-state N ≈ 12–39) | Per-row apparent Hill coefficient n(B) vs single-N MWC controls in the same window | Measured n(B) rises from ~1 to 17, same form and magnitude as single-N apparent n; everywhere within or below the envelope of the N = 6 and N = 39 models | No steepness beyond single-N MWC exists in the data; two N values are not required |

Three pipeline errors of our own were caught and repaired mid-course, all double-recorded: code 55 misaligned the β/amax unpacking order (exposed by an independent rerun); code 59 fed the total concentration sv to the model as the foreground F in the pooled pipeline (exposed by a single-point spot check; repaired as F = sv − B); code 57 clamped parameter bounds so that all multi-starts were infeasible (repaired). All three were repaired and independently rerun; no issued verdict moved.

### S8.4 The Weber line

The measured per-cell midpoints obey a textbook Weber line, K1/2 = 1.17·(1.95 + B) µM, across the seven backgrounds (maximum deviation 0.10 dex, mean deviation 0.03 dex; the B = 10 point lies on the line at −0.002 dex). Inverting the Weber slope c = e^{g*} − 1 against the model class requires cooperativity N ≈ 2–4, while the amplitude table requires N ≈ 12–39. The clean statement of the rupture is that the midpoint statistic sees a shallower dose-response than the amplitude statistic. A tail-statistics caveat is recorded with citation discipline: at B = 0 and B = 100 the largest observed response is only 0.44–0.48, so estimable cells cross half-amplitude with the help of the noise or upper tail, and 45–68% of cells per background are right-censored; the K1/2 median is therefore a sensitive-subpopulation statistic while the amplitude table is a whole-population median. This description stands as fact, but its standing as an explanation of the rupture was withdrawn by S8.5 (double-recorded retraction of explanatory status).

### S8.5 Same-source verdict (code 62)

Moore et al.^22^ is single-cell data, so the amplitude table can be rebuilt from strictly the same non-censored cells (32–55% per background) that yield estimable K1/2 values. Re-running the joint fit on this same-source construction leaves the Pareto front unmoved:

| λ | same-source (sensitive subpopulation) | whole-population control |
|---|---|---|
| 0 | R² = 0.971 / 0.573 dex | R² = 0.976 / 0.624 dex |
| 1 | R² = 0.868 / 0.205 dex | R² = 0.848 / 0.261 dex |
| 30 | R² = 0.599 / 0.071 dex | R² = 0.525 / 0.074 dex |

The rupture persists inside a strictly single population. It is not a population-averaging or statistic-sampling artefact; it points to a genuine structural omission in the standard MWC-plus-perfect-adaptation framework itself. This design preemptively closes the most common referee objection, population measurement bias.

### S8.6 Prospective mixed-ligand arm (code 70; pre-registration P1/P2)

Predictions for two unseen mixed-ligand arms were frozen on 2026-08-16, before downloading the seven .mat files of arms 230822/230823/231017/231018. The frozen inputs were the G10 empirical kernel R(F, B) = 0.2268·ln(1 + F/0.570) / (1 + (B/69.7)^0.778), the cross-ligand potency ratio ρ = 23 (1 µM L-Asp equivalent to 23 µM MeAsp; effective background B_eff = B_MeAsp + 23·B_L-Asp), and the Weber line of S8.4. Arm V1 (MeAsp foreground on 10 µM L-Asp + 100 µM MeAsp background, B_eff = 330 µM, 120 cells) and arm V2 (L-Asp foreground on 100 µM MeAsp + 0.1 µM L-Asp background, B_eff = 102.3 µM, 156 cells) were scored per clause. File identity was verified by the s/s_min_max fields against the condition list before any scoring; one extraction-convention error of ours (V1's s field records total MeAsp, since foreground and background are the same ligand) was caught by this check and repaired, double-recorded.

Midpoint clauses (P-2, total effective axis, frozen tolerance ±0.15 dex): all three scored readings hit. V1: prediction 388 µM, measured median 401.3 µM (IQR 388.5–407.6), deviation +0.014 dex, with the predicted high censoring observed (21/120 estimable, 82% censored). V2: prediction 122 µM, measured 133.2 µM (IQR 122.8–140.4), +0.038 dex (71/156 estimable, 54% censored). V2_alt, the ambiguity control in which s is read as total L-Asp: measured 130.8 µM, +0.031 dex. The foreground-axis point prediction was 0.86 µM L-Asp against measured F* of 1.342 µM (V2) and 1.241 µM (V2_alt), that is, 1.24–1.34 µM. The hit points B_eff = 330 and 102.3 lie outside the original fitting domain (B ≤ 100), so these are true extrapolations.

Amplitude clauses (P-1b, frozen tolerance 0.08): both primary arms missed. V1 maximum deviation 0.105 (measured below prediction at F = 5–20, converging at F = 40–80); V2 maximum deviation 0.126 (F = 0.2 measured 0.087 vs predicted 0.213; all four levels with F ≥ 0.6 within tolerance). The miss reproduces the known low-F shortfall of the G10 kernel, previously recorded in the same-source arms; the high-F ends of both arms converge. The ambiguity control V2_alt (maximum deviation 0.069) would pass the amplitude tolerance; it was recorded as a non-scoring control, since the G10 reconciliation alone cannot distinguish the two readings of s.

Structural clause (P-1a, decision line R²(logF) − R²(logr) ≥ +0.3): missed in all readings. V1: R²(logF) = 0.954 vs R²(logr) = 0.959, difference −0.005 (a tie; in the B_eff = 330, F ≤ 80 window the ratio coordinate is compressed into a 0.006–0.086 range and the two coordinates are mathematically near-degenerate, so the +0.3 line is physically unreachable in this regime; this pre-registration design limitation is double-recorded). V2: +0.177; V2_alt: +0.202, right direction, below the line.

Net account: P-2 hit 2/2 arms with combined deviation at most 0.04 dex; P-1b missed in the kernel's known weak region; P-1a missed because the test loses power in this regime. The fracture-side observation item (non-scored): the amplitude-versus-K1/2 coordinate structure persists under mixed backgrounds, with no evidence of a new artefact.

### S8.7 Registered open items

Three items are registered here; item (i) was adjudicated on archived evidence on 2026-09-27 as detailed below, while items (ii) and (iii) remain open. (i) Mechanism shortlist for the rupture's steepness source, narrowed to two candidates: active, methylation-state-dependent retuning of receptor-cluster effective parameters (a class-III paid-position candidate), versus passive consequences of MWC cluster-state equilibrium (class II); adjudicated on archived evidence on 2026-09-27 (verdict card `verdict_chemotaxis_rupture_two_candidates_2026-09-27.md`): the minimal form of candidate (a), active methylation-dependent retuning, is excluded twice over, by repair path 6 (best variant reaches 0.052 dex on the midpoint side but collapses the amplitude side to 0.68 with N0 and h driven to unphysical bounds) and by the cross-ligand arm, where a constant ρ ≈ 22–24 pure-competition MWC suffices (residual −0.045 at F = 5), leaving no room for ligand-specific retuning; candidate (b), passive MWC cluster-state equilibrium, is qualitatively compatible (the measured apparent n(B) curve lies within the single-N MWC envelope throughout) but quantitatively non-closing, so the rupture stands; the decisive experiment, a cheR/cheB fixed-methylation mutant pair measured on both statistics, is registered as follow-up item C3 and is outside the archived WT data. (ii) Citation debt: three references (Bellman–Åström; Hunt–Stein/Lehmann; Lapidoth–Shamai) must be repaid before line-1 unsealing. (iii) The P8 prospective cohort registration (sealed 2026-08-18): frozen parameters c = 0.279, Ki_ser = 0.182 µM, ρ = 23, and the G10 kernel; clauses P8-Q1 (serine shift curve T(B) = 1 + 0.279·(0.182 + B), tolerance ±0.1 dex; reference points B = 0.3 → 1.135, B = 3 → 1.889, B = 10 → 3.840; the control prediction T ≡ 1.051 separates from the shift prediction by more than 0.25 dex at B ≥ 3, so a Dryad serine arm with B ≥ 3 is the priority target and is awaited), P8-Q2 (cross-family zero shift, tolerance ±0.05 dex), P8-Q3 (literature affinity-ratio check, ρ in 15–35), and P8-Q4 (fracture-recurrence observation item, non-scored). Three already-judged items using seen data (competition-arm zero shift 0.051 = 0.051 µM, Δ = +0.002 dex, 273 cells; ser B = 1 shift T = 1.172 vs zero-fit prediction 1.330, Δ = +0.055 dex; L-Asp effective background B_eff ≈ 235 µM equivalent at ρ ≈ 23) are recorded as posterior checks and may not be cited as prospective evidence.

Source archive: `08_audit_trail/03_cell_line_3/results/report_cell_line3_scale_degeneracy_chemotaxis_rupture_detailed.md`; codes 51 (chemotaxis), 53–62, 70; `preregistration_mixed_ligand_validation_arm_P1P2.md`; `preregistration_P8_prospective_cohort_registry.md`.

---

---

## S9. GPCR meta-analysis (main-text Results 7)

### S9.1 Admission-ticket rule and the PTH1R rejection

A system is auditable if and only if three conditions hold: (i) redundant observables sharing upstream parameters exist; (ii) an absolute calibration anchor exists (an independent affinity measurement; the Moore et al.^22^ amax lesson is that without an anchor, marginal statistics are unfalsifiable); (iii) raw or single-cell-level data are public. The GPCR biased-signalling field passes: same-ligand multi-pathway concentration-response data supply redundancy, competitive-binding pKi values supply the anchor, and SI-level parameter tables supply data. The first candidate, PTH1R (Sachdev et al.^56^), was rejected for lack of an independent affinity anchor and the rejection double-recorded: its τ/KA fits hit the parameter boundary at +5.77/+9.58 dex, a numerical artefact recorded verbatim as the teaching case that anchorless systems are unfalsifiable. Three subjects were then admitted: D2R (Klein Herenbrink et al.^43^), AT1R (Wingler et al.^44^), and μOR (Gillis et al.^45^), plus one meta-level subject, the Biased Signaling Atlas.

### S9.2 D2R: identities, error propagation, and the three-ligand fracture (codes 63/64/68/69)

Code 63 (v1.1.0) reconciled binding pKi (spiperone tracer) against functionally inverted pKA with full error propagation across 7 ligands × 6 pathways. The fracture is localized to the three high-affinity partial agonists: bifeprunox Δ = −3.5 dex at 21–23σ, aripiprazole −2.6/−2.9 dex, cariprazine −1.2/−1.4 dex; fast ligands and full agonists are clean. Functional pKA is conserved across arms (0.08–0.31 dex) and serves as the positive control. The bias factor of bifeprunox relative to ropinirole reverses from −0.42 dex at 2 min to +1.93 dex at 90 min, a 2.35 dex swing with a sign change. The fracture position coincides with the load-bearing position of all "significant bias" reports in the study.

Code 64 is a negative record: an attempt to recover upstream rate ordering from single-concentration raw time series failed quality control (single concentration, low replication, baseline drift cannot support waveform normalisation); it is recorded as null and used as evidence in no direction.

Code 68 tested the kinetic explanation with zero free parameters: non-equilibrium occupancy ρ(L, t) built from the authors' own measured kon/koff predicts the pEC50 trajectories of 7 ligands × 3 assays × 8 time points, with the t → ∞ floor at pKd + log(1 + τ). The drift component is quantitatively explained for bifeprunox (model +1.4–1.5 vs measured +1.1–1.9 dex), giving the slow-ligand half of the bias reversal a mechanism. The fast ligands' negative drift (−0.9–−1.4 dex) is not explained by binding kinetics. The absolute-discrepancy component survives: in the cAMP/Gαo arms the model exceeds the measured values by 2.0–2.5 dex at every time point, robust across both tracers. Cariprazine was downgraded to undecidable, because the paper's two tracers differ 82-fold in kon and the existence of the fracture depends on which tracer is believed; this overturns part of code 63's confidence and is double-recorded.

Code 69 mapped pKA(arm, time) cell by cell: cAMP and Gαo are mutually consistent throughout; whole-cell impedance (CI) is systematically high by 1.4–2.2 dex for all 7 ligands; late-time CI nearly closes for aripiprazole and cariprazine, but bifeprunox remains 0.91 dex off under its most favourable combination, and S-3PPP overshoots in the reverse direction by 1.49 dex. The purified fracture statement: for the same ligand in the same cells, cAMP reads bifeprunox as a low-affinity strong partial agonist (pKA 7.8, Emax 88%) while CI reads it as a high-affinity weak partial agonist (pKA 9.45, Emax 44%); no single (KA, τ) reconciles the two. The final decomposition has four components: (i) slow-ligand binding kinetics (resolved); (ii) fast-ligand signal sag (unresolved); (iii) assay-system cross-arm displacement, CI high by 1.4–2.2 dex across all ligands (unresolved); (iv) the bifeprunox residual binding-function gap of 0.9–2.6 dex (unresolved, dual-tracer robust). The fracture has thereby been converted from a number into a structural map.

### S9.3 AT1R: cross-arm reconciliation (code 65)

For Wingler et al.^44^ no independent affinity anchor is usable as a decisive test, so the audit switched to direct cross-arm reconciliation: the same ligand's functional pKA must be conserved between the Gq arm and the β-arrestin arm. Verdict: TRV026 has pKA_Gq = 6.11 vs pKA_arr = 8.16, Δ = −2.05 dex at 5.0σ in wild type, replicated in the L112A mutant (−1.53 dex, 4.9σ); TRV055 is bounded at ≥ 1.6 dex by the full-agonist-arm inequality; the binding pKi = 7.50 sits exactly between the two arm readings. The D2R fracture is binding-state versus functional-state; the AT1R fracture is arm versus arm: two topologies, one criterion.

### S9.4 μOR: the GIRK-column fracture (code 66)

Gillis et al.^45^ main tables 1–3 plus SI table S1 report Emax, pEC50, log τ, and log(τ/KA) for the same ligand × assay units; the Black–Leff^27^ operational model supplies two independent fit-free identities within each cell. Verdict: the GIRK column is systematically fractured across all ligands. Table-3 log τ values sit 0.68–0.97 dex below the Emax-inverted values (3.9–8.4σ); identity 2 locates the GIRK column uniquely (4/5 ligands at 2.4–8.2σ); in the cross-assay pKA map GIRK is the outlier arm for every ligand (up to 12.4σ). An incidental observation, verified against the original PDF: the GPA and cAMP columns of SI table S1 are numerically identical down to the SEM. The load-bearing coincidence recurs: the paper's core conclusion (the GIRK-bias therapeutic window of SR-17018 and related ligands) rests on exactly the one fractured column. The relation to the Stahl & Bohn^46^ reanalysis is recorded but not adjudicated: same point, same direction; the audit's deliverable is the quantified inconsistency map, not a ruling on which side is correct.

### S9.5 Atlas meta layer (codes 67 and 71)

The Biased Signaling Atlas (9,041 rows, 214 papers) was audited at two levels. Its own arithmetic pipeline is exact: the Δlog(τ/KA) triangle closure holds at 3,931/3,931 with strictly zero violations (code 67). The archived literature parameters, however, collectively fail their own model identities (code 71, full-library identity pipeline): 3,096 entry rows grouped into 269 assay groups with 17,987 ligand pairs; 21.5% of pairs violate |Δlog(τ/KA) − Δlog(Emax/EC50)| by more than 0.3 dex and 4.2% by more than 1 dex; group medians exceed 0.3 dex in 70/269 groups; paper medians exceed 0.3 dex in 15 of 66 multi-ligand papers (23%). An Emax quality-control census flags 711 entries above 110%, 181 above 150%, 30 above 300%, and 8 negative. The meta-level finding for the Gillis case: the Atlas ingested only the β-CNA-processed GIRK column, so the fractured unprocessed column is silently absent from the database; the database's inclusion choices change the evidence base of the paper's conclusion. Practical consequence: any bias factor extracted from the Atlas should first pass cell-level identity reconciliation.

### S9.6 Robustness triple (code 73)

Three independent robustness layers. Module A, external blind negative control (haemoglobin oxygen binding, a healthy system): Hill inversion gives p50 = 26.69 mmHg and nH = 2.67; identity deviations are 0.0024 dex median and 0.0026 dex maximum with 0/5 alarms above 0.3 dex; the five-source p50 spread is 0.0065 dex; an injected positive control (half-points scaled by 2.5) gives maximum |Δ| = 0.615 dex with 3/5 alarms. The audit is silent on a healthy system and has teeth on a damaged one. Module B, group-level aggregation with common-mode correlation calibration: among 241 groups with n ≥ 3, high-confidence fractures (p < 0.001 at both ρ = 0 and ρ = 0.5) number 134/241 (55.6%, under the wider criterion that counts a group containing at least one unrepairable entry); under the most conservative spread assumption (SIG = 0.10 dex, ρ = 0.5) the count is 68/241 (28.2%), so the headline conclusion survives the legal-dispersion assumption. Module C, column-level jackknife: 65.7% of flagged groups keep their verdict under leave-one-out; the jackknife score is 1.0 in 195 groups, between 0.5 and 1 in 46 groups, and below 0.5 in none; for the Gillis GIRK column (n = 5) the jackknife score is 1.0, with 4/4 cells still above 0.3 dex under leave-one-out, so the fracture is not driven by a single ligand. The per-group confidence matrix is archived as atlas_confidence_matrix.csv.

### S9.7 Method false-positive self-check (code 72)

Synthetic null batteries were run before any verdict was issued. M1, Gillis-structured nulls at the paper's own SEM magnitudes (pEC50 ± 0.2, log τ ± 0.2, Emax ± 5, log(τ/KA) ± 0.15): identity-1 false-positive rate 0.0543, identity-2 rate 0.0091, against the observed GIRK column rates of 5/5 and 4/5 (whole-table |Δ| > 0.3 fractions 39%/38%). M2, legitimate Hill-coefficient variation (nH in [0.8, 1.25]): median spurious deviation 0.004 dex, with 0.000% above 0.3 dex; nH dispersion cannot produce the Atlas's 21.5%. M3, the zero-model criterion (pKi − pEC50 > 0.3 dex): false-positive rate 0.00003 over N = 200,000 draws, against the D2R observation of 3/7 ligands. M4, detection power (identity-2 structure with injected systematic fracture δ): per-cell detection rates 0.04 at δ = 0.2, 0.07 at 0.3, 0.24 at 0.5, 0.65 at 0.8, 0.87 at 1.0, and 1.00 at δ ≥ 1.5 dex; the criterion is conservative, missing most fractures below 0.5 dex and detecting reliably only at ≥ 0.8 dex or through joint same-direction column evidence. M5, threshold robustness on the measured Atlas distribution: group/paper exceedance rates are 39.03%/37.88% at a 0.2 dex threshold, 26.02%/22.73% at 0.3, 14.87%/13.64% at 0.5, 7.81%/6.06% at 0.8, and 5.20%/3.03% at 1.0; the qualitative conclusion that roughly one fifth to one quarter of groups and papers are inconsistent is stable across the 0.2–0.5 dex band. M6, joint column evidence and direction symmetry: the null probability of the GIRK column's 5/5 identity-1 alarms is 4.71 × 10⁻⁷, and of the 4/5 identity-2 alarms approximately 3.40 × 10⁻⁸ (binomial approximation); alarm directions split 857/899, consistent with the expected 50/50 under the null.

### S9.8 Repair layer (codes 77 and 80)

Code 77 (constraint projection plus fracture probability) maps every fractured record to its nearest identity-consistent point. Across the 269 groups the repair-cost RMS median is 0.072 dex, and 24.2% of groups require a repair larger than 0.3 dex. A two-component mixture fit by EM (converged at 174 steps) gives π_healthy = 0.745, σ_h = 0.0885 dex, σ_b = 0.683 dex; 15.0% of entries carry fracture probability above 0.9, against 14.5% with |deviation| > 0.3 dex. Repair cost and fracture probability agree at Spearman ρ = 0.962. Per-entry outputs are archived in atlas_repair_projection.csv.

Code 80 gives the μOR fracture a minimal, coordinate-aware repair. The extension M1 admits one shared reporting/function coordinate δ_arm per assay arm on top of the latent coordinate, δ_arm ~ N(0, σ_δ²), and is compared against the classic model M0 by marginal likelihood. On the primary dataset (66 residuals, duplicate cAMP logR excluded), empirical-Bayes σ_δ = 0.374 dex and ΔlogML = +128.2 (M0 logML = −121.638; M1 logML = +6.589). The posterior GIRK coordinate is δ_GIRK = −0.805 ± 0.049 dex (95% interval [−0.902, −0.708]), stable under leave-one-ligand-out within [−0.843, −0.732]. Residuals must be quoted with their scopes: for the GIRK column specifically, the median absolute residual falls from 0.811 to 0.107 dex (maximum 0.249 dex); for the whole table, the median absolute residual falls from 0.140 to 0.065 dex. A secondary cAMP coordinate, δ = +0.508 ± 0.163 dex, is recorded; its SI logR entry was excluded because of the GPA/cAMP duplication. The sensitivity analysis retaining the duplicate gives σ_δ = 0.320, ΔlogML = +123.8, δ_GIRK = −0.800 ± 0.049. The repair therefore works and names its own price: one additional coordinate per assay arm.

### S9.9 The load-bearing-coincidence statement

Definition: across cases, fracture positions coincide systematically with the positions on which the studies' quantitative conclusions load-bear. Per-case evidence: in D2R, the study's significant-bias reports concentrate on exactly the fractured ligands; in μOR, the therapeutic-window conclusion rests on the single fractured column; in the Atlas, one fifth of archived ligand pairs fail self-consistency and the fractured μOR column is absent by inclusion choice. The statement concerns the parameter foundation of the field's quantitative conclusions, not the biology of biased signalling itself, and it adjudicates no individual discrepancy. An interactive Atlas inconsistency map (per-group fracture flags, repair costs, and confidence labels) is deposited as an interactive supplement at `06_supplement_interactive/app/index.html`.

Source archive: `08_audit_trail/04_cell_line_4/_archive/results/report_cell_line4_GPCR_bias_fracture_audit_detailed.md`; codes 63–73, 77, 80 and their `_fulltext.txt` outputs; `08_audit_trail/04_cell_line_4/_archive/results/atlas_audit/` CSV tables.

---

## S10. Record trail

### S10.1 Artefact chain

Each adjudication produced the same artefact chain, each link a file: pre-registration document (hypotheses, clause list, decision lines, tolerances, seeds) → frozen decision table → scoring script (fixed seed SEED = 20260814 for P2/P3, 20260815 for P4-P7, bootstrap ×2000 unless otherwise stated) → verdict card (frozen at first execution) → errata where needed (appended, never overwriting). Verbatim samples of each link for P2–P7 are deposited with the archive; S7 cites them per adjudication.

### S10.2 Double-recorded failure inventory

The programme's failures are archived with the same discipline as its hits. The inventory includes: falsified adjudication clauses (P2-1, P2-2, P3-3, P4-1, P4-3, P5-1; S7); power-limited intermediates reported rather than resolved post hoc (P3-1, P4-2, P5-2, P6, P7); the two rejected NF-κB repair patches (codes 12, 15; S6.2); the ERK static-extrapolation falsification, which reversed an earlier claim of our own framework and is recorded in the version history (S6.3); three double-recorded chemotaxis pipeline errors (codes 55, 57, 59; S8.3); the PTH1R admission-ticket rejection, retained as a boundary-artefact teaching case (S9.1); and the framework's own bug lineage (six internal bugs and sixteen external-review items across cascade-framework versions v2.0–v4.2, plus Bug 7 and the chemotaxis analyst-normalisation Bug 8; S1.8, S1.9).

### S10.3 Version history as record trail

The absorbed NC manuscript lineage (v0.9–v0.11) carries dated change logs; the cascade-framework lineage (v2.0–v4.2) carries per-version bug and review-item appendices. Both lineages are unpublished elsewhere and are absorbed here in full (S11). The change logs are part of the deposit, so the path by which each claim entered, moved, or was retracted is inspectable.

### S10.4 Cloud re-execution protocol and results

Adjudications P2, P3 and P4 were re-executed on an independent cloud environment (2026-08-15) from the raw public datasets, with the frozen scoring scripts. Point estimates were bit-identical for P2 (all 21 promoter Δ values) and P4 (e.g., AUC_L = 0.5668941796791128) and verdicts were unchanged for all three. Records: `08_audit_trail/01_cell_line_1/results/P2P3_cloud_recheck_2026-08-15/`, `08_audit_trail/01_cell_line_1/results/P4_cloud_recheck_2026-08-15/`.

---

---

## S11. Cross-paper relationships (declarations)

**Assay-side preprint.** The companion assay-side account (TCS preprint^2^; ChemRxiv doi:10.26434/chemrxiv-2024-19rj6/v10) establishes, on the measurement side, that under multiplicative scale uncertainty only scale-free statistics survive transmission. The present paper inherits from it the definitions of the degeneracy group, the surviving-statistics set, and the payment theorem, and asks the complementary question on the cellular side. No results are shared between the two texts; the overlap is definitional and is cited.

**Companion machine-line manuscript.** The companion manuscript *Event coding in an uncalibrated world* (v16) develops the same stratification for artificial sensors and machine vision. Its scope boundary is the cell membrane: it does not claim results about intracellular signalling, and the present paper does not claim results about engineered systems. Three reverse predictions imported from the machine line are listed in S1.6 as prospective tests P-1/P-2/P-3.

**Absorption note.** Two internal lineages are absorbed into this paper and are not published elsewhere: the *Why cells count* manuscript lineage (v0.1–v0.11, the NC lineage) and the cascade-framework document lineage (v2.0–v4.2). All of their results that survive audit appear here, with their internal script numbers reconciled in S0.3.

**Prior-art positioning.** Two rounds of prior-art search were performed and archived. The first round (`08_audit_trail/05_programme_design/A3_paper_materials/prior_art_search_existing_work_vs_this_framework.md`) mapped the data-reconciliation, model-invalidation and identifiability literatures. The second, formal round (2026-09-23; `08_audit_trail/05_programme_design/paper_cellular_communication_audit_2026-09-23/prior_art_search_formal_round_2026-09-23/search_log_formal_round_v01.md` with per-query CSVs Q01–Q11) ran eleven structured queries on the Scholar index (top-20 relevant records each, title/author/year/abstract checked individually). Verdicts: no collision on any of the four novelty pillars: (i) exact degeneracy groups with reporting protocols for these four systems; (ii) pre-registered blind adjudication of scale-free-statistic sufficiency boundaries on public single-cell data; (iii) identity-level internal-consistency audit of a field database; (iv) the unifying claim that scale-free statistics are the currency of cellular communication. Known lineages (fold-change detection^16,47,57^; sloppiness and information-based sensitivity^32,63^; information optimisation^23,29–31,58–62^; dynamic encoding^12,64,65^; data reconciliation in engineering) are cited and bounded in the reference list. One residual item is flagged honestly: a WoS/Scopus pass requires institutional access not available at drafting time; the eleven query strings and result CSVs are archived so the pass can be re-run verbatim.

---

---

## S12. Extended Data figure inventory

- **Extended Data Fig. 1**: Framework schematic: layer model, degeneracy group action, survival set, three-leg criterion, information-loss identity.
- **Extended Data Fig. 2**: PdPC full panels: group action and Fisher spectrum; asymmetric-Michaelis triple; positive-feedback switch.
- **Extended Data Fig. 3**: p53 discrimination diagram with full signature grids (NF soft/spiky, EXC, canard note) and TrackC2 empirical validation.
- **Extended Data Fig. 4**: p53 channel-census six-panel audit figure (code 51r3).
- **Extended Data Fig. 5**: Dispersion signatures and noise-colour plane (code 51r4; code 52 TrackC2).
- **Extended Data Fig. 6**: Digital-limit theorem diagrams and conjecture decision tree.
- **Extended Data Fig. 7**: NF-κB dual-regime audit panels and the two rejected patches.
- **Extended Data Fig. 8**: Morris–Lecar duality instance table.
- **Extended Data Fig. 9**: Blind adjudications P2–P7: full unit tables and per-clause panels.
- **Extended Data Fig. 10**: Chemotaxis: full Pareto sweep, repair-path grid, Weber per-background residuals, prospective-arm panel.
- **Extended Data Fig. 11**: GPCR per-study fracture maps (D2R p_KA(arm, time); AT1R arms; μOR column table; code-80 extension).
- **Extended Data Fig. 12**: Atlas violation distributions, per-paper medians, constraint-projection repair map.
- **Extended Data Fig. 13**: Audit workflow: artefact chain and double-recorded failure ledger.
- **Extended Data Fig. 14**: Molecular-layer decomposition of the p53 timing constants (code 89): dominant Jacobian eigenvalue and period, phase-share attribution, K3/K4 closures.
- **Extended Data Fig. 15**: Independent closed-form cross-validation against the logistic-substitution DDE core (code 92): calibration reproduction, DDE integration check, state-consistency constraint.
- **Extended Data Fig. 16**: The cancer chain (P11 to P16): per-link verdict summary (regime switch; downstream misreading; channel capacity; resistance zero-direction walk; tissue-resolution boundary; frozen G13D prediction), colour-coded by verdict grade.

---

---

## S13. Evidence-strength register

*Added 2026-09-23, after the three-axis audit and after the mechanism-level stress test (code 82) described in S13.3. This register classifies every headline quantitative claim of the main text by the type of evidence beneath it, states the known soft spot of each claim, and points to the location where the defence lives. Its purpose is to make the attack surface of this paper explicit before a referee has to compute it.*

**Tier definitions.** Tier 1 (recompute-proof): exact mathematics verified numerically, or end-to-end recomputation from raw or entry-level data; survives hostile recomputation with the same inputs. Tier 2 (bounded): solid within the stated conventions; the soft spot (single dataset, digitised source, definition sensitivity, assumed dose conversion) is registered next to the claim. Tier 3 (calibrated wording): the direction of the claim is established, but the numerical magnitude is convention-dependent; the wording in the main text is bounded accordingly, or flagged here where a bound is still pending.

### S13.1 Register of headline claims

| ID | Claim (main-text location) | Tier | Basis | Soft spots, registered |
|---|---|---|---|---|
| E1 | PdPC four-dimensional degeneracy group; four parameter sets, one curve, max $|\Delta u| = 0$ (Fig. 2b) | 1 | Machine-precision algebraic identity of the model | A property of the PdPC model class, not of data; stated as such |
| E2 | Fisher conditioning $\chi \approx 30$ shape tier (Fig. 2c) | 1 | Direct eigendecomposition | – |
| E3 | Scale tier ${\sim}10^{300}$ (Fig. 2c) | 2 | Bounded extrapolation | Deepest directly computed point is $\chi \sim 10^{33}$; the $10^{300}$ figure is a labelled extrapolation (S2) |
| E4 | Hill-slope closed form $n_{\mathrm{H}} = 1 + 1/(2\kappa)$ (Fig. 2d) | 1 | Analytic; numerics agree to $6\times10^{-6}$ | – |
| E5 | NF Hopf onset, $A^2 \propto (s-s_H)$, $R^2 = 0.9926$ (Fig. 3a) | 1 | Continuation-grade numerics | – |
| E6 | Hysteresis width 0.275 vs normal-form 0.25 (Fig. 3b) | 1 | Direct measurement on the circuit | Binds only the subcritical escape hatch, not the main exclusion |
| E7 | Three-signature topology assignment, Theorem 4 (Results 3) | 1 | Thomas–Soulé theorem application | Genericity-based; subcritical exception handled by the hysteresis signature |
| E8 | p53 placed in excitable class; pure NF excluded (Fig. 3c) | 2 | Published statistics plus source-digitised slopes | $\alpha$ slopes come from digitised published figures, not raw traces; two $\sigma$ conventions (6.1 combined, 8.6 joint) both reported in S3.3 |
| E9 | Counting-law slope $\tau_r \approx 7.7$ h (Fig. 3e) | 2 | Digitised dose series | Same digitisation caveat as E8; shape corroborated at mechanism level (S13.3, M4) |
| E10 | Channel census 94.2% timing / 5.5% amplitude (Fig. 3d) | 3 | Fisher propagation on the seven-parameter excitable representative under a stated noise calibration | Mechanism-level recomputation (23-parameter Mönke model, code 82) shows the share swings 0.0–91.1% across defensible noise calibrations, nominal 45.1%; the number is a calibrated decomposition, not a structural invariant. Downgrade wording: see S13.4 |
| E11 | Two-layer death: CV(A) $0.002 \to 0.639$, CV(T) pinned 0.016 (Fig. 3d) | 3 | Ensemble statistics of the representative model | Direction corroborated at mechanism level; the pinned value 0.016 is a toy-ensemble statistic not directly remeasured on the mechanism model |
| E12 | Noise as exact null direction, $\lambda_7/\lambda_1 \leq 4.7\times10^{-19}$ (Results 4) | 1 | Exact zero direction, stable across three numerical configurations | A property of the first-moment statistics of the representative model; stated as such |
| E13 | Dispersion ratio: bare oscillator bounded 0.03–0.12; literature 5.8 requires an upstream stochastic sensing stage (Results 4) | 1 (bound), 2 (reproduction) | Analytic-style bound plus simulation | The reproduction point 5.7 at $\mu_d \approx T$ is an extrapolated point, labelled in S4.5. Independently corroborated at mechanism level: the Mönke core yields 0.01 (homogeneous) to 1.11 (Wip1-heterogeneous), below threshold (S13.3, M5) |
| E14 | Time-colour plane: measured $\mathrm{IPI}_{\mathrm{CV}} = 0.30$ implies fast intra-cellular stochasticity (Results 4) | 2 | Model-based inference | Classification, not a point estimate |
| E15 | Counting channel: $I(D;N) = 3.07$ of 3.17 bits (96.8%) (Fig. 3e) | 2 | Simulation under stated dose range and noise model | Range and noise conventions stated in S4; magnitude moves with them, conclusion (near-lossless counting) does not |
| E16 | Conditional digital-limit theorem, three predictions verified (S5) | 1 | Proof plus numerics on the circuit | Conditional form; uniqueness domain stated in S5 |
| E17 | Unconditional generalisation T′ and Theorem N (S5) | 1 | Proof | – |
| E18 | Digital-limit conjecture (S5) | – | Stated as a falsifiable conjecture, not a result | Correctly labelled |
| E19 | Chemotaxis rupture: no standard parameter set jointly consistent (Fig. 5a) | 2 | Pareto-front analysis of Moore et al.^22^, the best public single-cell dataset | Single dataset; the $N \approx 12$–$39$ lower bound 12 is soft (front segments with $R^2>0.9$ give 20–39), layered in S8; verdict is about a model class, not the bacterium (Limitations) |
| E20 | Weber line $K_{1/2} = 1.17\cdot(1.95+B)\,\mu$M, max deviation 0.10 dex across seven backgrounds (Fig. 5b) | 1 | Refit from raw data | – |
| E21 | Eight repair paths excluded (Fig. 5c) | 2 | Systematic patch enumeration within the model class | Exclusion is within the stated class; structural repairs outside it registered as open |
| E22 | Prospective mixed-ligand arm: scale-free midpoint hits $+0.014/+0.038/+0.031$ dex in the frozen $\pm0.15$ band; amplitude misses 0.105/0.126 vs 0.08 (Fig. 5e) | 1 | Pre-registered blind prediction on independent data | The strongest external evidence in the paper; failures recorded on the same card |
| E23 | Atlas-wide identity violations 21.5% (>0.3 dex) / 4.2% (>1 dex), 17,987 pairs (Fig. 6d) | 1 | Recomputed from entry-level table | Closure-triangle counts are pandas-version sensitive (3,167 vs 3,931); zero-violation conclusion invariant; quantifies inconsistency, not misconduct |
| E24 | D2R bifeprunox fracture $-3.5$ dex (Fig. 6a) | 2 | Case-level audit | Single-paper extraction |
| E25 | AT1R cross-arm fracture $-2.05$ dex, $5.0\sigma$ (Fig. 6b) | 2 | Case-level audit | Single-paper extraction |
| E26 | μOR GIRK systematic fracture and repair, $\Delta\log\mathrm{ML} = +128.2$, residual $0.811 \to 0.107$ dex (Fig. 6c) | 1 | code 80 end-to-end rerun | – |
| E27 | Failure concentration: 15 of 66 papers (Fig. 6d context) | 2 | Definition-dependent count | 15/66 stated in S9.5; 14/66 under the alternative pairing convention, stated here |
| E28 | Blind adjudications P2–P7, ledger 2/5/6/1 (Fig. 4) | 1 (procedure) | Pre-registered decision lines, independent datasets, falsifications recorded alongside hits | Individual clauses carry small-$n$ statistics; the P2 $k=4$ robustness arm is 6/7 (DCS2 at $-0.0102$ across the 0.01 line), reported as such everywhere |
| E29 | Cross-system consistency survey: Xenopus MAPK, ERK frequency modulation, NF-κB count encoding, tissue radiosensitivity (Discussion) | 2 | Literature synthesis | Qualitative consistency argument, not new measurement |
| E30 | Reverse predictions: three machine-side laws with frozen falsifiers (Discussion) | – | Frozen predictions, untested | Registered as open; the P2 adjudication has already contracted clause 1, stated in place |
| E31 | Hysteresis protocol proposal (Discussion) | – | Proposed experiment | A protocol, not a claim |
| E32 | Morris–Lecar duality check: event-time degradation 0.000% vs amplitude >70% under $0.3\times$–$3\times$ gain drift; amplitude leg fails at subcritical Type-II onset (Outlook) | 2 | Textbook working-point numerics | Domain bounded by bifurcation criticality, stated in place |

### S13.2 What the tiers add up to

Thirty-two headline entries: thirteen Tier 1, thirteen Tier 2, three Tier 3, three correctly labelled as conjecture or proposal. No claim requires retraction. Three claims (E10, E11, and the $\tau_r$ scale of E9) carry wording that must stay at the calibrated level; the required formulations are collected in S13.4.

### S13.3 Mechanism-level stress test (code 82, 2026-09-23)

After the three-axis audit, the published 23-parameter mechanistic model of the p53 pulse generator (Mönke et al.^10^, complete parameter table from the supplementary material, with the 2025 phase-resetting revision^66^ of the two mRNA production rates) was reconstructed and subjected to the audit stack. Results relevant to the register:

- M1. The mechanism model reproduces the archived signatures: sustained oscillation period 5.48 h (target 5.5 h), all-or-none pulses, amplitude insensitive to a 12-fold drive range (+0.9%), and the oscillation-to-sustained transition under high damage with the logarithmic but not the saturating DSB-sensing form. The topology verdict E8 is corroborated at mechanism level.
- M2. The Fisher spectrum over the 23 mechanistic parameters spans 10.5 orders of magnitude with 13 near-null directions; the stiffest direction is the Wip1 production axis (Tw, TW, P, dW), independently recovering the experimental conclusion that Wip1 sets the activation threshold.
- M3. The timing/amplitude information share is not invariant at mechanism level: 45.1% timing under the nominal calibration, 0.0–91.1% across the calibration grid. This drives the E10 downgrade (S13.4).
- M4. The induced-pulse counting law is hyperbolic at mechanism level ($R^2 = 0.995$ vs 0.492 linear), corroborating the shape of E9; the absolute dose scale depends on the assumed DSB-per-Gy conversion ($D_c$ differs by a factor 3.8 under the stated convention), so E9's slope is quoted as digitised-scale only.
- M5. The dispersion-ratio signature E13 is corroborated: the mechanism core produces Var$(t_1)$/Var$(T)$ of 0.01 (homogeneous) to 1.11 (Wip1 heterogeneity), below the literature value 5.8, supporting the claim that the first-pulse dispersion originates upstream of the oscillator core.
- M6. Two model-experiment boundaries registered: Wip1 heterogeneity at CV 0.5 produces amplitude CV 26.1% against the 10–15% target, and deep Wip1 knockdown locks the model into a sustained high state rather than increasing pulse counts as in siRNA experiments. Both are logged in the code 82 verdict card as mechanism-level boundaries, not paper-claim failures.
- M7. Post-audit route ledger (codes 83–88, 2026-09-25; verdict cards in `08_audit_trail/03_cell_line_3/results/`). (i) The literature dispersion ratio was pinned to its measurement context: variance-ratio convention, 5 Gy, MCF7 (see S4.5). (ii) A hazard-rate trigger driven by the repair dynamics (code 84) reproduces the counting law and recruitment but is falsified by the dose-independence of first-pulse timing (Lahav et al.^4^, constraint C1): hazard timing necessarily scales as 1/D. (iii) ATM gain heterogeneity anchored independently by the measured foci-count Fano factor (CV $\approx$ 0.2; code 85) closes the heterogeneity account but moves the dispersion ratio in the wrong direction, a falsified prediction preserved verbatim. (iv) A comparison reported on 2026-09-25 morning (code 86, "6.2 vs 5.8") was retracted the same day when dose verification exposed a convention mismatch (SD-ratio with within-cell IPI versus the archived variance ratio with pooled IPI); the retraction is recorded in the code 85–86 verdict card and stands as an audit-trail case. (v) The C1-consistent architecture (dose enters only a response gate and a persistence-termination channel; the clock is dose-independent; codes 87/88) meets C1–C4 simultaneously: $t_1$ flat at 2.5–3.0 h across 0.1–10 Gy, SD$(t_1) \approx 4.0$ h, pooled IPI SD $\approx$ 1.5 h, variance ratio 6.7–7.3 against the archived 5.8. Its residual failures are registered: sustained oscillation at 0.3 Gy (35% of cells) cannot be produced by any damage-driven termination (literature-internal tension), and the end-to-end dose information under the C1-consistent architecture is 0.55–0.63 bits per cell, revised down from the internal generative-model upper bound of 1.325 bits, part of which rested on dose-dependent timing that C1 excludes. The main-text claims (E8–E15) do not use that upper bound and are unaffected.
- M8. Molecular-layer decomposition of the emergent timing constants (codes 89–93, 2026-09-25; verdict cards in `08_audit_trail/03_cell_line_3/results/`). The four timing observables of the p53 pulse train were decomposed onto the kinetic constants of the reconstructed 23-parameter mechanism model: the oscillation period is set by the dominant Jacobian eigenvalue (5.95 h against the simulated 5.48 h, 8.5%), the interpulse interval and first-pulse delay close arithmetically against their measured bands, and the pulse width (FWHM 2.78 h) decomposes into a plateau segment plus an Mdm2-recovery-dominated collapse tail. An independent closed-form cross-validation against Belgacem's logistic-substitution DDE model (arXiv 2605.23722) reproduces his period to 0.07% and locates the state-consistency constraint that reconciles his parameter set with the measured 5.5 h. Full detail in S14.

### S13.4 Required calibrated formulations (wording ledger)

| Entry | Current strength | Required formulation |
|---|---|---|
| E10 | "carries 94.2% of its parameter information on timing channels" | The 94.2%/5.5% split is the decomposition of the seven-parameter representative under the stated calibration; at mechanism level the timing share remains dominant in direction but ranges with calibration (45.1% nominal). Main text, channel-census paragraph and Fig. 3d legend to carry the calibrated wording |
| E11 | "CV(T) stays pinned at 0.016" | Pinning is a representative-ensemble statistic; direction corroborated at mechanism level, exact value quoted with its model scope |
| E9 | "$\tau_r \approx 7.7$ h" | Digitised-scale slope; mechanism-level shape corroboration noted, absolute dose scale convention-dependent |

*All three formulations were applied to the main text on 2026-09-24 (abstract, channel-census paragraph, Fig. 3d and 3e legends, and the Results topology paragraph); the main text and this register are now consistent by construction.*

*The register is maintained as a living document; entries M1–M6 were added after the frozen audit of 2026-09-23 and are cross-referenced to the code 82 verdict card and report (`08_audit_trail/03_cell_line_3/results/verdict_card_code82_p53_mechanism_model_2026-09-23.md`). Entries M7–M8 were added on 2026-09-25 and are cross-referenced to the verdict cards of codes 83–93.*

## S14. Molecular-layer decomposition of the emergent timing constants (codes 89–93, 2026-09-25)

This section closes the p53 line at the molecular layer. The four timing observables used throughout the main text (oscillation period $T$, interpulse interval IPI, first-pulse delay $t_1$, pulse width $w$) are decomposed onto the kinetic constants of the reconstructed 23-parameter mechanism model (Mönke et al.^10^ calibration with the 2025 revision^66^; modelling lineage^13,48,49,67,68^), each decomposition is checked against its measured band, and the residual boundaries are registered. All verdict cards and result JSONs are in `08_audit_trail/03_cell_line_3/results/`; all scripts are deterministic (seed 20260925).

### S14.1 Constant register

Each kinetic constant entering the decompositions was classified by provenance: K (measured directly in the cited source), M (model-calibrated by Mönke et al.^10^), D (derived from K and M entries by the model's own steady-state or phase constraints), G (grouped, identifiable only as a combination). The register shows that every constant entering the closed timing expressions below is either K or D anchored on K; no timing conclusion rests on a free fitting parameter of the present work.

### S14.2 Decomposition of the four timing observables (codes 89–91)

**K1, period.** Linearising the six-species model at its operating point gives a dominant Jacobian eigenvalue $0.627 + 1.057i$, i.e. $T = 2\pi/1.057 = 5.95$ h against the simulated 5.48 h (8.5% deviation, registered as boundary T5: the linear prediction inherits the model's own width bias, see S14.5). Phase-share attribution over the oscillation cycle: p53 30%, ATM* 24%, mdm2 mRNA 18%, wip1 mRNA 15%, Wip1 protein 10%, Mdm2 protein 2%. Prior art is acknowledged: oscillation-period direction from molecular constants was established for simplified cores by Monk^13^ and by Wang et al.^14^, and a closed-form period for a logistic-substitution DDE core was given by Belgacem^15^. The contribution here is the complete six-species phase-share attribution and the audit stack around it, not the direction itself.

**K3, interpulse interval.** Phase-space decomposition gives IPI $1.70$ h (phase estimate) / $2.00$ h (naive constant-sum estimate) against the measured $2.0 \pm 0.5$ h.

**K4, first-pulse delay.** $t_1 = 1.13$ h (phase) / $1.27$ h (naive) against the measured $1.25$ h anchor.

K3 and K4 are registered as *consistency closures*, not independent predictions: the model-calibrated constants used in the sums were themselves anchored on the same measured bands, so the closure is partly circular by construction (audit item A2). The non-circular content is that no inconsistency exists between the molecular constants and the timing observables at the stated precision.

**K2, pulse width (code 91).** Simulated FWHM of $2.78$ h, decomposing into a plateau segment of 0.79 h plus a collapse tail of 1.99 h dominated by Mdm2 recovery. The width is dose-independent (CV of $0.055$ across the dose grid). Perturbation boundaries registered: Wip1 at 85% widens pulses by +30% (against the invariance wording in Mönke et al. S9, registered as T6), and Wip1 at 70% locks the model into the sustained state.

### S14.3 Full-cascade analysis and first-pulse timing budget (code 90)

Fourteen loops of the DSB $\to$ ATM $\to$ p53 $\to$ Mdm2/Wip1 cascade were audited for closure between the molecular constants and the cascade-level observables: 8 closed, 4 in registered tension, 3 gaps at the time of audit; gap L13 was subsequently closed by code 91, and L12 belongs to the same family as boundary T6. The first-pulse delay decomposes as a sensing segment of 0.83 h ($\gamma$H2AX formation 0.5 h + ATM activation 0.19 h + p53 threshold crossing 0.14 h) plus a ramp segment of 2.38 h, totalling 3.21 h against the measured 2.5–3.0 h band: a borderline closure, registered as such. The temporal ordering of cascade events is consistent with the measured sequence.

### S14.4 Independent closed-form cross-validation (code 92)

Belgacem's logistic-substitution DDE model of the p53–Mdm2 core^15^ was reproduced and used as an independent check. (i) His closed-form quantities reproduce exactly under his stated calibration ($AB = 1.7204$ against his 1.72; $T = 5.642$ h against his stated $\approx 5.6$ h). (ii) An independent Python DDE integration of his equations at $\tau = 1.0$ h gives $T = 5.646$ h against his closed form 5.642 h (0.07%). (iii) State-consistency constraint: his basal-state half-lives combined with the measured 2 h delay yield $T = 8.40$ h, failing against experiment; his stressed-state parameter set (p53 degradation $\gamma = \ln 2 / 0.25$ h, Mdm2 $d_M = 2$/h, $\tau = 2.0$ h) yields $T = 5.60$ h, a three-way agreement between the measured 5.5 h, his closed form 5.60 h, and our six-species simulation 5.48 h. The constraint is registered in the constant register (section F): the oscillatory regime parameters must be read as stressed-state quantities, not basal-state ones.

### S14.5 Registered boundaries and literature-internal conflicts

T5. The linear-eigenvalue period (5.95 h) overshoots the simulation (5.48 h) by 8.5%, and the model width (2.78 h) sits below the archived 3.5 h band; both directions are consistent with the mechanism core running slightly fast and narrow under the nominal calibration.

T6. Wip1 perturbation widens pulses by +30% at 85% and locks the system at 70%, whereas Mönke et al. S9 reports approximate width invariance; this is a literature-internal conflict between Batchelor et al. Fig. 6/7 and Mönke et al. S9E/F (audit item A6), registered rather than resolved.

L5 closure attempt (2026-09-27, negative). The $D_c$ tension (archived 0.261 Gy vs mechanism-level 0.99 Gy, a factor 3.8) is conducted entirely by the single conversion assumption of L0, 35 DSB per Gy. Surveying the experimental literature on DSB yield per Gy in mammalian cells bounds the plausible conversion at 20–40 DSB per Gy per cell: 35 per Gy from Rothkamm and Löbrich's γH2AX calibration, 30–36 per Gy from PFGE-based assays on G1 diploid cells, and about 19–20 per Gy from live-cell 53BP1 focus counting. Closing the tension would require roughly 9 DSB per Gy, far outside this band. L5 is therefore not closable by the conversion assumption and is promoted from registered tension to established tension: the counting-law slope discrepancy between the digitised archive scale and the mechanism-level scale is real, and E9's slope continues to be quoted as digitised-scale only.

### S14.6 Structural correspondence with the scale-invariance framework (code 93)

Three correspondences between Belgacem's^15^ logistic substitution and the scale-invariance framework of this paper were verified at machine precision or better. (i) Hill-to-logistic identity: $\mathrm{Hill}_n(x) = x^n/(x^n + \theta^n)$ is exactly $\mathrm{logistic}_n(\ln x)$; the logit linearisation is the same algebraic operation as the log-map lemma of the S5.1 proof sketch (Lemma DL1; max error $2.2\times10^{-16}$ over three parameter sets). (ii) Concentration-scale freedom of the clock: under the joint rescaling $(\kappa, \theta, x) \to s(\kappa, \theta, x)$, $\lambda \to \lambda/s$, the characteristic frequency $\omega_c$ and closed-form period are invariant to machine precision while the equilibrium occupancy $f(x^*)$ is scale-identical; the p53 clock pins time-constant ratios, not absolute concentrations. (iii) Event-time functionals survive arbitrary monotone amplitude rescaling: applying $y = x^{1.7} + 0.3\sqrt{x}$ to the simulated p53 trace shifts all seven pulse-peak times by exactly $0.00$ h and leaves the period at 5.480166666666666 h, while the amplitude CV moves from 0.308 to 0.438. Amplitude statistics are a coordinate convention; event times are the physical content. This is the same separation as the main-text timing-versus-amplitude channel census, here shown to be a mathematical property of the dynamics rather than a property of the estimator.

## S15. Morphogen-gradient scaling analysis (Bicoid to gap genes, 2026-09-27)

**Data and provenance.** Per-boundary precision and scaling tables of Morton de Lachapelle and Bergmann^71^ (PMC2858443, supplementary Datasets S1 and S2; 84 boundaries across Kr, Gt, Hb and Eve stripes 1–7, at 1$\times$, 2$\times$ and 4$\times$ bcd dosage, more than $150$ embryos) were retrieved from the publisher archive on 2026-09-27 and parsed verbatim to CSV; row-level correspondence between the two tables was asserted before analysis (analysis script and archived CSVs in the project repository, morphogen line).

**Pre-registered checks and outcomes.** Four checks were registered before opening the scaling table. P1 (anterior hyper-scaling at 1$\times$ bcd): confirmed, median scaling coefficient $S = 1.728$ over the ten anterior boundaries (below $35\%$ EL), all ten point estimates above 1.58, five with 68% CI excluding 1. P2 (mid-embryo scaling near 1): essentially confirmed, median $S = 1.175$ at 1$\times$ and $1.067$ at 2$\times$, with 53–59% of mid-embryo boundaries compatible with $S = 1$ within CI. P3 (monotone dose dependence): confirmed globally, CI-weighted mean $S = 1.069, 0.955, 0.841$ at 1$\times$, 2$\times$, 4$\times$; at 4$\times$ the mid-embryo domain turns hypo-scaling (median 0.782, 0/22 hyper-scaling). P4 (precision does not collapse posteriorly): confirmed, $\sigma(x/L)$ median 2.11% EL over all 84 boundaries (range 0.93–4.09%), with no positional breakdown for the posterior Eve stripes.

**Interpretation within the framework.** Two deviations from perfect scaling are exactly where the physics says they should be. Because boundary position is read off threshold crossings of the Bcd concentration profile rather than off geometry directly, the scaling coefficient inherits the profile's own scaling: the dose dependence of $S$ (P3) is the signature of concentration-threshold readout. Among the four candidate mechanisms compared in Ref.^71^, the observed $S(x)$ shape (anterior hyper-scaling, mid-embryo near-perfect scaling, mild posterior hypo-scaling) selects the nuclear-density degradation model, in which Bcd degradation tracks nuclear density $\propto N/L^2$. Nuclear density is a geometric counting object: counting nuclei measures embryo length without any maintained reference, so length normalisation is supplied by a class-I carrier and the readout is calibration-free over the mid-embryo.

**Registered limitations and open items.** (i) Dataset S1 reports aggregate $\sigma(x/L)$ only; without per-embryo lengths the absolute (µm) precision arm could not be computed, and the scaling coefficient was used as the decisive statistic instead, as pre-registered. (ii) The mild posterior hypo-scaling (median $S = 0.889$ beyond 56% EL) is consistent with a terminal-system anchor that does not scale; candidate paid reference, pending verification against the original discussion. (iii) One point anomaly registered: Hb boundary at 2$\times$ bcd has $S = 1.225$ (CI 0.115), above its 1$\times$ value 1.020; unresolved. (iv) Per-embryo datasets (Gregor et al. 2007; Petkova et al. 2019) are the designated follow-up source for the absolute-precision arm.

### S15b. Per-embryo Bcd anchor test (P9; Liu et al. 2013 profiles via the Nikolic et al. repository, 2026-09-27)

The absolute-precision arm designated in S15(iv) was executed as pre-registered adjudication P9 (preregistration frozen before any data values were inspected; script `p9_bcd_anchor_audit.py`, SEED = 20260927). Per-embryo exponential length constants were fitted in absolute coordinates ($x = x_s \cdot L$, window $x/L \in [0.10, 0.50]$, log-linear, 581/582 valid embryos, median fit $R^2 = 0.985$) on the 2XA reference-line profiles of Liu et al.^73^ distributed by the scale-invariance repository of Nikolic et al.^74^.

**P9-1 (main clause): hit.** The scaling exponent of the length constant, $\lambda \propto L^{\alpha}$, is $\alpha = -0.182$ with 95% bootstrap CI $[-0.476, +0.109]$ (embryo-level resampling $\times 2000$; Theil–Sen arm $-0.199$), against frozen lines (anchor: point $\leq 0.7$ and CI upper below $1$; perfect scaling: point $\geq 0.9$). The length constant is absolute within the species: natural $\pm 10\%$ length variation is absorbed without retuning, and perfect scaling is excluded at ${\sim}7\sigma$. The median $\lambda = 87.6$ µm sits in the literature band 80–100 µm, a pipeline-correctness control. This is the per-embryo, absolute-coordinate, explicit-exponent counterpart of the information-theoretic non-scaling statement of Ref.^74^.

**P9-2 (control clause): consistent at low precision, extraction-sensitive.** The Hb anterior boundary scaling exponent from raw profiles returned $\alpha_b = -0.64$ with a wide CI (SE 0.50; power-limited at $\pm 10\%$ length spread), formally triggering the registered divergence rule. The designated power upgrade on the time-corrected N301 subset returned $\alpha_b = 1.94$, CI $[0.43, 3.41]$, covering 1 while P9-1 excludes 1: per the frozen rule the "profile anchored, readout rescaled by the network" reading stands, but the two extractions bracket the estimate and precision is low; the verdict is registered as consistent, low-precision, extraction-sensitive.

**P9-3 (exploratory, no line).** Anterior amplitude $B_0$ scales weakly if at all ($\gamma = 0.38$, CI $[-0.19, 0.94]$).

**Attribution statement (three layers, per audit discipline).** (i) Literature facts: anterior localisation of bcd mRNA is an active, maintained process (dynein-mediated transport plus anterior anchoring, with continual active transport in late oogenesis), and the terminal (Torso) system is maternally deposited with pole-localised activation; both are textbook-established. (ii) This paper's measurements: the absolute length constant (P9-1), the dose-dependent scaling-coefficient profile with anterior hyper- and posterior hypo-scaling (S15), and the network-level near-perfect scaling.^74^ (iii) Framework interpretation: the two maternal systems are read as the class-III paid references of the axis, their positions legible in the scaling deviations; this reading is ours and is falsifiable: disruption of the anterior anchoring machinery (exu/swa class) is predicted to move the anterior scaling behaviour selectively, sparing the mid-embryo. Qualitative anterior deletions in those mutants are compatible; the scaling-signature test has not been performed.

**Registered tensions.** Cross-species $\lambda \propto L$ over a five-fold length range (Gregor et al. 2005) and partial within-species retuning in abnormally large embryos ($\lambda$ 102.6 to 141.8 µm; Liu et al. 2013) define a two-regime picture: the anchor is absolute under natural within-line variation and retunable across evolutionary or extreme-size regimes. The large-embryo per-embryo data required to resolve the quantitative tension are not publicly available (reconnaissance 2026-09-27; related datasets are author-on-request); registered as the decisive follow-up.

**Cross-regime reconciliation (verified against primary sources, 2026-09-27).** (i) Across species, the relative length constant is a structural constant: in Gregor et al.^75^ the per-species distributions of $\lambda/L$ (27 L. sericata, 35 D. melanogaster, 18 D. busckii embryos) have means agreeing within 2%, while dextran diffusion constants vary only slightly and the nuclear count is fixed ($\log_2 N_{\mathrm{nuc}} = 12.8 \pm 0.2$). The same source proposes the mechanism the framework expects: if degradation occurs dominantly within nuclei, the effective lifetime tracks nuclear density, which scales with embryo size because the nuclear count is fixed; the inferred lifetime retuning (about 3 min in D. busckii to 32 min in L. sericata) is therefore supplied by counting geometry, not by per-species molecular recalibration. (ii) Within species, the mean-level check on the size-selected inbred lines of Cheung et al.^76^ is directionally consistent: with $L$ ratio $1.236 \pm 0.068$, the nuclear-density prediction $\lambda \propto L$ gives 126.8 µm against the observed $141.8 \pm 19.1$ µm ($0.78\sigma$), and the reported $\lambda/L$ values are indistinguishable ($0.22 \pm 0.03$ vs $0.22 \pm 0.05$). The residual $+12\%$ is in the direction of the documented confound in that line (an abnormally broadened bcd mRNA source, threefold signal area, not generic to large embryos), and a $\lambda \propto L^2$ alternative sits equally close ($-0.78\sigma$); the check is registered as directional, mean-level only, confound on record. (iii) The two regimes are thereby compatible with one statement: molecular constants (diffusivity, per-nucleus degradation) are fixed; what rescales the anchor across size regimes is nuclear-counting geometry, a class-I carrier.

## S16. Wnt/β-catenin fold-change analysis (Goentoro 2009, 2026-09-27)

**Data and provenance.** Goentoro and Kirschner^72^ (Mol. Cell 36, 872–881; PMC2921914) provide no machine-readable supplementary data; the audit digitised Fig. 4F (fold-change of β-catenin signalling versus LiCl dose, normalised to control) from the publisher figure. Green marker centroids were extracted programmatically from the 468×599 archive image (panel upscaled 4×) and mapped through a two-point log-linear calibration (y: 1 and 0.1 gridlines; x: 0 and 80 mM ticks). Registered precision: ±0.05 relative fold-change, ±2 mM; sufficient for plateau and breakdown claims, not for sub-10% statements. Digitised values are archived in the project repository (Wnt line, `fig4F_digitized.csv`; verdict card `verdict_Goentoro2009_foldchange_2026-09-27.md`).

**Pre-registered checks and outcomes.** P1 (plateau): all four plateau points at ≤50 mM LiCl fall within ±5% of unity (0.96, 1.04, 1.02, 0.98), passing the ±20% band. P2 (breakdown): fold-change falls from 0.98 at 49 mM to 0.61 at 58 mM and to 0.25–0.27 at 70–82 mM; the breakdown onset localises to 50–58 mM, compatible within digitisation precision with the original qualitative statement. Corroborating panels read at lower precision: Axin1 overexpression leaves fold-change at ≈1 (Fig. 4D); 80 mM LiCl reduces the Wnt3a response fold-change from ≈6 to ≈1.5 (Fig. 4E); the same dose-shape holds in β-catenin-overexpressing cells (Fig. 4H).

**Interpretation within the framework.** The fold-change (FCD) readout is confirmed quantitatively: cells read stimulus fold over background, not absolute level, across a plateau spanning at least 0–50 mM LiCl. The normalisation mechanism identified in the source literature, balance of synthesis and destruction-complex turnover, is a synthesis–degradation ratio of class-I type and requires no maintained paid reference. Wnt/β-catenin is thereby added to the ratio-law camp as its third system (main-text Table 6).

**Registered limitations.** Single-source, single readout arm (293T reporter arm of one study); high-resolution figure retrieval was blocked by the publisher's anti-crawl challenge and the audit used the standard-resolution archive image; independent cross-laboratory replication of the FCD plateau is the designated follow-up.

## S17. Two-stage cascade analysis: M3R-ACh to Ca2+ (Keshelava et al. 2018, 2026-09-27)

**Data and provenance.** Supplementary Data 3 of Keshelava et al.^20^ (PMC5830429; `peak_values.csv`, per-cell Fura-2 ratio peaks, baseline-subtracted by the authors' pipeline) was parsed verbatim: 27 experiments, 433 cells each stimulated at seven acetylcholine concentrations with five repetitions; 353 cells with complete 7x5 matrices entered analysis. Row-to-dose mapping was verified empirically: the population mean peak declines monotonically across rows (0.215 to 0.018 ratio units), the last row sits at baseline consistent with the paper's statement that single cells typically do not respond to the lowest ACh concentrations, and within-row repetitions drift upward (ratio 1.03-1.22), excluding an adaptation artefact. Absolute dose values follow the Methods ladder [10, 3, 1.5, 0.75, 0.5, 0.25, 0.1] uM; a provenance discrepancy is registered: the authors' Supplementary Data 2 MATLAB script carries a legacy ladder [4, 2, 1, 0.5, 0.25, 0.13, 0.06]. All statistics reported here except absolute EC50 values are invariant under monotone relabelling. Analysis script and outputs: project repository, cascade line (`cascade_audit_gpcr_ca.py`, `cascade_audit_results.json`).

**A. Population-level dose-response (fracture propagation).** Per-experiment Hill fits (26 experiments passing): median top 0.202 ratio units (CV 0.24), median EC50 0.63 uM with CV 1.16 and range 0.16-5.17 uM (32-fold spread, no bound collisions), median Hill coefficient 1.84 (CV 0.34). The parameter drift documented at the receptor-pharmacology layer (S9) propagates into the downstream population dose-response.

**B. Single-cell precision.** Median within-cell coefficient of variation across the five identical repetitions: 3.4-3.7% at the two highest doses, 5.9% at 1.5 uM, rising to about 20% at the three lowest doses (near the detection threshold). Across-cell CV of the per-cell mean: 28-36% at high dose, exceeding 50% at low dose. Per-cell rank preservation across adjacent doses: Spearman rho 0.90, 0.86, 0.85, 0.87, 0.73, 0.58 descending the ladder. The fracture is confined to the across-cell gain distribution; the intracellular measurement is precise.

**C. Carrier adjudication.** Adjacent-dose discrimination AUC from raw amplitude (353 cells pooled): 0.556, 0.573, 0.665, 0.577, 0.733, 0.705. Self-normalisation to each cell's own maximal-dose response gives 0.572, 0.547, 0.698, 0.579, 0.735, 0.685, no significant improvement; this negative result is registered with its interpretation: the five-repetition averaging already compresses within-cell noise to 3.5%, so the decoding bottleneck is the across-cell gain distribution, not intracellular noise. These values are quantitatively compatible with the source paper's channel-capacity estimate of 2.06 +/- 0.31 bits (lower bound 1.65 +/- 0.18).

**D. Event-presence carrier.** Fraction of cells with mean peak above 0.02 ratio units, descending the ladder: 1.00, 1.00, 0.98, 0.88, 0.83, 0.57, 0.24. At the low-dose end, dose information is carried by whether the cell responds at all.

**Interpretation.** Because a cell cannot know its own position in the across-cell gain distribution, absolute amplitude is not a reportable statistic downstream of the first stage; self-referenced statistics and event presence are what survive. This is the cascade-level instance of the free-first principle, complementary to the p53 cascade in which amplitude is discarded already at the sensing stage (S3, S4, S14). Verdict card: `verdict_GPCR_Ca_cascade_Keshelava2018_2026-09-27.md`.

**Registered limitations.** Single ligand (ACh), single receptor (M3R), single cell system; a second cascade instance is registered as follow-up (candidate: the three-stage chemotaxis chain, receptor cluster to CheY to motor, archive in hand). The legacy-ladder discrepancy noted above moves no conclusion except absolute EC50 quoting.

### S17b. Chemotaxis cascade relay: baseline-restoration analysis (Moore 2024 archive, 2026-09-27)

**Data and provenance.** The Moore et al. FRET archive (Dryad doi:10.5061/dryad.nvx0k6dzz; 48 .mat files) stores per-cell data as 35 events x 20 samples: each event spans one 40 s attractant pulse onset, with 10 samples (5 s) of pre-onset OFF baseline and 10 samples (5 s) of post-onset response. The intra-pulse adaptation plateau is not part of the export (per-cell traces carry a 35 s gap over the sustained-pulse period); what the archive supports is a history-dependence test of baseline restoration across the 5 s washout between pulses. 2,217 cells and 77,595 pulse events passed structural checks. One engineering trap is registered: the stimulus channel is stored as unsigned 8-bit integers, so downward stimulus transitions wrap around in naive differencing; the corrected pipeline is `chemotaxis_adaptation_audit.py` (cascade line, project repository).

**Result.** The pre-pulse baseline median varies over only 0.056 activity units across previous-pulse doses from 0 to 180 uM, with no monotone trend (0.485 after dose 0; 0.436 after 80 uM); the baseline difference between pulses following high (40/80 uM) and low (2/10 uM) predecessors is -0.0135, i.e. 2.2% of the median response amplitude (0.61). Baseline restoration by the methylation adaptation relay is therefore complete within the washout window and independent of dose history across a 90-fold range.

**Interpretation.** In the three-stage chemotaxis cascade, the parameter rupture lives at the receptor-cluster layer (S8), the adaptation relay resets the operating point with 2.2% residual (this section), and the motor output is a digital fold-change switch; the system's single paid reference (methylation-maintained effective affinity) sits at exactly the middle stage that must absorb the upstream rupture. Registered limitation: the measurement is baseline restoration after a 5 s washout, not the intra-pulse adaptation curve, which requires the original long trajectories.

### S17c. Cross-system confinement test (P10; ERK and NF-kB pulse trains, 2026-09-27)

Whether the gain-distribution confinement of S17 generalises was tested under a frozen preregistration (P10; identical pipeline on both systems: per-cell baseline at the 20th percentile, peaks at prominence 0.10 and minimum distance 8 frames, cells with $\geq 3$ pulses; confinement ratio $R$ = across-cell CV of per-cell mean pulse amplitude over median within-cell CV; split-half rank preservation; bootstrap $\times 2000$, SEED = 20260927; frozen lines: confinement $R \geq 3$, falsification $R \leq 1.5$; rank $\rho \geq 0.5$, falsification below $0.2$).

**Verdict: both clauses falsified, the generalisation is withdrawn to the M3R case.** ERK-KTR pulse trains under sustained ligand (Chavez-Abiega et al.^19^; histamine and UK14304 arms) returned $R = 0.59$ and $1.08$; NF-$\kappa$B sequential-stimulation trains (Wang et al.^21^; 11,265 cells) returned $R = 1.20$ (CI [0.77, 1.61]) with split-half $\rho = -0.40$. Two structural registrations accompany the falsification. First, the NF-$\kappa$B design is a memory paradigm in which the second stimulus is attenuated by design, not an identical-repetition design; the verdict stands but the regime differs from S17. Second, the Akt readout channel produced no qualifying pulses at all: it is a sustained, non-pulsing channel, itself carrier-consistent. A supplementary S1P-ERK arm returned $R = 2.17$ (intermediate by the frozen bands). The registered interpretation, marked post hoc: amplitude in ERK/NF-$\kappa$B dies by within-cell pulse-to-pulse stochasticity rather than by across-cell gain dispersion; both death modes deny amplitude an anchor, and the carrier-level conclusion (frequency and count codes; Table 6) is unchanged. Main-text wording is unaffected: the cascade section claimed the rule for the three audited cascades only.

## S18. Complete ledger of spin-off discoveries (2026-09-27)

Eleven findings produced by the audits that are independent of the main theorem. Each is stated with its evidence anchor and status. Main-text Table 7 presents the five sharpest; this ledger is the complete record.

**Methodology level.**

### S18.1 Repair timescale from the counting-law slope
Lemma 3's slope inverts the DSB repair timescale: $\tau_r \approx \nu \cdot T = 1.4 \times 5.5\ \mathrm{h} \approx 7.7\ \mathrm{h}$, consistent with the fast DSB-repair component. Any pulse-count dataset thereby yields a repair rate without a $\gamma$H2AX assay. Anchor: S14.

### S18.2 A computable topology diagnostic
The dispersion ratio $t_1/\mathrm{IPI}$ is bounded at 0.03-0.12 for any bare oscillator; the literature value 5.8, reproduced at 5.7, certifies an independent stochastic sensing stage upstream of the oscillator core, and its magnitude inverts the upstream delay dispersion. Anchor: S4, S14.

### S18.3 The hysteresis protocol
The second signature of the discrimination scheme has never been measured for p53 because all published dose series are one-directional scans; the two-way scan protocol is frozen with readout criteria, power analysis and registered risks. Anchor: S3.5.

### S18.4 Analysis-with-calibration method
The GPCR meta-audit measured its own false-positive rate on synthetic null batteries before issuing verdicts (code 72), survived common-mode calibration and jackknife (code 73), and ships a constraint-projection repair map (code 77) plus a one-coordinate model extension that resolves the muOR fracture outright (code 80). Anchor: S9.

**Mechanism level.**

### S18.5 Rupture confined to the gain distribution
Across-experiment EC50 scatter in the M3R-to-Ca2+ cascade (32-fold) decomposes into across-cell gain dispersion against 3.4-3.7% within-cell repeatability and rank preservation rho 0.58-0.90 across the full ladder (0.85-0.90 above the lowest two doses): the population-level parameter drift and the single-cell heterogeneity literature are two views of one object. Anchor: S17.

### S18.6 Paid-reference positioning rule
The chemotaxis chain's single paid reference (methylation-maintained effective affinity) sits at exactly the relay that absorbs the upstream rupture; baseline restoration measured at 2.2% residual, independent of dose history over 90-fold. Anchor: S17b, S8.

### S18.7 The chemotaxis rupture itself
No standard MWC-plus-perfect-adaptation parameter set is simultaneously consistent with the amplitude table and the per-cell midpoint distribution of Moore et al.; eight minimal repair paths excluded; the minimal form of the active-retuning candidate excluded twice over. Anchor: S8, S8.7(i).

### S18.8 Wnt breakdown boundary quantified
The fold-change plateau of beta-catenin signalling holds within +/-5% across 0-50 mM LiCl and breaks down at 50-58 mM, localising the boundary the source study stated qualitatively. Anchor: S16.

### S18.9 Nuclear density as a free geometric counter
In the Drosophila gap-gene readout, the observed scaling-coefficient profile selects the nuclear-density degradation model: length normalisation supplied by counting nuclei, a class-I carrier with zero maintained reference. Anchor: S15.

### S18.10 The p53 clock pins time-constant ratios, not absolute concentrations
Under joint rescaling the characteristic frequency and closed-form period are machine-precision invariant; event-time functionals survive arbitrary monotone amplitude rescaling exactly. Anchor: S14.6.

**Field level.**

### S18.11 Database-level identity failures
21.5% of 17,987 archived ligand pairs in the Biased Signaling Atlas violate the operational model's cell-level identity by more than 0.3 dex (4.2% by more than 1 dex); 15 of 66 multi-ligand papers exceed 0.3 dex median deviation; inclusion choices silently dropped the fractured muOR column, changing the evidence base of a published conclusion. Anchor: S9.


## S19. The cancer chain (P11-P16): pre-registration, verdicts and registered limitations

This section is the SI home of the Results subsection "Hijacked, not interrupted". Every link was pre-registered before data access; cards and code are archived in the project ledger (tumour-prediction line). Verdicts below quote the frozen cards verbatim in substance.

### S19.1 P11: the p53 counting law across cell lines

Data: digitised Stewart-Ornstein & Lahav 2017^77^ (Figs. 3F/3G/5F/6B) and Finzel 2021^82^ (Fig. 30). Verdicts: P11-1 falsified with direction reversal (rho = +0.857, one-sided p = 0.0068, seven oscillatory lines): slow repair does not speed counting; it collapses the oscillatory regime into a sustained one. The counting law is contracted to the oscillatory regime. P11-2 (ATM IC50 stratifies dose-dependent pulse broadening): undecidable at 6/12 concordance against a frozen 40-75% indecision band. P11-3 (non-transformed MCF10A tighter than cancer lines): weakly falsified (MCF10A CV 0.57 versus cancer median 0.50; asymmetric sample sizes registered). Open items: raw single-cell trajectories of ref. 77 are not public (G1); a second non-transformed control line is needed (G2); cross-lab absolute FWHM values are not comparable (G3).

### S19.2 P12: gain distributions and the resistance zero direction

P12a (mutations widen the across-cell ERK gain distribution): falsified; baseline CV 0.339/0.353/0.353 across three lines, |delta CV| <= 0.014. P12b (clinical resistance mutations preserve catalytic competence): confirmed on the ABL axis (clinical n = 5 versus predicted-not-clinical n = 8, median |delta log10(kcat/KM)| 0.237 versus 0.705, fold ratio 2.94, one-sided Mann-Whitney exact p = 0.0326). EGFR axis (saturation mutagenesis, GSE305057): our self-built enrichment classification was discarded after registration because low-count denominators inflate enrichment ratios (registered amendment A5); on author-validated variant lists, resistance variants sit above the missense activity background (median ratio 2.2, one-sided p = 0.002, n = 17), with a registered counterexample subset (L718Q/V near the nonsense floor; C797S at 0.003): the resistance manifold splits into an activity-preserving branch and a binding-disrupting branch that appears clinically only on compensatory backgrounds. BRAF third axis: undecidable, no same-study dual measurement found (registered).

### S19.3 P13: ligand-to-ERK channel capacity under oncogenic Ras/Raf

Data: Gillies et al. 2020^78^ Source Data Fig. 2 (EKAR3, 8 isogenic MEF lines, 6 ligands, 4 dose levels including ligand-free; about 400 cells per condition). Locked output: 0-30 min post-stimulus peak minus pre-stimulus median. Estimator: equiprobable 8-bin discretisation per line, plug-in MI with Miller-Madow correction, 2,000 paired cell-level bootstraps. Main arm (EGF), delta I versus KRASWT (0.426 bits): G12C -0.201 [-0.259, -0.143]; G12D -0.228 [-0.282, -0.184]; G12V -0.050 [-0.109, +0.001]; Q61R -0.335 [-0.387, -0.292]; BRAF-V600E -0.413 [-0.458, -0.378]. Cross-ligand consistency: Q61R and BRAF-V600E negative in 6/6 ligands; G12V direction mixed (significantly positive on PDGF and IGF). Robustness: deltas stable across 4/6/8/12/16 bins and under a continuous KSG estimator (300 bootstraps). Registered limitations: EKAR3 FRET rather than KTR readout (amendment A1); four dose levels cap MI at 2 bits, so absolute compression is underestimated; MEF is not human tumour tissue. Verdict: confirmed for strong alleles; contracted to an allele-strength-graded claim.

### S19.4 P14: downstream decoding of the p53 regime switch

Data: Jimenez-Asins et al. 2022^79^ public repository (MCF7; oscillatory 10 Gy X-ray versus rising Nutlin-3 arms; RNA-seq TPM, 0-9 h hourly plus 12/24/48 h, two replicates per arm). Our independent audit of the archived TP53 protein integrals (TMT, 0-9 h) found the rising arm carries 1.65 times the pulsatile arm's exposure (20.78/19.12 versus 36.70/29.18); the AUC correction is a registered log-linear assumption. At t = 6-9 h on 9,103 induced genes (max fold change >= 2 in either arm): 29.0% differ more than twofold upward in the rising arm before correction, 20.4% after correction; 21.3% differ more than twofold in the pulsatile direction after correction; 47.7% within twofold. Direction: CDKN1A +2.30, MDM2 +2.25, TP53I3 +2.60 corrected bits (arrest/terminal fate amplified by the sustained regime); DDB2 +0.18, XPC +0.45 (repair genes track exposure, shape-insensitive); PTEN null control flat. Verdict: confirmed with the registered correction; two decoder classes (counting-type and level-type) coexist.

### S19.5 P15: bulk-tissue projection boundary

Data: TCGA Pan-Cancer Atlas via cBioPortal (RPPA phospho-ERK = mean of MAPK1/3 PT202/Y204; total EGFR; MC3 mutation calls; COADREAD, LUAD, SKCM, PAAD; about 1,000 tumours). Verdict: falsified at bulk resolution. EGFR abundance does not predict phospho-ERK even in wild-type tumours (|r| <= 0.18 in all four types); a registered correction repaired a COADREAD/LUAD figure transposition in the first version of the verdict card (verdict direction unchanged). Registered asymmetry: this null does not overturn P13; it marks the resolution at which the chain becomes invisible. Side observation consistent with the chain: LUAD KRAS-G12 tumours carry higher mean phospho-ERK than wild type (+0.30 versus -0.07).

### S19.6 P16: frozen prospective prediction

KRAS-G13D, an allele of intermediate biochemical strength untested in any comparable public dataset, is predicted to show delta I in [-0.25, -0.05] bits versus KRASWT under the P13 pipeline, with ordering G12V > G13D > G12D. Frozen 2026-09-28; awaiting data.

### S19.7 Fig6 decay side-evidence (exploratory)

ERK deactivation half-times from the archived decay fits of ref. 78 (Fig. 6; EGF arm, r2 > 0.8 replicates): KrasWT 6.2 min versus G12D 4.7, G12V 2.8, Q61R 2.8, BRAF-V600E 5.3 min. Mutant lines deactivate faster, consistent with phosphatase/feedback renormalisation compensating the compressed channel. Six replicates per line; registered as directional side evidence only.
