# -*- coding: utf-8 -*-
"""
code80_Gillis_arm_coordinate_extended_model.py
version v1.0.0 · 2026-08-16

Purpose: test the minimal viable embryo of a "next-generation model" — not inventing new equations from
      scratch, but adding one constrained assay-arm coordinate delta_arm onto the old skeleton to see whether the Gillis muOR rupture is legitimately absorbed.

Model:
  old model M0: latent log-tau and pKA per ligand x arm, back-computed from Emax and pEC50;
             main-text Table 3 log-tau and SI Table S1 log(tau/KA) should scatter zero-mean around that latent coordinate.
  extended model M1: beyond the same latent coordinate, each assay arm may carry one shared reporter/functional coordinate delta_arm:
             log-tau_obs = log-tau_latent + delta_arm + eps
             logR_obs = log-tau_latent + pKA_latent + delta_arm + eps
             delta_arm ~ N(0, sigma_delta^2)
  compare M0/M1 by marginal likelihood; sigma_delta is both sensitivity-scanned and empirical-Bayes estimated.

Data-handling discipline:
  the SI logR columns of GPA and cAMP are identical in the original paper (double-entered by code 66); the primary
  analysis drops cAMP-logR to avoid counting duplicated table values as two independent evidences; a no-drop sensitivity analysis is also run.

Primary-analysis results (2026-08-16, confirmed by an independent re-run):
- primary analysis: 66 residuals; empirical-Bayes sigma_delta = 0.374 dex.
- M0 log marginal likelihood = -121.638; M1 = +6.589; delta logML = +128.227.
- GIRK posterior coordinate: delta_GIRK = -0.805 ± 0.049 dex, 95% interval [-0.902, -0.708].
- after adding delta_arm, the GIRK median residual |r| drops from 0.811 dex to 0.107 dex, max 0.249 dex.
- leave-one-ligand-out keeps delta_GIRK within [-0.843, -0.732]; not driven by a single ligand.
- secondary signal: cAMP log-tau coordinate delta=+0.508 ± 0.163 dex; its SI logR was dropped for the GPA/cAMP duplication,
  so the evidence is weaker than GIRK and is flagged as to-be-rechecked rather than a main conclusion.
- sensitivity (keeping duplicate cAMP logR): sigma_delta=0.320, delta logML=+123.844, delta_GIRK=-0.800 ± 0.049.

Verdict:
  the minimal embryo of the next-generation model holds, but its identity must be stated clearly — it is a
  "measurement/coordinate-aware" extended model, not a new receptor-biophysics theory. It proves that the Gillis
  GIRK rupture can be systematically absorbed by an arm-level coordinate term of about -0.8 dex; whether that term
  corresponds to real coupling structure, data-processing shift, or a table-production accident, current data cannot adjudicate.

Outputs:
- /mnt/agents/output/04_细胞线4/结果/code80_arm_coordinate_residuals.csv
- /mnt/agents/output/04_细胞线4/结果/code80_arm_coordinate_posterior.csv
- /mnt/agents/output/04_细胞线4/结果/code80_model_comparison.csv
- /mnt/agents/output/04_细胞线4/结果/code80_leave_one_out_girk.csv
"""

import os

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import norm

OUTDIR = "/mnt/agents/output/04_细胞线4/结果"
LIG6 = ["Morphine", "Oxycodone", "Oliceridine", "PZM21", "SR-17018", "Buprenorphine"]
LN10 = np.log(10)

# Gillis 2020 main-text/SI table values, checked value by value against code 66
EMAX = {
    "Nb33": [71, 3, 70, 4, 42, 8, 38, 3, 20, 6, 26, 3],
    "mGsi": [80, 7, 75, 3, 51, 7, 52, 5, 35, 5, 36, 3],
    "GPA": [98, 5, 103, 4, 85, 5, 86, 7, 61, 13, 79, 6],
    "cAMP": [97, 4, 106, 7, 86, 5, 84, 7, 62, 12, 86, 9],
    "GIRK": [88, 4, None, None, 74, 4, 86, 3, 78, 7, 54, 5],
    "GRK2rec": [66, 6, 64, 8, 40, 3, 36, 7, 41, 9, 28, 5],
    "BarrGRK2": [70, 3, 72, 5, 58, 3, 59, 4, 49, 7, 34, 2],
}
PEC50 = {
    "Nb33": [6.67, .03, 5.93, .12, 7.19, .28, 6.88, .36, 7.48, .46, 8.28, .21],
    "mGsi": [6.94, .05, 6.22, .11, 7.46, .18, 7.49, .14, 7.15, .33, 8.65, .14],
    "GPA": [7.72, .16, 6.94, .23, 8.38, .19, 8.16, .19, 7.66, .18, 8.85, .06],
    "cAMP": [8.10, .26, 7.01, .18, 8.66, .13, 8.64, .25, 7.67, .33, 9.61, .37],
    "GIRK": [7.98, .10, None, None, 8.43, .12, 8.42, .07, 6.43, .13, 7.97, .21],
    "GRK2rec": [6.78, .08, 6.28, .09, 7.26, .51, 7.58, .19, 6.99, .19, 7.95, .28],
    "BarrGRK2": [7.31, .13, 6.22, .06, 7.71, .08, 7.56, .07, 6.48, .46, 8.50, .21],
}
LOGTAU = {
    "Nb33": [0.42, .08, 0.56, .12, -0.34, .12, -0.33, .07, -0.86, .21, -0.62, .09],
    "mGsi": [0.70, .07, 0.65, .09, -0.12, .05, -0.08, .08, -0.37, .16, -0.40, .05],
    "GPA": [1.74, .22, 1.86, .26, 1.15, .25, 1.18, .28, 0.66, .37, 0.63, .25],
    "cAMP": [2.00, .31, 1.87, .22, 1.29, .23, 1.44, .39, 1.04, .28, 1.35, .39],
    "GIRK": [0.09, .05, None, None, -0.24, .05, -0.18, .04, -0.28, .12, -0.61, .10],
    "GRK2rec": [0.21, .07, 0.22, .14, -0.30, .18, -0.37, .11, -0.30, .09, -0.57, .08],
    "BarrGRK2": [0.34, .07, 0.35, .09, 0.13, .03, 0.11, .06, -0.04, .15, -0.30, .09],
}
LOGR = {
    "Nb33": [6.5, .04, 5.8, .11, 6.7, .17, 6.9, .18, 7.1, .37, 7.6, .2],
    "mGsi": [6.8, .06, 6.1, .1, 7.4, .08, 7.2, .1, 6.6, .38, 8.0, .14],
    "GPA": [7.8, .19, 7.0, .24, 8.2, .24, 8.0, .27, 6.1, .58, 8.6, .13],
    "cAMP": [7.8, .19, 7.0, .24, 8.2, .24, 8.0, .27, 6.1, .58, 8.6, .13],
    "GIRK": [6.9, .03, None, None, 7.5, .04, 7.3, .03, 5.5, .07, 7.0, .06],
    "BarrGRK2": [7.2, .14, 6.1, .07, 7.5, .07, 7.3, .07, 6.7, .5, None, None],
}


def logtau_from_emax(em, eme):
    tau = em / (100 - em)
    return np.log10(tau), (1 / LN10) * 100 * eme / (em * (100 - em))


def s_log1ptau(logtau, sd_logtau):
    tau = 10 ** logtau
    return (tau / (1 + tau)) * sd_logtau


def build_residuals(exclude_duplicate=True):
    rows = []
    for arm in EMAX:
        for j, ligand in enumerate(LIG6):
            em, eme = EMAX[arm][2*j:2*j+2]
            pe, pee = PEC50[arm][2*j:2*j+2]
            if em is None or pe is None or em >= 100:
                continue

            lt_latent, sd_lt = logtau_from_emax(em, eme)
            pka = pe - np.log10(1 + 10 ** lt_latent)
            sd_pka = np.sqrt(pee ** 2 + s_log1ptau(lt_latent, sd_lt) ** 2)

            if arm in LOGTAU:
                lt_obs, lt_se = LOGTAU[arm][2*j:2*j+2]
                if lt_obs is not None:
                    rows.append(dict(
                        arm=arm, ligand=ligand, channel="logtau",
                        residual=lt_obs - lt_latent,
                        sd=np.sqrt(lt_se ** 2 + sd_lt ** 2),
                    ))

            if arm in LOGR:
                # GPA/cAMP SI logR columns are exact duplicates in the original: primary analysis drops cAMP-logR
                if exclude_duplicate and arm == "cAMP":
                    continue
                lr_obs, lr_se = LOGR[arm][2*j:2*j+2]
                if lr_obs is not None:
                    pred = lt_latent + pka
                    sd = np.sqrt(lr_se ** 2 + sd_lt ** 2 + sd_pka ** 2)
                    rows.append(dict(
                        arm=arm, ligand=ligand, channel="logR",
                        residual=lr_obs - pred, sd=sd,
                    ))

    return pd.DataFrame(rows)


def loglik_m0(df):
    return float(norm.logpdf(df.residual, 0, df.sd).sum())


def loglik_m1(df, sigma_delta):
    total = 0.0
    for _, group in df.groupby("arm"):
        r = group.residual.values
        sd = group.sd.values
        n = len(r)
        cov = np.diag(sd ** 2) + sigma_delta ** 2 * np.ones((n, n))
        sign, logdet = np.linalg.slogdet(cov)
        total += -0.5 * (r @ np.linalg.solve(cov, r) + logdet + n * np.log(2 * np.pi))
    return float(total)


def posterior_by_arm(df, sigma_delta):
    rows = []
    for arm, group in df.groupby("arm"):
        weights = 1 / group.sd ** 2
        precision = 1 / sigma_delta ** 2 + weights.sum()
        mean = float((group.residual * weights).sum() / precision)
        sd = float(np.sqrt(1 / precision))
        rows.append(dict(
            arm=arm,
            n=len(group),
            delta=mean,
            delta_sd=sd,
            ci95_lo=mean - 1.96 * sd,
            ci95_hi=mean + 1.96 * sd,
        ))
    return pd.DataFrame(rows)


def add_m1_residuals(df, posterior):
    delta = dict(zip(posterior.arm, posterior.delta))
    out = df.copy()
    out["residual_m1"] = out.apply(lambda row: row.residual - delta[row.arm], axis=1)
    out["delta_arm"] = out.arm.map(delta)
    return out


def leave_one_out_girk(df, sigma_delta):
    rows = []
    for ligand in LIG6:
        group = df[(df.arm == "GIRK") & (df.ligand != ligand)]
        weights = 1 / group.sd ** 2
        precision = 1 / sigma_delta ** 2 + weights.sum()
        mean = float((group.residual * weights).sum() / precision)
        sd = float(np.sqrt(1 / precision))
        rows.append(dict(leave_out=ligand, delta=mean, delta_sd=sd))
    return pd.DataFrame(rows)


def fit_model(df, label):
    ll0 = loglik_m0(df)
    opt = minimize_scalar(
        lambda s: -loglik_m1(df, s), bounds=(0.01, 2.0), method="bounded"
    )
    sigma_hat = float(opt.x)
    ll1 = float(-opt.fun)
    posterior = posterior_by_arm(df, sigma_hat)
    residuals = add_m1_residuals(df, posterior)
    comparison = pd.DataFrame([
        dict(model="M0_classic", sigma_delta=0.0, log_marginal_likelihood=ll0),
        dict(model="M1_arm_coordinate", sigma_delta=sigma_hat,
             log_marginal_likelihood=ll1),
    ])
    comparison["delta_logml_vs_m0"] = comparison.log_marginal_likelihood - ll0
    comparison["dataset"] = label
    return residuals, posterior, comparison, sigma_hat, ll0, ll1


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    primary = build_residuals(exclude_duplicate=True)
    sensitivity = build_residuals(exclude_duplicate=False)

    res, post, comp, sigma_hat, ll0, ll1 = fit_model(primary, "primary_exclude_cAMP_logR_duplicate")
    res_s, post_s, comp_s, sigma_s, ll0_s, ll1_s = fit_model(sensitivity, "sensitivity_keep_duplicate")

    loo = leave_one_out_girk(primary, sigma_hat)

    # fixed-prior sensitivity: 0.089 comes from the code-77 healthy-literature noise scale
    fixed_rows = []
    for sigma in [0.05, 0.089, 0.2, 0.5, 1.0]:
        fixed_rows.append(dict(
            sigma_delta=sigma,
            logml_m1=loglik_m1(primary, sigma),
            delta_logml_vs_m0=loglik_m1(primary, sigma) - ll0,
        ))
    fixed = pd.DataFrame(fixed_rows)

    res.to_csv(os.path.join(OUTDIR, "code80_arm_coordinate_residuals.csv"), index=False)
    post.to_csv(os.path.join(OUTDIR, "code80_arm_coordinate_posterior.csv"), index=False)
    pd.concat([comp, comp_s], ignore_index=True).to_csv(
        os.path.join(OUTDIR, "code80_model_comparison.csv"), index=False
    )
    fixed.to_csv(os.path.join(OUTDIR, "code80_sigma_delta_sensitivity.csv"), index=False)
    loo.to_csv(os.path.join(OUTDIR, "code80_leave_one_out_girk.csv"), index=False)

    girk = post[post.arm == "GIRK"].iloc[0]
    girk_res = res[res.arm == "GIRK"]
    print("code80 Gillis assay-arm coordinate extended model done")
    print(f"primary N={len(primary)}, sigma_delta={sigma_hat:.3f} dex")
    print(f"M0 logML={ll0:.3f}, M1 logML={ll1:.3f}, delta={ll1-ll0:.3f}")
    print(f"GIRK delta={girk.delta:.3f} ± {girk.delta_sd:.3f} dex, 95%CI=({girk.ci95_lo:.3f},{girk.ci95_hi:.3f})")
    print(f"GIRK |r| median: {girk_res.residual.abs().median():.3f} -> {girk_res.residual_m1.abs().median():.3f} dex; max={girk_res.residual_m1.abs().max():.3f}")
    print(f"sensitivity (duplicates kept): sigma_delta={sigma_s:.3f}, delta={ll1_s-ll0_s:.3f}")


if __name__ == "__main__":
    main()
