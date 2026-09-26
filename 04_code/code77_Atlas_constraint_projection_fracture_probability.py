# -*- coding: utf-8 -*-
"""
Code 77 — Atlas constraint projection and fracture probability — BpC
Version v1.0.0 · 2026-08-16

First step "from auditor to architect" (plan revised to match data reality):
The GPCR Atlas has no raw dose-response curves, so the original "constraint-embedding fit" plan has no raw material;
revised into two constructive operators that need no raw data:

  B′ | Constraint projection (minimum-cost repair):
     The within-group identity requires q_i ≡ logRA_i − tc_i to be equal across the group (= −log E_sys realized value).
     Projection = moving each entry's (tc, logRA) in least squares onto the identity manifold q_i = q̄:
       under equal-split repair each parameter moves (q_i−q̄)/2 dex; group total cost = Σ(q_i−q̄)²/2.
     Outputs: per-group repair cost, per-parameter move, largest-gap ligand ("who pays the most for consistency").
     Position statement: the projection output is "the most economical explanation assuming the framework is true", not "repaired ground truth".

  C | Fracture-probability mixture model:
     dev_entry modeled as a two-component zero-mean normal mixture (healthy N(0,σ_h) + fractured N(0,σ_b)), EM estimation,
     threshold-free per-entry P(fracture). σ_h is also a measured scale of "healthy literature noise".

  Cross-validation: Spearman ρ of group-level repair cost vs within-group mean P(fracture).

Input: atlas_audit_entries.csv (product of code 71)
Output: atlas_repair_projection.csv (group level), atlas_fracture_probability.csv (entry level)

Rulings (run of 2026-08-16):
  B′: 24.2% of 269 groups need >0.3 dex-level repair — converging with code 71 (26%) and code 73 (28.2% conservative convention);
      three methods agree.
  C: EM converged healthy fraction 74.5% (σ_h=0.089 dex, σ_b=0.68 dex) — threshold-free independent confirmation
      of the "about 3/4 healthy" conclusion; P>0.9 entries 15.0% ≈ the 14.5% with |dev|>0.3.
      σ_h=0.089 dex measured the reporting-noise scale of healthy literature (larger than the legitimate mechanism spread
      0.004–0.05 of code 72 M2; the difference is reporting/digitization noise).
  Cross: ρ=0.962, the two operators strongly corroborate each other.
"""
import numpy as np, pandas as pd

ENT = '/mnt/agents/output/04_细胞线4/结果/atlas_audit/atlas_audit_entries.csv'
OUT = '/mnt/agents/output/04_细胞线4/结果/atlas_audit'

def main():
    E = pd.read_csv(ENT).dropna(subset=['dev_entry'])
    E['grp'] = (E['doi'].astype(str)+'|'+E['receptor'].astype(str)+'|'+
                E['Measured process'].astype(str)+'|'+E['Pathway level'].astype(str)+'|'+
                E['Cell line'].astype(str)+'|'+E['Primary effector subtype'].astype(str))
    E['q'] = E['logRA'] - E['tc']

    # ---- B′ constraint projection ----
    recs = []
    for gname, d in E.groupby('grp'):
        q = d['q'].to_numpy(); n = len(q)
        if n < 2: continue
        qstar = q.mean(); move = q - qstar
        i_w = np.argmax(np.abs(move))
        recs.append({'grp':gname,'n':n,'cost_total_dex2':(move**2).sum()/2,
                     'cost_rms_dex':np.sqrt(((move/2)**2).mean()),
                     'worst_ligand':d['ligand'].iloc[i_w],'worst_move_dex':move[i_w]})
    BJ = pd.DataFrame(recs)
    E['repair_move'] = E['q'] - E['grp'].map(E.groupby('grp')['q'].mean())
    E['repair_per_param'] = E['repair_move']/2

    # ---- C mixture-model EM ----
    x = E['dev_entry'].to_numpy()
    pi, s_h, s_b = 0.75, 0.05, 0.8
    for it in range(2000):
        lh = pi*np.exp(-x**2/(2*s_h**2))/s_h
        lb = (1-pi)*np.exp(-x**2/(2*s_b**2))/s_b
        pb = lb/(lh+lb)
        pi_n = 1-pb.mean()
        sh_n = np.sqrt(((1-pb)*x**2).sum()/(1-pb).sum())
        sb_n = np.sqrt((pb*x**2).sum()/pb.sum())
        if max(abs(pi_n-pi), abs(sh_n-s_h), abs(sb_n-s_b)) < 1e-10:
            pi,s_h,s_b = pi_n,sh_n,sb_n; break
        pi,s_h,s_b = pi_n,sh_n,sb_n
    E['p_broken'] = pb

    BJ.to_csv(f'{OUT}/atlas_repair_projection.csv', index=False)
    E.to_csv(f'{OUT}/atlas_fracture_probability.csv', index=False)

    from scipy.stats import spearmanr
    chk = E.groupby('grp').agg(pmean=('p_broken','mean')).reset_index().merge(BJ, on='grp')
    rho = spearmanr(chk.cost_rms_dex, chk.pmean).statistic

    print(f"[B′] groups {len(BJ)}; median repair-cost RMS {BJ['cost_rms_dex'].median():.3f} dex; "
          f"groups needing >0.3 dex-level repair {(BJ['cost_rms_dex']>0.15).mean():.1%}")
    print(f"[C] EM converged @{it} steps: pi_healthy={pi:.3f}, σ_h={s_h:.4f}, σ_b={s_b:.3f} dex; "
          f"P>0.9 entries {(pb>0.9).mean():.1%} (cf. |dev|>0.3: {(np.abs(x)>0.3).mean():.1%})")
    print(f"[cross] repair cost vs fracture probability Spearman ρ = {rho:.3f}")
    print("[B′ top-5 costliest groups]")
    print(BJ.nlargest(5,'cost_rms_dex')[['grp','n','cost_rms_dex','worst_ligand','worst_move_dex']].to_string(index=False))

if __name__ == '__main__':
    main()
