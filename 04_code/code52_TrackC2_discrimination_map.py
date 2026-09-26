# -*- coding: utf-8 -*-
"""
code52 Track C2: real-data (alpha, nu) points on the topology discrimination map
Data sources (levels annotated faithfully):
  D1 in-house digitization: Batchelor2011 Fig1G/H (NCS amplitude 3 doses, UV amplitude 5 doses) -> alpha computed directly + bootstrap
  D2 original-paper statements (in-house fulltext.xml anchor): DSB pulse amplitude/duration are dose-independent, pulse count rises with dose
  D3 literature text values: Lahav 2004 (gamma-IR 0.1-10 Gy, N:1->~6); Moenke 2017 NCS nu~2 (via TCS-project analysis report)
Model anchors: NF = memo model-A table values; EXC = memo FHN + code 51 (nu=tau_r/T, parameter-dependent)
Outputs: (alpha,nu) discrimination map + tau_r inverted from the counting-law slope
"""
import sys, json
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))

WS = Path(r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911")
DATA = WS / "01_细胞线" / "公开数据" / "Batchelor2011" / "digitized_Fig1GH.json"
OUT_PNG = WS / "03_细胞线3" / "结果" / "代码52_TrackC2_alpha-nu判别图.png"
OUT_JSON = WS / "03_细胞线3" / "结果" / "代码52_TrackC2_结果.json"

RNG = np.random.default_rng(11)


def alpha_boot(dose, val, err, n=20000):
    """log-log slope bootstrap: each point resampled independently as N(val, err)."""
    dose = np.asarray(dose, float); val = np.asarray(val, float); err = np.asarray(err, float)
    slopes = []
    ld = np.log(dose)
    for _ in range(n):
        v = val + err * RNG.standard_normal(len(val))
        if np.any(v <= 0):
            continue
        slopes.append(np.polyfit(ld, np.log(v), 1)[0])
    slopes = np.array(slopes)
    return float(slopes.mean()), float(slopes.std(ddof=1)), float(np.percentile(slopes, 2.5)), float(np.percentile(slopes, 97.5))


def main():
    d = json.loads(DATA.read_text(encoding="utf-8"))
    a_ncs = alpha_boot(d["NCS_amp"]["dose"], d["NCS_amp"]["val"], d["NCS_amp"]["err"])
    a_uv = alpha_boot(d["UV_amp"]["dose"], d["UV_amp"]["val"], d["UV_amp"]["err"])
    print(f"α_NCS = {a_ncs[0]:.3f} ± {a_ncs[1]:.3f} (95%CI [{a_ncs[2]:.3f},{a_ncs[3]:.3f}])")
    print(f"α_UV  = {a_uv[0]:.3f} ± {a_uv[1]:.3f} (95%CI [{a_uv[2]:.3f},{a_uv[3]:.3f}])")

    # data points (alpha, nu, alpha_err, nu_err, label, source level)
    pts = [
        dict(name="γ-IR（DSB）", alpha=0.05, aerr=0.10, nu=1.4, nerr=0.4,
             src="D2+D3: Batchelor2011 original text 'amplitude dose-independent'; Lahav2004 N:1->~6/0.1-10Gy"),
        dict(name="NCS（DSB）", alpha=a_ncs[0], aerr=a_ncs[1], nu=2.0, nerr=0.5,
             src="D1: in-house digitized Fig1G; nu via Moenke2017/TCS report"),
        dict(name="UV", alpha=a_uv[0], aerr=a_uv[1], nu=0.0, nerr=0.2,
             src="D1: in-house digitized Fig1H; original text 'single pulse' -> nu~0"),
    ]
    # model anchors
    nf_alpha = np.log(9.47 / 1.24) / np.log(0.90 / 0.15)
    models = [
        dict(name="model anchor NF (memo model A)", alpha=float(nf_alpha), nu=0.0),
        dict(name="model anchor EXC (FHN, nu=tau_r/T)", alpha=0.04, nu=None),
    ]
    # counting-law slope inversion: real nu_gamma ~1.4 -> tau_r = nu*T ~ 1.4x5.5h
    tau_r_inv = 1.4 * 5.5
    print(f"counting-law inversion: tau_r ~ nu*T = 1.4x5.5h ~ {tau_r_inv:.1f} h (DSB-repair timescale magnitude)")

    out = dict(alpha_NCS=dict(zip(["mean", "sd", "ci_lo", "ci_hi"], a_ncs)),
               alpha_UV=dict(zip(["mean", "sd", "ci_lo", "ci_hi"], a_uv)),
               points=pts, models=models, tau_r_inversion_h=tau_r_inv)
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    # ---------- figure ----------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from daimon_runtime import setup_plot
    setup_plot()

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.8))

    # (a) (alpha, nu) discrimination plane
    ax = axes[0]
    ax.axvspan(-0.15, 0.1, color="#d5f5e3", alpha=0.7)
    ax.axvspan(0.5, 1.5, color="#fadbd8", alpha=0.7)
    ax.axvline(0.1, color="#27ae60", ls="--", lw=1)
    ax.axvline(0.5, color="#c0392b", ls="--", lw=1)
    ax.text(-0.02, 10.6, "EXC zone (digital signature)\nalpha<~0.1, nu>0", fontsize=9, color="#1e8449")
    ax.text(0.86, 10.6, "NF zone (analog signature)\nalpha>~0.5", fontsize=9, color="#922b21")
    colors = {"γ-IR（DSB）": "#2471a3", "NCS（DSB）": "#8e44ad", "UV": "#ca6f1e"}
    for p in pts:
        ax.errorbar(p["alpha"], p["nu"], xerr=p["aerr"], yerr=p["nerr"],
                    fmt="o", ms=9, capsize=5, color=colors[p["name"]],
                    label=f'{p["name"]}  α={p["alpha"]:.2f}±{p["aerr"]:.2f}, ν={p["nu"]:.1f}±{p["nerr"]:.1f}')
    ax.plot(nf_alpha, 0, "s", ms=11, color="#c0392b", mec="k",
            label=f"model anchor NF (memo table A) alpha={nf_alpha:.2f}, nu~0")
    ax.plot(0.04, 10.5, "^", ms=11, color="#1e8449", mec="k",
            label="model anchor EXC (code 51, nu=tau_r/T=10.5, parameter-dependent)")
    ax.annotate("", xy=(0.04, 10.5), xytext=(0.04, 0),
                arrowprops=dict(arrowstyle="-", color="#1e8449", ls=":", lw=1))
    ax.set_xlabel("alpha = dlnA/dlnD (amplitude-dose exponent)")
    ax.set_ylabel("nu = dN/dlnD (count-dose exponent)")
    ax.set_xlim(-0.15, 1.5); ax.set_ylim(-0.6, 11.5)
    ax.set_title("(a) Track C2: real data on the (alpha, nu) discrimination plane")
    ax.legend(fontsize=8, loc="center left")

    # (b) counting-law slope inversion
    ax = axes[1]
    lnD = np.linspace(-2.5, 2.5, 50)
    ax.plot(lnD, 1.4 * lnD + 3.0, "-", color="#2471a3", lw=2,
            label="real slope nu~1.4 (gamma-IR, Lahav)")
    ax.plot(lnD, 10.5 * lnD + 3.0, "--", color="#1e8449", lw=2,
            label="code-51 model slope nu=tau_r/T=10.5")
    ax.set_xlabel("ln D (log dose)"); ax.set_ylabel("N (pulse count, schematic shift)")
    ax.set_title("(b) counting-law slope = tau_r/T -> inverted repair timescale\n"
                 "tau_r ~ nu*T = 1.4x5.5h ~ 7.7h (~DSB repair timescale)")
    ax.legend(fontsize=9)
    ax.text(0.02, 0.95, "Lemma 3: N ~ (tau_r/T)*ln(D0/D_c)\nthe slope itself is an identifiable estimator of the repair time constant",
            transform=ax.transAxes, fontsize=9, va="top",
            bbox=dict(boxstyle="round", fc="#fef9e7", ec="#b7950b"))

    fig.suptitle("code 52 · Track C2: real (alpha, nu) signatures vs topology-selection-theorem discrimination zones", y=0.98)
    fig.tight_layout()
    fig.savefig(OUT_PNG, bbox_inches="tight")
    plt.close(fig)
    print("written:", OUT_PNG.name, OUT_JSON.name)


if __name__ == "__main__":
    main()
